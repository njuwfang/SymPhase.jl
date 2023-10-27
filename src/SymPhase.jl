"""
A high-performance module for sampling Stabilizer circuits by symbolizing the phases of Pauli strings. 
"""
module SymPhase

using PrecompileTools

using LoopVectorization
using SparseArrays
using BitIntegers

export SymStabilizer, cX, cY, cZ, cPhase, cHadamard, cCNOT
export all_zeros, apply!, mul_left!, isX, isY, isZ, zero!, rowcopy!, rowswap!, parse_stim, Sampler
export transpose_p!, transpose_d!, _transpose_16x_!, _transpose__x16!, projectZ!, _transposed_index, _isone
#export one, zero, copy

const _num_packed_bits_ = 512
const _shift3_ = Int(log2(_num_packed_bits_))

@inline _div4(j) = (j-1)>>9+1
@inline _div3(j) = (j-1)>>_shift3_+1
@inline _div2(j) = ((j-1)&511)>>3+1
@inline _div1(j) = (j-1)&(_num_packed_bits_-1)+1
@inline _shift(j) = ((j-1)&511)&7
@inline _pow(j) = one(UInt8)<<_shift(j)
@inline _div8(j) = (j-1)>>3+1
@inline _div32(j) = (j-1)>>5+1
@inline _div64(j) = (j-1)>>6+1
@inline _div128(j) = (j-1)>>7+1
@inline _pow32(j) = one(UInt32)<<((j-1)&31)
@inline _pow64(j) = one(UInt64)<<((j-1)&63)

struct SymStabilizer
    nq::Int
    ns::Int
    len3::Int
    len4::Int
    enable_T::Bool
    max_ns::Array{Int64,2}
    min_ns::Array{Int64,2}
    xzs::Array{UInt8,4}
    phases::Array{UInt8,2}
    temp_phases::Array{Bool,2}
    symbols::Union{Array{UInt32,2},Array{UInt32,3}}
    T::Array{UInt32,3}
    T_inv::Array{UInt32,3}
    temp::Array{UInt8,2}
    pow32::Array{UInt32,1}
    temp_s::Array{Bool,1}
end

function Base.zero(::Type{SymStabilizer}, nq, ns;enable_T=false)
    len3 = Int(ceil((nq+1)/_num_packed_bits_))
    len4 = len3
    lens = Int(ceil(ns/32))
    if enable_T
        T = zeros(UInt32, (_num_packed_bits_*len3)>>4, _num_packed_bits_, len3<<1)
        T_inv = zeros(UInt32, (_num_packed_bits_*len3)>>4, _num_packed_bits_, len3<<1)
        @inbounds for row in 1:nq
            j3 = _div3(row)
            j1 = _div1(row)
            d = _div32(row)
            pow = _pow32(row)
            T[d,j1,j3] = pow
            T_inv[d,j1,j3] = pow
        end
        @inbounds for row in 1+len3<<_shift3_:nq+len3<<_shift3_
            j3 = _div3(row)
            j1 = _div1(row)
            d = _div32(row)
            pow = _pow32(row)
            T[d,j1,j3] = pow
            T_inv[d,j1,j3] = pow
        end
        symbols = zeros(UInt32, (_num_packed_bits_*len3)>>4, ns)
    else
        T = zeros(UInt32, 0,0,0)
        T_inv = zeros(UInt32, 0,0,0)
        symbols = zeros(UInt32, lens, _num_packed_bits_, len3<<1)
    end
    SymStabilizer(
        nq, ns, len3, len4, enable_T,
        zeros(Int64,_num_packed_bits_,len3<<1),
        fill(typemax(Int64), (_num_packed_bits_,len3<<1)),
        zeros(UInt8, _num_packed_bits_, 64, len3<<1, len4<<1),
        zeros(UInt8, _num_packed_bits_, len3<<1),
        zeros(Bool, _num_packed_bits_, len3<<1),
        symbols,
        T,
        T_inv,
        zeros(UInt8, _num_packed_bits_, 64),
        zeros(UInt32, 32),
        zeros(Bool, ns)
    )
end

function all_zeros(::Type{SymStabilizer}, nq, ns;enable_T=false)
    q = zero(SymStabilizer, nq, ns;enable_T=enable_T)
    id64 = 0x8040201008040201
    xzs64 = reinterpret(UInt64, q.xzs)
    ss = size(xzs64, 1)÷size(xzs64, 2)
    for j4 in axes(xzs64, 4)
        for j2 in axes(xzs64, 2)
            xzs64[((j4-1)%ss)*size(xzs64, 2)+j2,j2,(j4-1)÷ss+1,j4] = id64
        end
    end
    q
end

# index for measurement
@inline function Base.getindex(q::SymStabilizer, row, col)
    dx4 = _div4(col)
    dz4 = dx4+q.len4
    dr3 = _div3(row)
    dx2 = _div2(col)
    dr1 = _div1(row)
    pow = _pow(col)

    n1,n2,_,_ = size(q.xzs)
    j = (dr1-1)*n2+dx2-1
    j1 = j%n1+1
    j2 = j÷n1+1

    q.xzs[j1,j2,dr3,dx4]&pow!=0, q.xzs[j1,j2,dr3,dz4]&pow!=0
end

# assign for measurement
@inline function Base.setindex!(q::SymStabilizer, (x, z)::Tuple{Bool, Bool}, row, col)
    dx4 = _div4(col)
    dz4 = dx4+q.len4
    dr3 = _div3(row)
    dx2 = _div2(col)
    dr1 = _div1(row)
    pow = _pow(col)

    n1,n2,_,_ = size(q.xzs)
    j = (dr1-1)*n2+dx2-1
    j1 = j%n1+1
    j2 = j÷n1+1

    if x
        q.xzs[j1,j2,dr3,dx4] |= pow
    else
        q.xzs[j1,j2,dr3,dx4] &= ~pow
    end

    if z
        q.xzs[j1,j2,dr3,dz4] |= pow
    else
        q.xzs[j1,j2,dr3,dz4] &= ~pow
    end

    nothing
end

# set a row to I for measurement
function zero!(q::SymStabilizer, row;zero_T=false)
    n1,_,n3,n4 = size(q.xzs)
    xzs = reshape(reinterpret(UInt512, q.xzs), n1,n3,n4)

    dr3 = _div3(row)
    dr1 = _div1(row)

    @inbounds for j4 in axes(xzs, 3)
        xzs[dr1,dr3,j4] = zero(UInt512)
    end

    q.phases[dr1,dr3] = zero(UInt8)
    
    if q.enable_T
        if zero_T
            @turbo for j1 in axes(q.T, 1)
                q.T[j1,dr1,dr3] = zero(UInt32)
            end
        else
            # setting symbols to 0 is too time-consuming 
            low = q.min_ns[dr1,dr3]
            high = q.max_ns[dr1,dr3]
            @inbounds for j in low:high
                cnt = zero(UInt32)
                @turbo for k in axes(q.symbols, 1)
                    cnt += count_ones(q.T[k,dr1,dr3] & q.symbols[k,j])
                end
                if cnt&1!=0
                    @turbo for k1 in axes(q.symbols, 1)
                        q.symbols[k1,j] ⊻= q.T_inv[k1, dr1, dr3]
                    end
                end
            end
        end
    else
        l = _div32(q.min_ns[dr1,dr3])
        h = _div32(q.max_ns[dr1,dr3])
        @turbo for j in l:h
            q.symbols[j,dr1,dr3] = zero(UInt32)
        end 
    end

    q.min_ns[dr1,dr3] = typemax(Int64)
    q.max_ns[dr1,dr3] = 0
    
    nothing
end

# swap two rows for measurement
function rowswap!(q::SymStabilizer, t, s)
    n1,_,n3,n4 = size(q.xzs)
    xzs = reshape(reinterpret(UInt512, q.xzs), n1,n3,n4)

    dt3 = _div3(t)
    dt1 = _div1(t)
    ds3 = _div3(s)
    ds1 = _div1(s)

    @inbounds for j4 in axes(xzs, 3)
        xzs[dt1,dt3,j4], xzs[ds1,ds3,j4] = xzs[ds1,ds3,j4], xzs[dt1,dt3,j4]
    end

    q.phases[dt1,dt3], q.phases[ds1,ds3] = q.phases[ds1,ds3], q.phases[dt1,dt3]

    if q.enable_T
        @turbo for j1 in axes(q.T, 1)
            a1, b1 = q.T[j1,dt1,dt3], q.T[j1,ds1,ds3]
            q.T[j1,ds1,ds3] = a1
            q.T[j1,dt1,dt3] = b1
            a2, b2 = q.T_inv[j1,dt1,dt3], q.T_inv[j1,ds1,ds3]
            q.T_inv[j1,ds1,ds3] = a2
            q.T_inv[j1,dt1,dt3] = b2
        end
    else
        l = _div32(min(q.min_ns[dt1,dt3], q.min_ns[ds1,ds3]))
        h = _div32(max(q.max_ns[dt1,dt3], q.max_ns[ds1,ds3]))
        @turbo for j in l:h#axes(q.symbols, 1)
            a, b = q.symbols[j,ds1,ds3], q.symbols[j,dt1,dt3]
            q.symbols[j,ds1,ds3] = b
            q.symbols[j,dt1,dt3] = a
        end
    end

    q.max_ns[dt1,dt3], q.max_ns[ds1,ds3] = q.max_ns[ds1,ds3], q.max_ns[dt1,dt3]
    q.min_ns[dt1,dt3], q.min_ns[ds1,ds3] = q.min_ns[ds1,ds3], q.min_ns[dt1,dt3]

    nothing
end

function _isone(q::SymStabilizer, i1, i3, l)
    if q.enable_T
        cnt = zero(UInt32)
        @inbounds @simd for j1 in axes(q.symbols, 1)
            cnt ⊻= q.T[j1,i1,i3] & q.symbols[j1,l]
        end
        return count_ones(cnt)&1!=0
    else
        d = _div32(l)
        pow = _pow32(l)
        return q.symbols[d,i1,i3]&pow!=0
    end
end

include("transpose.jl")

include("clifford_gates.jl")

include("mul_left.jl")

include("measurement.jl")

include("circuit.jl")

include("parse_stim.jl")

include("bit_array.jl")

include("sampler.jl")

include("misc.jl")


@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    initfile = joinpath(@__DIR__, "../test/precompile.stim")
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        circuit = parse_stim(initfile)
        sampler = Sampler(circuit)
    end
end

end
