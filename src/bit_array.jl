using Random: rand!, randexp

mutable struct OpSymbol
    len::Int
    nops::Int
    ops::Vector{GateType}
    args::Vector{Float32}
end

OpSymbol() = OpSymbol(0,0,Vector{GateType}(undef, 128),Vector{Float32}(undef, 128))

@eval function Base.push!(s::OpSymbol, op::GateType, arg::Float32=0.0f0)
    if s.nops == size(s.ops, 1)
        s.ops = _double_(s.ops)
        s.args = _double_(s.args)
    end
    s.nops += 1
    s.ops[s.nops] = op
    s.args[s.nops] = arg
    if $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_M_], " || ")))
        s.len += 1
    elseif $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_E_], " || ")))
        s.len += 1
    elseif op == DEPOLARIZE1
        s.len += 2
    elseif op == DEPOLARIZE2
        s.len += 4
    end

    s
end

@eval function biased_random_bits11(s::OpSymbol, shots)
    C = BitMatrix(undef, shots, s.len)
    temp = BitVector(undef, shots)
    temp1 = BitMatrix(undef, shots, 2)
    temp2 = BitMatrix(undef, shots, 4)
    n = 0
    for j in 1:s.nops
        op = s.ops[j]
        if $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_M_], " || ")))
            n += 1
            rand!(temp)
            C[:,n] = temp
        elseif $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_E_], " || ")))
            n += 1
            set_biased_random_bits!(s.args[j], temp)
            C[:,n] = temp
        elseif op == DEPOLARIZE1
            n += 2
            set_biased_random_bits!(s.args[j], temp1)
            C[:,n-1:n] = temp1
        elseif op == DEPOLARIZE2
            n += 4
            set_biased_random_bits!(s.args[j], temp2)
            C[:,n-3:n] = temp2
        end
    end
    C
end

@eval function biased_random_bits(s::OpSymbol, shots)
    lens = _div64(shots)
    C = Matrix{UInt64}(undef, lens, s.len)
    temp = Vector{UInt64}(undef, lens)
    temp1 = Matrix{UInt64}(undef, lens, 2)
    temp2 = Matrix{UInt64}(undef, lens, 4)
    n = 0
    for j in 1:s.nops
        op = s.ops[j]
        if $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_M_], " || ")))
            n += 1
            rand!(temp)
            C[:,n] = temp
        elseif $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_E_], " || ")))
            n += 1
            set_biased_random_bits!(s.args[j], temp)
            C[:,n] = temp
        elseif op == DEPOLARIZE1
            n += 2
            set_biased_random_bits!(s.args[j], temp1)
            C[:,n-1:n] = temp1
        elseif op == DEPOLARIZE2
            n += 4
            set_biased_random_bits!(s.args[j], temp2)
            C[:,n-3:n] = temp2
        end
    end
    C
end

"""
Generate random BitArray with biased probability.

See Stim, https://algassert.com/post/2200
"""
function biased_random_bits(p, N)
    (p >= 0.0 && p <= 1.0) || error("Invalid probability: $p")
    C = BitVector(undef, N)
    biased_random_bits!(p, C)
end

function biased_random_bits!(p, C::BitVector)
    if p > 0.5
        biased_random_bits!(1-p, C)
        Cc = C.chunks
        @inbounds @simd for j in axes(Cc, 1)
            Cc[j] ⊻= typemax(UInt64)
        end
    elseif p == 0.5
        rand!(C)
    elseif p < 0.02#0390625 # 1/2^8
        set_biased_random_bits!(p, C)
    else
        Cc = C.chunks
        COIN_FLIPS = 8
        BUCKETS = float(1 << COIN_FLIPS)
        raised = p * BUCKETS
        raised_floor = floor(raised)
        p_left = (raised - raised_floor) / (BUCKETS - raised_floor)
        p_top_bits = UInt64(raised_floor)

        for j in eachindex(Cc)
            result = zero(UInt64)
            alive = rand(UInt64)
            for k_bit = (COIN_FLIPS-2):-1:0
                shoot = rand(UInt64)
                result ⊻= shoot & alive & -((p_top_bits>>k_bit) & 1)
                alive &= ~shoot
            end
            Cc[j] = result
        end

        or_biased_random_bits!(p_left, C)
    end
end

function set_biased_random_bits!(p, C::BitVector)
    Cc = C.chunks
    @inbounds @simd for j in axes(Cc, 1)
        Cc[j] = zero(UInt64)
    end
    if p == 0.0
        return
    end
    N = length(C)
    k = geo(p) + 1
    while k <= N
        C[k] = true
        k += geo(p) + 1
    end
end

function set_biased_random_bits!(p, C::Vector{UInt64})
    @inbounds @simd for j in axes(C, 1)
        C[j] = zero(UInt64)
    end
    if p == 0.0
        return
    end
    N = size(C,1)<<6
    k = geo(p) + 1
    while k <= N
        C[_div64(k)] ⊻= _pow64(k)
        k += geo(p) + 1
    end
end

function or_biased_random_bits!(p, C::BitVector)
    if p == 0.0
        return
    end
    N = length(C)
    k = geo(p) + 1
    while k <= N
        C[k] |= true
        k += geo(p) + 1
    end
end

function set_biased_random_bits!(p, C::BitMatrix)
    Cc = C.chunks
    @inbounds @simd for j in axes(Cc, 1)
        Cc[j] = zero(UInt64)
    end
    if p == 0.0
        return
    end
    N, M = size(C)
    k = geo(p) + 1
    if M == 2    
        while k <= N
            r = 1+rand(UInt8)%3
            C[k,1] = r&1!=0
            C[k,2] = r&2!=0
            k += geo(p) + 1
        end
    elseif M == 4
        while k <= N
            r = 1+rand(UInt8)%15
            C[k,1] = r&1!=0
            C[k,2] = r&2!=0
            C[k,3] = r&4!=0
            C[k,4] = r&8!=0
            k += geo(p) + 1
        end
    end
    return
end

function set_biased_random_bits!(p, C::Matrix{UInt64})
    @inbounds @simd for j in eachindex(C)
        C[j] = zero(UInt64)
    end
    if p == 0.0
        return
    end
    N, M = size(C)
    N = N<<6
    k = geo(p) + 1
    if M == 2    
        while k <= N
            r = 1+rand(UInt8)%3
            dk = _div64(k)
            pk = _pow64(k)
            C[dk,1] ⊻= (r&1!=0)*pk
            C[dk,2] ⊻= (r&2!=0)*pk
            k += geo(p) + 1
        end
    elseif M == 4
        while k <= N
            r = 1+rand(UInt8)%15
            dk = _div64(k)
            pk = _pow64(k)
            C[dk,1] = (r&1!=0)*pk
            C[dk,2] = (r&2!=0)*pk
            C[dk,3] = (r&4!=0)*pk
            C[dk,4] = (r&8!=0)*pk
            k += geo(p) + 1
        end
    end
    return
end

geo(p) = floor(UInt64, -randexp() / log1p(-p))