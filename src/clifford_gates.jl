macro SymGate1(name)
    prefixname = Symbol(:c, name)
    quote
        struct $(esc(prefixname))
            q::Int
            $(esc(prefixname))(q) = new(q)
        end
    end
end

macro SymGate2(name)
    prefixname = Symbol(:c, name)
    quote
        struct $(esc(prefixname))
            q1::Int
            q2::Int
            $(esc(prefixname))(q1, q2) = new(q1, q2)
        end
    end
end

@SymGate1 X
@SymGate1 Z
@SymGate1 Y
@SymGate1 Hadamard
@SymGate1 Phase

@SymGate1 X_error
@SymGate1 Z_error
@SymGate1 Y_error
@SymGate1 Depolarize1

@SymGate2 CNOT
@SymGate2 SWAP
@SymGate2 Depolarize2

function apply!(q::SymStabilizer, gate::cHadamard)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    @turbo for j3 in axes(q.xzs, 3)
        for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4]&pow, q.xzs[j1,dz2,j3,dz4]&pow
            q.xzs[j1,dx2,j3,dx4] ⊻= x⊻z
            q.xzs[j1,dz2,j3,dz4] ⊻= x⊻z
            q.phases[j1,j3] ⊻= x&z
        end
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cPhase)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    @turbo for j3 in axes(q.xzs, 3)
        for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4]&pow, q.xzs[j1,dz2,j3,dz4]&pow
            q.xzs[j1,dz2,j3,dz4] ⊻= x
            q.phases[j1,j3] ⊻= x&z
        end
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cX)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4]&pow, q.xzs[j1,dz2,j3,dz4]&pow
            q.phases[j1,j3] ⊻= z
        end
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cX, symbol_index;do_transpose=true)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    ds1 = _div1s(symbol_index)
    ds4 = _div4s(symbol_index)
    shift = _shift(symbol_index) - _shift(a)
    n1,n2,n3,n4 = size(q.xzs)
    m1,m2,m3,m4 = size(q.symbols)
    if do_transpose
        transpose_symbols_d!(q, ds4, ds4)
    end
    for j3 in axes(q.xzs, 3)
        offsetz = (dz4-1)*n3*n2*n1 + (j3-1)*n2*n1 + (dz2-1)*n1
        offsets = (ds4-1)*m3*m1*m2 + (j3-1)*m1*m2 + (ds1-1)*m2
        @inbounds @simd for j1 in axes(q.xzs, 1)
            q.symbols[j1+offsets] = (q.xzs[j1+offsetz]&pow)<<shift
        end
    end
    if do_transpose
        transpose_symbols_p!(q, ds4, ds4)
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cY)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4]&pow, q.xzs[j1,dz2,j3,dz4]&pow
            q.phases[j1,j3] ⊻= x⊻z
        end
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cY, symbol_index;do_transpose=true)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    ds1 = _div1s(symbol_index)
    ds4 = _div4s(symbol_index)
    shift = _shift(symbol_index) - _shift(a)
    n1,n2,n3,n4 = size(q.xzs)
    m1,m2,m3,m4 = size(q.symbols)
    if do_transpose
        transpose_symbols_d!(q, ds4, ds4)
    end
    for j3 in axes(q.xzs, 3)
        offsetx = (dx4-1)*n3*n2*n1 + (j3-1)*n2*n1 + (dx2-1)*n1
        offsetz = (dz4-1)*n3*n2*n1 + (j3-1)*n2*n1 + (dz2-1)*n1
        offsets = (ds4-1)*m3*m1*m2 + (j3-1)*m1*m2 + (ds1-1)*m2
        @inbounds @simd for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1+offsetx], q.xzs[j1+offsetz]
            q.symbols[j1+offsets] = ((x⊻z)&pow)<<shift
        end
    end
    if do_transpose
        transpose_symbols_p!(q, ds4, ds4)
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cZ)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4]&pow, q.xzs[j1,dz2,j3,dz4]&pow
            q.phases[j1,j3] ⊻= x
        end
    end
    nothing
end

function apply!(q::SymStabilizer, gate::cZ, symbol_index;do_transpose=true)
    a = gate.q
    dx4 = _div4(a)
    dx2 = _div2(a)
    pow = _pow(a)
    ds1 = _div1s(symbol_index)
    ds4 = _div4s(symbol_index)
    shift = _shift(symbol_index) - _shift(a)
    n1,n2,n3,n4 = size(q.xzs)
    m1,m2,m3,m4 = size(q.symbols)
    if do_transpose
        transpose_symbols_d!(q, ds4, ds4)
    end
    for j3 in axes(q.xzs, 3)
        offsetx = (dx4-1)*n3*n2*n1 + (j3-1)*n2*n1 + (dx2-1)*n1
        offsets = (ds4-1)*m3*m1*m2 + (j3-1)*m1*m2 + (ds1-1)*m2
        @inbounds @simd for j1 in axes(q.xzs, 1)
            q.symbols[j1+offsets] = (q.xzs[j1+offsetx]&pow)<<shift
        end
    end
    if do_transpose
        transpose_symbols_p!(q, ds4, ds4)
    end

    nothing
end

function apply!(q::SymStabilizer, gate::cCNOT)
    a1 = gate.q1
    a2 = gate.q2

    d1x4 = _div4(a1)
    d1z4 = d1x4+q.len4
    d1x2 = _div2(a1)
    d1z2 = d1x2
    pow1 = _pow(a1)
    d2x4 = _div4(a2)
    d2z4 = d2x4+q.len4
    d2x2 = _div2(a2)
    d2z2 = d2x2
    pow2 = _pow(a2)

    @turbo for j3 in axes(q.xzs, 3)
        for j1 in axes(q.xzs, 1)
            x1, z1 = q.xzs[j1,d1x2,j3,d1x4]&pow1!=0, q.xzs[j1,d1z2,j3,d1z4]&pow1!=0
            x2, z2 = q.xzs[j1,d2x2,j3,d2x4]&pow2!=0, q.xzs[j1,d2z2,j3,d2z4]&pow2!=0
            q.xzs[j1,d2x2,j3,d2x4] ⊻= x1*pow2
            q.xzs[j1,d1z2,j3,d1z4] ⊻= z2*pow1
            q.phases[j1,j3] ⊻= (x1&z1&x2&z2)|(x1&z2&(~(z1|x2)))
        end
    end

    nothing
end

function apply!(q::SymStabilizer, gate::cSWAP)
    a1 = gate.q1
    a2 = gate.q2

    d1x4 = _div4(a1)
    d1z4 = d1x4+q.len4
    d1x2 = _div2(a1)
    d1z2 = d1x2
    pow1 = _pow(a1)
    d2x4 = _div4(a2)
    d2z4 = d2x4+q.len4
    d2x2 = _div2(a2)
    d2z2 = d2x2
    pow2 = _pow(a2)

    if d1x2 == d2x2 && d1x4 == d2x4
        @turbo for j3 in axes(q.xzs, 3)
            for j1 in axes(q.xzs, 1)
                x1, z1 = q.xzs[j1,d1x2,j3,d1x4]&pow1!=0, q.xzs[j1,d1z2,j3,d1z4]&pow1!=0
                x2, z2 = q.xzs[j1,d2x2,j3,d2x4]&pow2!=0, q.xzs[j1,d2z2,j3,d2z4]&pow2!=0

                q.xzs[j1,d1x2,j3,d1x4] ⊻= (x1⊻x2)*(pow1⊻pow2)
                q.xzs[j1,d1z2,j3,d1z4] ⊻= (z1⊻z2)*(pow1⊻pow2)
                
            end
        end
    else
        @turbo for j3 in axes(q.xzs, 3)
            for j1 in axes(q.xzs, 1)
                x1, z1 = q.xzs[j1,d1x2,j3,d1x4]&pow1!=0, q.xzs[j1,d1z2,j3,d1z4]&pow1!=0
                x2, z2 = q.xzs[j1,d2x2,j3,d2x4]&pow2!=0, q.xzs[j1,d2z2,j3,d2z4]&pow2!=0
                
                q.xzs[j1,d1x2,j3,d1x4] ⊻= (x1⊻x2)*pow1
                q.xzs[j1,d2x2,j3,d2x4] ⊻= (x1⊻x2)*pow2

                q.xzs[j1,d1z2,j3,d1z4] ⊻= (z1⊻z2)*pow1
                q.xzs[j1,d2z2,j3,d2z4] ⊻= (z1⊻z2)*pow2
            end
        end
    end

    nothing
end

apply!(q::SymStabilizer, gate::cX_error, symbol_index::Int) = apply!(q, cX(gate.q), symbol_index)
apply!(q::SymStabilizer, gate::cY_error, symbol_index::Int) = apply!(q, cY(gate.q), symbol_index)
apply!(q::SymStabilizer, gate::cZ_error, symbol_index::Int) = apply!(q, cZ(gate.q), symbol_index)

function apply!(q::SymStabilizer, gate::cDepolarize1, i1::Int, i2::Int;do_transpose=true)
    if do_transpose
        l = _div4(min(i1,i2))
        h = _div4(max(i1,i2))

        transpose_symbols_d!(q, l, h)

        apply!(q, cX(gate.q), i1;do_transpose=false)
        apply!(q, cZ(gate.q), i2;do_transpose=false)

        transpose_symbols_p!(q, l, h)
    else
        apply!(q, cX(gate.q), i1;do_transpose=false)
        apply!(q, cZ(gate.q), i2;do_transpose=false)
    end

    nothing
end

function apply!(q::SymStabilizer, gate::cDepolarize2, i1::Int, i2::Int, i3::Int, i4::Int;do_transpose=true)
    if do_transpose
        l = _div4(min(i1,i2,i3,i4))
        h = _div4(max(i1,i2,i3,i4))

        transpose_symbols_d!(q, l, h)

        apply!(q, cX(gate.q1), i1;do_transpose=false)
        apply!(q, cZ(gate.q1), i2;do_transpose=false)
        apply!(q, cX(gate.q2), i3;do_transpose=false)
        apply!(q, cZ(gate.q2), i4;do_transpose=false)

        transpose_symbols_p!(q, l, h)
    else
        apply!(q, cX(gate.q1), i1;do_transpose=false)
        apply!(q, cZ(gate.q1), i2;do_transpose=false)
        apply!(q, cX(gate.q2), i3;do_transpose=false)
        apply!(q, cZ(gate.q2), i4;do_transpose=false)
    end
    
    nothing
end