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

function apply!(q::SymStabilizer, gate::cX, symbol_index)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            q.temp_phases[j1,j3] = q.xzs[j1,dz2,j3,dz4]&pow!=0
        end
    end

    if q.enable_T
        @inbounds for j3 in axes(q.T_inv, 3)
            for j2 in axes(q.T_inv, 2)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                @turbo for j1 in axes(q.symbols, 1)
                    q.symbols[j1,symbol_index] ⊻= q.T_inv[j1,j2,j3]
                end
            end
            end
        end
    else
        ds = _div32(symbol_index)
        pows = _pow32(symbol_index)
        @inbounds for j3 in axes(q.temp_phases, 2)
            for j2 in axes(q.temp_phases, 1)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                    q.symbols[ds,j2,j3] ⊻= pows
                end
            end
        end
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

function apply!(q::SymStabilizer, gate::cY, symbol_index)
    a = gate.q
    dx4 = _div4(a)
    dz4 = dx4+q.len4
    dx2 = _div2(a)
    dz2 = dx2
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            x, z = q.xzs[j1,dx2,j3,dx4], q.xzs[j1,dz2,j3,dz4]
            q.temp_phases[j1,j3] = (x⊻z)&pow!=0
        end
    end

    if q.enable_T
        @inbounds for j3 in axes(q.T_inv, 3)
            for j2 in axes(q.T_inv, 2)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                @turbo for j1 in axes(q.symbols, 1)
                    q.symbols[j1,symbol_index] ⊻= q.T_inv[j1,j2,j3]
                end
            end
            end
        end
    else
        ds = _div32(symbol_index)
        pows = _pow32(symbol_index)
        @inbounds for j3 in axes(q.temp_phases, 2)
            for j2 in axes(q.temp_phases, 1)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                    q.symbols[ds,j2,j3] ⊻= pows
                end
            end
        end
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

function apply!(q::SymStabilizer, gate::cZ, symbol_index)
    a = gate.q
    dx4 = _div4(a)
    dx2 = _div2(a)
    pow = _pow(a)
    for j3 in axes(q.xzs, 3)
        @inbounds @simd for j1 in axes(q.xzs, 1)
            q.temp_phases[j1,j3] = q.xzs[j1,dx2,j3,dx4]&pow!=0
        end
    end

    if q.enable_T
        @inbounds for j3 in axes(q.T_inv, 3)
            for j2 in axes(q.T_inv, 2)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                @turbo for j1 in axes(q.symbols, 1)
                    q.symbols[j1,symbol_index] ⊻= q.T_inv[j1,j2,j3]
                end
            end
            end
        end
    else
        ds = _div32(symbol_index)
        pows = _pow32(symbol_index)
        @inbounds for j3 in axes(q.temp_phases, 2)
            for j2 in axes(q.temp_phases, 1)
                if q.temp_phases[j2,j3]
                    q.min_ns[j2,j3] = min(q.min_ns[j2,j3], symbol_index)
                    q.max_ns[j2,j3] = max(q.max_ns[j2,j3], symbol_index)
                    q.symbols[ds,j2,j3] ⊻= pows
                end
            end
        end
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

apply!(q::SymStabilizer, gate::cX_error, symbol_index::Int) = apply!(q, cX(gate.q), symbol_index)
apply!(q::SymStabilizer, gate::cY_error, symbol_index::Int) = apply!(q, cY(gate.q), symbol_index)
apply!(q::SymStabilizer, gate::cZ_error, symbol_index::Int) = apply!(q, cZ(gate.q), symbol_index)

function apply!(q::SymStabilizer, gate::cDepolarize1, i1::Int, i2::Int)
    apply!(q, cX(gate.q), i1)
    apply!(q, cZ(gate.q), i2)
end

function apply!(q::SymStabilizer, gate::cDepolarize2, i1::Int, i2::Int, i3::Int, i4::Int)
    apply!(q, cX(gate.q1), i1)
    apply!(q, cZ(gate.q1), i2)
    apply!(q, cX(gate.q2), i3)
    apply!(q, cZ(gate.q2), i4)
end