using LoopVectorization
using Random: bitrand, rand
using SparseArrays

@eval function Sampler(c::Circuit)
    q = all_zeros(SymStabilizer, c.nqubits, c.nsymbols)
    ns = 0
    nm = 0

    IG = _Vector{Int}()
    JG = _Vector{Int}()
    VG = _Vector{UInt8}()
    b = zeros(Bool, c.nmeasurements)
    op_symbol = OpSymbol()

    op = c.ops[1]
    $(Meta.parse(
    """
        if false
    $(join([
        """
        elseif op == $(key)
            for k in 1:c.target_inds[1]
                apply!(q, $(value)(c.targets[k]))
            end
        """
        for (key, value) in [(X, cX), (Y, cY), (Z, cZ), (H, cHadamard), (S, cPhase)]
        ])
    )
        end
    """
    ))

    if op == DEPOLARIZE1
        l = _div4s(ns+1)
        h = _div4s(ns+c.target_inds[1]<<1)
        transpose_symbols_d!(q, l, h)
        for k in 1:c.target_inds[1]
            apply!(q, cDepolarize1(c.targets[k]), ns+1, ns+2;do_transpose=false)
            ns += 2
            push!(op_symbol, DEPOLARIZE1, c.args[1])
        end
        transpose_symbols_p!(q, l, h)
        update_symbols_range!(q, ns+1-c.target_inds[1]<<1, ns)
    elseif op == CX
        for k in 1:2:c.target_inds[1]
            apply!(q, cCNOT(c.targets[k], c.targets[k+1]))
        end
    elseif op == DEPOLARIZE2
        l = _div4s(ns+1)
        h = _div4s(ns+(c.target_inds[1])<<2)
        transpose_symbols_d!(q, l, h)
        for k in 1:2:c.target_inds[1]
            apply!(q, cDepolarize2(c.targets[k], c.targets[k+1]), ns+1, ns+2, ns+3, ns+4;do_transpose=false)
            ns += 4
            push!(op_symbol, DEPOLARIZE2, c.args[1])
        end
        transpose_symbols_p!(q, l, h)
        update_symbols_range!(q, ns+1-c.target_inds[1]<<2, ns)
    elseif op == M
        transpose_p!(q)
        for k in 1:c.target_inds[1]
            nm += 1
            isrand, phase, i1,i3 = projectZ!(q, c.targets[k], ns+1;do_transpose=false)
            b[nm] = phase
            if isrand
                ns += 1
                push!(op_symbol, M)
                push!(IG, ns)
                push!(JG, nm)
                push!(VG, 1)
            else
                @inbounds for l in q.min_ns[i1,i3]:q.max_ns[i1,i3]
                    if _isone(q, i1, i3, l)
                        push!(IG, l)
                        push!(JG, nm)
                        push!(VG, 1)
                    end
                end
            end
        end
        transpose_d!(q)
    end
    
    for j in 2:c.nops
        op = c.ops[j]
        $(Meta.parse(
        """
            if false
        $(join([
            """
            elseif op == $(key)
                for k in c.target_inds[j-1]+1:c.target_inds[j]
                    apply!(q, $(value)(c.targets[k]))
                end
            """
            for (key, value) in [(X, cX), (Y, cY), (Z, cZ), (H, cHadamard), (S, cPhase)]
            ])
        )
            end
        """
        ))

        if op == DEPOLARIZE1
            l = _div4s(ns+1)
            h = _div4s(ns+(c.target_inds[j]-c.target_inds[j-1])<<1)
            transpose_symbols_d!(q, l, h)
            for k in c.target_inds[j-1]+1:c.target_inds[j]
                apply!(q, cDepolarize1(c.targets[k]), ns+1, ns+2;do_transpose=false)
                ns += 2
                push!(op_symbol, DEPOLARIZE1, c.args[j])
            end
            transpose_symbols_p!(q, l, h)
            update_symbols_range!(q, ns+1-(c.target_inds[j]-c.target_inds[j-1])<<1, ns)
        elseif op == CX
            for k in c.target_inds[j-1]+1:2:c.target_inds[j]
                apply!(q, cCNOT(c.targets[k], c.targets[k+1]))
            end
        elseif op == DEPOLARIZE2
            l = _div4s(ns+1)
            h = _div4s(ns+(c.target_inds[j]-c.target_inds[j-1])<<2)
            transpose_symbols_d!(q, l, h)
            for k in c.target_inds[j-1]+1:2:c.target_inds[j]
                apply!(q, cDepolarize2(c.targets[k], c.targets[k+1]), ns+1, ns+2, ns+3, ns+4;do_transpose=false)
                ns += 4
                push!(op_symbol, DEPOLARIZE2, c.args[j])
            end
            transpose_symbols_p!(q, l, h)
            update_symbols_range!(q, ns+1-(c.target_inds[j]-c.target_inds[j-1])<<2, ns)
        elseif op == M
            transpose_p!(q)
            for k in c.target_inds[j-1]+1:c.target_inds[j]
                nm += 1
                isrand, phase, i1, i3 = projectZ!(q, c.targets[k], ns+1;do_transpose=false)
                b[nm] = phase
                #@show isrand
                if isrand
                    ns += 1
                    push!(op_symbol, M)
                    push!(IG, ns)
                    push!(JG, nm)
                    push!(VG, 1)
                else
                    @inbounds for l in q.min_ns[i1,i3]:q.max_ns[i1,i3]
                        if _isone(q, i1, i3, l)
                            push!(IG, l)
                            push!(JG, nm)
                            push!(VG, 1)
                        end
                    end
                end
            end
            transpose_d!(q)
        end
    end
    
    G = sparse((@view IG.values[1:IG.len]), (@view JG.values[1:JG.len]), (@view VG.values[1:VG.len]), ns, nm)
    #@show IG.len, nm, ns

    sampler = shots -> _Sampler(G, b, op_symbol, shots)
    precompile(sampler, (Int,))

    sampler
end

function _Sampler(G, b, op_symbol, shots)
    nm = size(G, 2)
    xz = biased_random_bits(op_symbol, shots)
    samples = zeros(UInt64, _div64(shots), nm)
    bits_mul_s_uint8s!(samples, xz, G, b)
    samples
end

function bits_mul_s_uint8s!(C::Matrix{UInt64}, A::Matrix{UInt64}, B::AbstractSparseMatrix{UInt8, Int}, b::AbstractVector{Bool})
    rv = rowvals(B)
    M, N = size(C)
    @inbounds for n ∈ 1:N
        temp = ~zero(UInt64)*b[n]
        @turbo for m in 1:M
            C[m, n] = temp
        end
        for k in nzrange(B, n)
            rvk = rv[k]
            @turbo for m in 1:M
                C[m, n] ⊻= A[m, rvk]
            end
        end
    end
end