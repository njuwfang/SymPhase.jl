function _ip(q::SymStabilizer, l_T, r_T)
    cnt = zero(UInt32)
    dr3 = _div3(r_T)
    dr1 = _div1(r_T)
    dl3 = _div3(l_T)
    dl1 = _div1(l_T)
    @inbounds @simd for j1 in axes(q.T, 1)
        cnt ⊻= q.T[j1,dl1,dl3] & q.T_inv[j1,dr1,dr3]
    end

    count_ones(cnt)&1!=0
end

function _display_T(q::SymStabilizer)
    a = Matrix{Bool}(undef, q.nq<<1, q.nq<<1)
    for j2 in 1:q.nq
        for j1 in 1:q.nq
            a[j1,j2] = _ip(q, j1,j2)
            a[j1+q.nq,j2] = _ip(q, j1+q.len3<<_shift3_,j2)
            a[j1,j2+q.nq] = _ip(q, j1,j2+q.len3<<_shift3_)
            a[j1+q.nq,j2+q.nq] = _ip(q, j1+q.len3<<_shift3_,j2+q.len3<<_shift3_)
        end
    end
    display(a)
end

# prime show
function Base.show(io::IO, q::SymStabilizer)
    for row in 1:q.nq
        for col in 1:q.nq
            dx4 = _div4(col)
            dz4 = dx4+q.len4
            dr3 = _div3(row)
            dx2 = _div2(col)
            dz2 = dx2
            dr1 = _div1(row)
            pow = _pow(col)
            x = q.xzs[dr1,dx2,dr3,dx4]&pow!=0
            z = q.xzs[dr1,dz2,dr3,dz4]&pow!=0
            if x && z
                print(io, "Y")
            elseif x
                print(io, "X")
            elseif z
                print(io, "Z")
            else
                print(io, "I")
            end
        end
        println(io)
    end

    for row in 1:q.nq
        for col in 1:q.nq
            dx4 = _div4(col)
            dz4 = dx4+q.len4
            dr3 = _div3(row)
            dx2 = _div2(col)
            dz2 = dx2
            dr1 = _div1(row)
            pow = _pow(col)
            x = q.xzs[dr1,dx2,dr3+q.len3,dx4]&pow!=0
            z = q.xzs[dr1,dz2,dr3+q.len3,dz4]&pow!=0
            if x && z
                print(io, "Y")
            elseif x
                print(io, "X")
            elseif z
                print(io, "Z")
            else
                print(io, "I")
            end
        end
        println(io)
    end
end