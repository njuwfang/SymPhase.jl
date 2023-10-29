function mul_left!(q::SymStabilizer, r, l)
    n1,_,n3,n4 = size(q.xzs)
    xzs = reshape(reinterpret(UInt512, q.xzs), n1,n3,n4)

    dl1 = _div1(l)
    dl3 = _div3(l)
    dr1 = _div1(r)
    dr3 = _div3(r)
    cnt1 = zero(UInt512)
    cnt2 = zero(UInt512)

    @inbounds for j4 in 1:q.len4
        x1, z1 = xzs[dr1, dr3, j4], xzs[dr1, dr3, j4+q.len4]
        x2, z2 = xzs[dl1, dl3, j4], xzs[dl1, dl3, j4+q.len4]
        xzs[dr1, dr3, j4] = newx1 = x1 ⊻ x2
        xzs[dr1, dr3, j4+q.len4] = newz1 = z1 ⊻ z2
        x1z2 = x1 & z2
        anti_comm = (x2 & z1) ⊻ x1z2
        cnt2 ⊻= (cnt1 ⊻ newx1 ⊻ newz1 ⊻ x1z2) & anti_comm
        cnt1 ⊻= anti_comm
    end

    extra_phase = ((count_ones(cnt1)⊻(count_ones(cnt2)<<1))&0x3)==0x2
    
    q.phases[dr1,dr3] ⊻= extra_phase ⊻ q.phases[dl1, dl3]

    for j4 in _div4s(q.min_ns[dl1,dl3]):_div4s(q.max_ns[dl1,dl3])
        @turbo for j1 in axes(q.symbols, 1)
            q.symbols[j1,dr1,dr3,j4] ⊻= q.symbols[j1,dl1,dl3,j4]
        end
    end

    q.min_ns[dr1,dr3] = min(q.min_ns[dr1,dr3], q.min_ns[dl1,dl3])
    q.max_ns[dr1,dr3] = max(q.max_ns[dr1,dr3], q.max_ns[dl1,dl3])

    nothing
end