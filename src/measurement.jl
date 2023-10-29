@inline function isX(q, row, col)
    dx4 = _div4(col)
    dr3 = _div3(row)
    dx2 = _div2(col)
    dr1 = _div1(row)
    pow = _pow(col)

    n1,n2,_,_ = size(q.xzs)
    j = (dr1-1)*n2+dx2-1
    j1 = j%n1+1
    j2 = j÷n1+1

    q.xzs[j1,j2,dr3,dx4]&pow!=0
end

@inline function isZ(q, row, col)
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

    q.xzs[j1,j2,dr3,dz4]&pow!=0
end
@inline function isY(q, row, col)
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

    q.xzs[j1,j2,dr3,dx4]&pow!=0 && q.xzs[j1,j2,dr3,dz4]&pow!=0
end
function isXorZ(q, row, col)
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

    q.xzs[j1,j2,dr3,dx4]&pow!=0 || q.xzs[j1,j2,dr3,dz4]&pow!=0
end


function project_cond!(q::SymStabilizer, qubit, symbol_index, cond::Val{IS};do_transpose::Bool) where IS
    anticommutes = 0
    nq = q.nq

    if do_transpose
        transpose_p!(q)
    end

    for j in 1:nq
        if IS(q, j+q.len3<<_shift3_, qubit)
            anticommutes = j
            break
        end
    end

    if anticommutes == 0
        zero!(q, nq+1)
        for j in 1:nq
            if IS(q, j, qubit)
                mul_left!(q, nq+1, j+q.len3<<_shift3_)
            end
        end

        if do_transpose
            transpose_d!(q)
        end

        return false, count_ones(q.phases[_div1(nq+1),_div3(nq+1)])&1!=0, _div1(nq+1),_div3(nq+1)
    else
        rowcopy!(q, anticommutes, q.len3<<_shift3_+anticommutes)
        for j in 1:anticommutes-1
            if IS(q, j, qubit)
                mul_left!(q, j, anticommutes)
            end
        end
        for j in anticommutes+1:nq
            if IS(q, j, qubit)
                mul_left!(q, j, anticommutes)
            end
        end
        for j in anticommutes+1:nq
            if IS(q, j+q.len3<<_shift3_, qubit)
                mul_left!(q, j+q.len3<<_shift3_, anticommutes)
            end
        end
        
        #q[q.len3<<_shift3_+anticommutes, qubit] = (false, true)
        #j1,j2 = _transposed_index(_div1(anticommutes),_div2(symbol_index))
        #q.symbols[j1,j2,_div3(anticommutes)+q.len3,_div4(symbol_index)] ⊻= _pow(symbol_index)
        da3 = _div3(anticommutes)+q.len3
        da1 = _div1(anticommutes)

        ds1 = _div1s(symbol_index)
        ds4 = _div4s(symbol_index)
        pows = _pow(symbol_index)

        zero!(q, anticommutes+q.len3<<_shift3_)
        q[anticommutes+q.len3<<_shift3_, qubit] = (false, true)

        q.symbols[ds1,da1,da3,ds4] = pows

        q.min_ns[da1,da3] = symbol_index
        q.max_ns[da1,da3] = symbol_index

        if do_transpose
            transpose_d!(q)
        end

        return true, false, da1,da3
    end
end

projectX!(q::SymStabilizer, qubit, symbol_index;do_transpose::Bool=true) =
project_cond!(q, qubit, symbol_index, Val(isZ);do_transpose)
projectZ!(q::SymStabilizer, qubit, symbol_index;do_transpose::Bool=true) =
project_cond!(q, qubit, symbol_index, Val(isX);do_transpose)
projectY!(q::SymStabilizer, qubit, symbol_index;do_transpose::Bool=true) =
project_cond!(q, qubit, symbol_index, Val(isXorZ);do_transpose)