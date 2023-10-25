# const _num_packed_bits_ = 128

primitive type UInt256 256 end

@inline function _transposed_index(j1,j2)
    j = (j1-1)<<4+j2-1
    j%_num_packed_bits_+1, j÷_num_packed_bits_+1
end

function transpose_p!(q::SymStabilizer)
    for j4 in axes(q.xzs, 4)
        for j3 in axes(q.xzs, 3)
            _transpose__x16!(q.temp, (@view q.xzs[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_d!(q::SymStabilizer)
    for j4 in axes(q.xzs, 4)
        for j3 in axes(q.xzs, 3)
            _transpose_16x_!(q.temp, (@view q.xzs[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_symbols_p!(q::SymStabilizer, min_lens, max_lens)
    for j4 in min_lens:max_lens
        for j3 in axes(q.symbols, 3)
            _transpose__x16!(q.temp, (@view q.symbols[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_symbols_d!(q::SymStabilizer, min_lens, max_lens)
    for j4 in min_lens:max_lens
        for j3 in axes(q.symbols, 3)
            _transpose_16x_!(q.temp, (@view q.symbols[:,:,j3,j4]))
        end
    end
    nothing
end

@inline function _transpose__x16!(A, B)
    n = _num_packed_bits_<<1#length(B)>>3
    Ar, Br = A, B
    for j in 0:7
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:3
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:1
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:0
        interleave!(Br, Ar, j*n, n)
    end
    nothing
end

@inline function _transpose_16x_!(A, B)
    n = _num_packed_bits_<<1#length(B)>>3
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:0
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:1
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:3
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = A, B
    for j in 0:7
        interleave_i!(Ar, Br, j*n, n)
    end
    nothing
end

@inline function interleave!(A, B, offset, n)
    step = n>>1
    @inbounds for j in 1:step
        A[offset+(j-1)<<1+1] = B[offset+j]
        A[offset+j<<1] = B[offset+j+step]
    end
    nothing
end

@inline function interleave_i!(A, B, offset, n)
    step = n>>1
    @inbounds for j in 1:step
        B[offset+j] = A[offset+(j-1)<<1+1]
        B[offset+j+step] = A[offset+j<<1]
    end
    nothing
end