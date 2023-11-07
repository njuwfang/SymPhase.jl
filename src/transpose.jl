# const _num_packed_bits_ = 128

#primitive type UInt256 256 end

@inline function _transposed_index(j1,j2)
    j = (j2-1)<<8+j1
    j&1023+1, j>>10+1
end

function transpose_p!(q::SymStabilizer)
    for j4 in axes(q.xzs, 4)
        for j3 in axes(q.xzs, 3)
            _transpose_512x64!(q.temp1, (@view q.xzs[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_d!(q::SymStabilizer)
    for j4 in axes(q.xzs, 4)
        for j3 in axes(q.xzs, 3)
            _transpose_64x512!(q.temp1, (@view q.xzs[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_symbols_p!(q::SymStabilizer, l, h)
    for j4 in l:h
        for j3 in axes(q.symbols, 3)
            _transpose_512x1024!(q.temp2, (@view q.symbols[:,:,j3,j4]))
        end
    end
    nothing
end

function transpose_symbols_d!(q::SymStabilizer, l, h)
    for j4 in l:h
        for j3 in axes(q.symbols, 3)
            _transpose_1024x512!(q.temp2, (@view q.symbols[:,:,j3,j4]))
        end
    end
    nothing
end

"""
SIMD for matrix transposition
http://parallelproject15s.appspot.com/2015/05/08/Project-Description.html
"""

@inline function _transpose_512x64!(A, B)
    n = _num_packed_bits_<<1#length(B)>>3
    Ar, Br = A, B
    for j in 0:31
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:15
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:7
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:3
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt128, A), reinterpret(UInt128, B)
    for j in 0:1
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt256, A), reinterpret(UInt256, B)
    for j in 0:0
        interleave!(Br, Ar, j*n, n)
    end
    nothing
end

@inline function _transpose_64x512!(A, B)
    n = _num_packed_bits_<<1#length(B)>>3
    Ar, Br = reinterpret(UInt256, A), reinterpret(UInt256, B)
    for j in 0:0
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt128, A), reinterpret(UInt128, B)
    for j in 0:1
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:3
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:7
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:15
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = A, B
    for j in 0:31
        interleave_i!(Ar, Br, j*n, n)
    end
    nothing
end

@inline function _transpose_512x1024!(A, B)
    n = 1024
    Ar, Br = A, B
    for j in 0:511
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:255
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:127
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:63
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt128, A), reinterpret(UInt128, B)
    for j in 0:31
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt256, A), reinterpret(UInt256, B)
    for j in 0:15
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt512, A), reinterpret(UInt512, B)
    for j in 0:7
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt1024, A), reinterpret(UInt1024, B)
    for j in 0:3
        interleave!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt2048, A), reinterpret(UInt2048, B)
    for j in 0:1
        interleave!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt4096, A), reinterpret(UInt4096, B)
    for j in 0:0
        interleave!(Br, Ar, j*n, n)
    end
    nothing
end

@inline function _transpose_1024x512!(A, B)
    n = 1024
    Ar, Br = reinterpret(UInt4096, A), reinterpret(UInt4096, B)
    for j in 0:0
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt2048, A), reinterpret(UInt2048, B)
    for j in 0:1
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt1024, A), reinterpret(UInt1024, B)
    for j in 0:3
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt512, A), reinterpret(UInt512, B)
    for j in 0:7
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt256, A), reinterpret(UInt256, B)
    for j in 0:15
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt128, A), reinterpret(UInt128, B)
    for j in 0:31
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt64, A), reinterpret(UInt64, B)
    for j in 0:63
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = reinterpret(UInt32, A), reinterpret(UInt32, B)
    for j in 0:127
        interleave_i!(Ar, Br, j*n, n)
    end
    Ar, Br = reinterpret(UInt16, A), reinterpret(UInt16, B)
    for j in 0:255
        interleave_i!(Br, Ar, j*n, n)
    end
    Ar, Br = A, B
    for j in 0:511
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
    end
    @inbounds for j in 1:step
        B[offset+j+step] = A[offset+j<<1]
    end
    nothing
end