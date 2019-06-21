function build_circulant_matrix(v::Vector{Tv}) where {Tv<:Number}

    # length of the vector
    N = length(v)

    # allocate space for circulant matrix N × N
    C = Matrix{Tv}(undef, N, N)

    # loop over column
    for i2 = 1 : N

        # loop over row, first part
        for i1 = 1 : i2-1
            idx = N - i2 + i1 + 1
            C[i1,i2] = v[idx]
        end

        # loop over row, second part
        for i1 = i2 : N
            idx = i1 - i2 + 1
            C[i1,i2] = v[idx]
        end
    end

    return C
end

# test fast multiplication
N = 111;
s = randn(N);
v = randn(N);
C = build_circulant_matrix(s);

r1 = C * v
r2 = real(ifft( fft(s) .* fft(v) ))

function build_toeplitz_matrix(v::Vector{Tv}) where {Tv<: Number}

    # length of the vector
    N = length(v)

    # compute L, K
    L = floor(Int64, N/2) + 1
    K = N + 1 - L

    # allocate space for Toeplitze matrix
    T = Matrix{Tv}(undef, L, K)

    # loop over column
    for i2 = 1 : K

        # loop over row
        for i1 = 1 : L
            idx = K + i1 - i2
            T[i1,i2] = v[idx]
        end
    end

    return T
end

# compute toeplitz times a vector
function toeplitz_times_vector(c::Vector{Tv}, v::Vector{Tv}) where {Tv<:Number}

    # length of vector
    N = length(c)

    # compute L, K
    L = floor(Int64, N/2) + 1
    K = N + 1 - L

    # check the length of v
    length(v) == K || DimensionMismatch()

    # padding zeros to v
    v_pad = vcat(v, zeros(Tv,L-1))

    # fast computation
    r_pad = real(ifft( fft(c) .* fft(v_pad) ) )

    return r_pad[K:end]

    # # the second way to compute
    # c_hat = Vector{Tv}(undef, N)
    # c_hat[L+1:N] .= c[1:K-1]
    # c_hat[1:L]   .= c[K:N]
    #
    # # fast computation
    # r_pad = real(ifft( fft(c_hat) .* fft(v_pad) ) )
    #
    # return r_pad[1:L]
end

L = 51; K = 50;
N = L+K-1;
c = randn(N); v = randn(K);
r = toeplitz_times_vector(c, v);

T  = build_toeplitz_matrix(c);
r1 = T * v;
norm(r-r1) / norm(r)
