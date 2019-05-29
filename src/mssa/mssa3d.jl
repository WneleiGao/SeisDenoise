"""
   build a level-2 block hankel matrix from a complex matrix
"""
function hankelization(d::Matrix{Complex{Tv}}) where {Tv<:AbstractFloat}

    # size of data
    (n1, n2) = size(d)

    # k is number of columns, l is number of rows
    k1 = floor(Int64, n1/2) + 1; l1 = n1 - k1 + 1;
    k2 = floor(Int64, n2/2) + 1; l2 = n2 - k2 + 1;

    # size of block matrix
    H = zeros(Complex{Tv}, l1*l2, k1*k2)

    for j2 = 1 : k2
        col2 = (j2-1) * k1

        for i2 = 1 : l2
            row2 = (i2-1) * l1
            m2   = i2 + j2 -1

            for j1 = 1 : k1
                col1 = col2 + j1

                for i1 = 1 : l1
                    row1 = row2 + i1
                    m1 = i1 + j1 - 1

                    H[row1, col1] = d[m1, m2]
                end
            end

        end
    end

    return H
end

"""
   recover a matrix by anti-diagonal averaging a block hankel matrix
"""
function anti_diagonal_average(H::Matrix{Complex{Tv}}, n1::Ti, n2::Ti,
         l1::Ti, k1::Ti, l2::Ti, k2::Ti) where {Ti<:Integer, Tv<:AbstractFloat}

    D     = zeros(Complex{Tv}, n1, n2)
    count = zeros(Int64, n1, n2)

    # level2
    for j2 = 1 : k2
        col2 = (j2-1) * k1

        for i2 = 1 : l2
            row2 = (i2-1) * l1
            idx2 = i2 + j2 - 1

            # level 1
            for j1 = 1 : k1
                col1 = col2 + j1

                for i1 = 1 : l1
                    row1 = row2 + i1
                    idx1 = i1 + j1 - 1

                    count[idx1, idx2] = count[idx1, idx2] + 1
                    D[idx1, idx2]     = D[idx1, idx2] + H[row1, col1]
                end
            end
        end
    end

    for j = 1 : n2
        for i = 1 : n1
            D[i,j] = D[i,j] / count[i,j]
        end
    end
    return D
end

"""
   mssa for a 3D cube
"""
function mssa(d::Array{Tv,3}, dt::Tv, rk::Ti; flow=2.0, fhigh=65.0) where{Tv<:AbstractFloat, Ti<:Integer}

    # fourier transform along time axis
    (n1, n2, n3) = size(d)
    k1 = floor(Int64, n2/2)+1; l1 = n2 - k1 + 1
    k2 = floor(Int64, n3/2)+1; l2 = n3 - k2 + 1

    df = fft(d, 1);
    dp = zeros(Complex{Tv}, n1, n2, n3)

    # the start and end index of fequency slice
    ilow  = floor(Int64, flow *dt*n1) + 1
    ihigh = floor(Int64, fhigh*dt*n1) + 1

    inyq  = floor(Int64, n1/2) + 1
    ilow  = ilow  < 2      ? 2      : ilow
    ihigh = ihigh > inyq-1 ? inyq-1 : ihigh

    for i = ilow : ihigh

        # build hankel matrix
        H = hankelization(df[i, :, :])

        # truncated svd
        (u, s, v) = svd(H)
        H = u[:,1:rk] * diagm(0=>s[1:rk]) * (v[:,1:rk])'

        # anti-diagonal averaging
        tmp= anti_diagonal_average(H, n2, n3, l1, k1, l2, k2)
        dp[i,:,:] .= tmp

        # symmetrical property
        dp[n1-i+2,:,:] .= conj(tmp)
    end

    # inverse fft and take the real part
    dp = real(ifft(dp, 1))
    return dp
end
