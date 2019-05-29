function fxy_eigen(d::Array{Tv,3}, dt::Tv, rk::Ti; flow=2.0, fhigh=65.0) where{Tv<:AbstractFloat, Ti<:Integer}

    # fourier transform along time axis
    (n1, n2, n3) = size(d)

    df = fft(d, 1)
    dp = zeros(Complex{Tv}, n1, n2, n3)

    # the start and end index of fequency slice
    ilow  = floor(Int64, flow *dt*n1) + 1
    ihigh = floor(Int64, fhigh*dt*n1) + 1

    inyq  = floor(Int64, n1/2) + 1
    ilow  = ilow  < 2      ? 2      : ilow
    ihigh = ihigh > inyq-1 ? inyq-1 : ihigh

    for i = ilow : ihigh

        H = df[i, :, :]

        (U, S, V) = svd(H)
        H = u[:,1:rk] * diagm(0=>s[1:rk]) * (v[:,1:rk])'

        dp[i,:,:] = H
        dp[n1-i+2,:,:] = conj(H)
    end

    dp = real(ifft(dp, 1))
    return dp
end
