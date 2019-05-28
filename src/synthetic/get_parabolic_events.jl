function seis_parab_3D(x1::Vector{Tv}, x2::Vector{Tv};
                     static1=0.0, static2=0.0,
                     ot=0.0, dt =0.004, nt =500, f0=20.0,
                     v1 = [1100, 1500, 1800],
                     v2 = [1100, 1500, 1800],
                     apx1 = [200.0, 300.0, 500.0],
                     apx2 = [ 33.0,  22.0,   0.0],
                     tau  = [0.3,  0.9, 1.3],
                     amp  = [1.0, -1.4, 0.6]) where {Tv<:Real}

    # generate ricker  wavelet
    w = Ricker(dt=dt, f0=f0)

    # number of frequency slice
    nf  = nextpow2(nt)
    # half of frequency slice
    nfh = round(Int64, floor(nf/2)) + 1
    # discrete radian frequency
    wrs = collect(0:1:nfh-1)*2*pi/(nf*dt)

    # allocate space
    nx1 = length(x1)
    nx2 = length(x2)
    D   = zeros(Complex128, nf, nx1, nx2)

    # assimulate statics in shot gather
    if static1 == 0.0
       static1 = zeros(nx1)
    end
    if static2 == 0.0
       static2 = zeros(nx2)
    end
    if length(static1) != nx1 || length(static2) != nx2
      errors("length of static does not match the number of trace")
    end

    # make the wavelet zeros phase
    nw = length(w)
    t_delay = (nw-1)*dt/2

    # pading zeros to wavelet and fft
    w = vcat(w, zeros(nf-nw))
    W = fft(w)

    # loop over events
    nevents = length(amp)
    for ie = 1 : nevents
        time_shift3 = tau[ie] - t_delay
        p1 = sign(v1[ie]) / v1[ie]^2;
        p2 = sign(v2[ie]) / v2[ie]^2;

        for ix2 = 1 : nx2
            time_shift2 = time_shift3 + p2*(x2[ix2]-apx2[ie])^2 + static2[ix2]

            for ix1 = 1 : nx1
                time_shift1 = time_shift2 + p1*(x1[ix1]-apx1[ie])^2 + static1[ix1]

                for iw = 2 : nfh-1
                    phase_shift        = wrs[iw] * time_shift1
                    D[iw,ix1,ix2]     += W[iw] * amp[ie] * exp(-im * phase_shift)
                    D[nf-iw+2,ix1,ix2] = conj(D[iw,ix1,ix2])
                end
            end
        end
    end

    d = real(ifft(D,1))[1:nt,:,:]
    return d

end
