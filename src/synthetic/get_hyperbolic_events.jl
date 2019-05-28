function SeisHyp3D(; ot =0.0, dt =0.004, nt =500,
                           ox1=0.0, dx1=40.0 , nx1=51,
                           ox2=0.0, dx2=40.0 , nx2=51,
                           apx1=[0.,0.,0.], apx2=[0.,0.,0.],
                           v1  =[4000, 2000, 1500], v2  =[2000, 2200,2500],
                           tau =[0.3,0.9,1.2], amp=[1.0,-1.0,0.7], f0=20.0)

    w = Ricker(dt=dt,f0=f0)
    nf = nextpow2(nt)
    nw = length(w)
    t_delay = (nw-1)*dt/2     # can show the full wavelet
    w = vcat(w, zeros(nf-nw))
    W = fft(w)
    x1 = ox1 + collect(0:1:nx1-1)*dx1
    x2 = ox2 + collect(0:1:nx2-1)*dx2
    nevents = length(amp)
    D = zeros(Complex128, nf, nx1, nx2)
    nfh = round(Int64, floor(nf/2)) + 1
    wrs = collect(0:1:nfh-1)*2*pi/(nf*dt)     # Frequency in rad/sec
    for ie = 1:nevents
        p1 = sign(v1[ie]) / v1[ie]^2; p2 = sign(v2[ie]) / v2[ie]^2;
        for ix2 = 1:nx2
            for ix1 = 1:nx1
                for iw = 2:nfh-1
                    phase = wrs[iw] * (sqrt(tau[ie]^2 + p1*(x1[ix1]-apx1[ie])^2
                                                      + p2*(x2[ix2]-apx2[ie])^2
                                                      - t_delay))
                    D[iw,ix1,ix2] += W[iw]*amp[ie]*exp(-im * phase)
                    D[nf-iw+2,ix1,ix2] = conj(D[iw,ix1,ix2])
                end
           end
        end
    end

    d = real(ifft(D,1))[1:nt,:,:]
    return d
end
