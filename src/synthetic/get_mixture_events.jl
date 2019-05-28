"""
   get ricker wavelet with two input arguments (fdom: dominant frequency and time sampling rate dt)
"""
function ricker(fdom, dt)

 	  nw = 2.2 / fdom / dt         # decide number of samples
	  nw = 2*floor(Int64, nw/2)+1  # gaurantee the number of samples is odd

	  nc = floor(Int64, nw/2)
	  w  = zeros(nw)               # allocate space for wavelet
	  k  = collect(1:nw)

    for i = 1 : nw
        alpha = (nc-i+1) * fdom * dt * pi
        beta  = alpha^2
        w[i]  = (1 - beta * 2) * exp(-beta)
    end
	  return w
end

"""
   get seismic cube (up to 5D) consisted of mixtured events(linear, hyperbolic, parabolic)
"""
function get_mixture_events(; ot =0.0, dt =0.004, nt =500,
                              ox1=0.0, dx1=10.0 , nx1=200,
                              ox2=0.0, dx2=10.0 , nx2=200,
                              ox3=0.0, dx3=10.0 , nx3=1,
                              ox4=0.0, dx4=10.0 , nx4=1,
                              apx1=[0.,0.,0.], apx2=[0.,0.,0.],
                              v1  =[6000 , 2000, -2700],
                              v2  =[14000, 4000,  4000],
                              tau =[0.1, 0.4 , 0.9],
                              amp =[1.0, -1.0, 0.7],
                              event_type=['l', 'p', 'h'],
                              static_flag = true, sigma=0.075, L=13,
                              f0=20.0, data_format=Float64)

    # check the length of smoother
    if mod(L, 2) == 0
       @warn "L must be an odd number"
       L = L + 1
    end

    # source wavelet
    w  = ricker(f0, dt)
    nw = length(w)

    # number of frequency slice
    nf = nextpow(2, nt)
    if nf < nw
       error("number of samples is too small")
    end
    df = 1.0 / dt / nf

    # padding zeros to source wavelet
    w = vcat(w, zeros(nf-nw))

    # fourier transform
    fw = fft(w)

    # spatial axis
    x1 = ox1 .+ collect(0 : 1 : nx1-1) * dx1
    x2 = ox2 .+ collect(0 : 1 : nx2-1) * dx2

    # number of events
    nevents = length(amp)
    if nevents != length(apx1) || nevents != length(apx2) || nevents != length(v1) || nevents != length(v2) || nevents != length(tau) || nevents != length(event_type)
       error("input parameter does not match")
    end

    # allocate memory for data in frequency domain
    d = zeros(Complex{Float64}, nf, nx1, nx2)

    # index of nyquist frequency
    nfh = floor(Int64, nf/2) + 1

    # frequency upper limite
    fhigh = f0 * 2.5
    iw_upper = floor(Int64, fhigh/df) + 1
    if iw_upper > nfh
       iw_upper = nfh
    end

    # frequency lower limite (default to be df)
    iw_lower = 2

    #  random statics
    pt = zeros(nevents, nx1, nx2)
    if static_flag
       smoother = hamming(L)
       smoother = smoother / sum(smoother)
       drop     = floor(Int64, L/2)*2

       for ie = 1 : nevents
           tmp = conv2(smoother, smoother, randn(nx1+drop, nx2+drop))[drop+1:end-drop, drop+1:end-drop] * sigma
           pt[ie,:,:] .= tmp[:,:]
       end
    end

    # loop over all the traces
    for ix2 = 1 : nx2
        for ix1 = 1 : nx1

            # loop over frequency
            for iw = iw_lower : iw_upper

                # radian frequency
                omega = 2.0 * pi * (iw-1)*df

                # loop over events
                for ie = 1 : nevents

                    # linear event
                    if event_type[ie] == 'l'
                       p1 = 1.0 / v1[ie]
                       p2 = 1.0 / v2[ie]
                       s1 = abs(x1[ix1]-apx1[ie])
                       s2 = abs(x2[ix2]-apx2[ie])

                    # curve event
                  elseif event_type[ie] == 'p' || event_type[ie] == 'h'
                       p1 = sign(v1[ie]) / v1[ie]^2
                       p2 = sign(v2[ie]) / v2[ie]^2
                       s1 = (x1[ix1]-apx1[ie])^2
                       s2 = (x2[ix2]-apx2[ie])^2
                    end

                    # phase shift corresponding to time delay
                    if event_type[ie] == 'l' || event_type[ie] == 'p'
                       phase = omega * (tau[ie] + p1*s1 + p2*s2 + pt[ie,ix1,ix2])
                    elseif event_type[ie] == 'h'
                       phase = omega * (sqrt(abs(tau[ie] + p1*s1 + p2*s2)) + pt[ie,ix1,ix2])
                    end

                    d[iw,ix1,ix2] += fw[iw] * amp[ie] * exp(-im * phase)
                end

                # conjugacy
                d[nf-iw+2,ix1,ix2] = conj(d[iw,ix1,ix2])
            end
        end
    end

    # get the data in time domain
    d = ifft(d, 1)
    return convert(Array{data_format,3}, real(d[1:nt,:,:]))

end

# testing the code
# d = get_mixture_events(nt=250, nx1=200, nx2=200, dt=0.008, dx1=10, dx2=10,
#                               v1=[6000. ,2000,-2700.,-2700.,6000. ],
#                               v2=[14000.,4000, 4000., 4000.,14000.],
#                               tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
#                               amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
#                               apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
#                               apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
#                               event_type=['l', 'h', 'p', 'p', 'l' ],
#                               f0=20.0, static_flag=true, L=13, sigma=0.1)
