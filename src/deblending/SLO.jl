function apply_time_shift(d::Matrix{Tv}, dither::Vector, dt) where {Tv<:Real}

    # to eliminate the influence of periodic property of DFT
    # maybe worth to padding zeros to both ends, no implemented at here.

    # get dimension of data
    (n1, n2) = size(d)
    n2 == length(dither) || error("check the length of dither")

    # fourier transform along time axis
    D = fft(d, 1)

    # radian frequency interval
    delta_w = 2 * pi * (1.0 / dt / n1)

    # compute half of the radian frequency axis
    # if n1 is even, nyquist frequency is computed, otherwise not
    nw = floor(Int64, n1/2) + 1
    w  = zeros(Tv, nw)
    for i = 1 : nw
        w[i] = (i-1) * delta_w
    end

    # apply time shift to each trace
    for i2 = 1 : n2
        for i1 = 2 : nw  # start from first non-zero frequency component

            # the conjugate property of real series
            j1 = n1 - i1 + 2
            D[i1,i2] = D[i1,i2] * exp(-im * w[i1] * dither[i2])
            D[j1,i2] = conj(D[i1,i2])
        end
    end

    # transform back to time domain
    return real(ifft(D,1))

end

"""
   From second dimension determine the number of boat
"""
function deblending_by_inversion(d::Matrix{Tv}, dither::Array{Tv1}, trace_length, dt;
         max_iter::Ti=100, alpha=0.9, truncate=[]) where {Tv<:AbstractFloat, Ti<:Int64, Tv1<:Real}

     # number of simultaneous source
     num_src = size(dither, 2)

     # get the dimension of super gather
     (n1, n2) = size(d)

     # recording length (number of samples)
     ns = floor(Int64, trace_length/dt) + 1

     # determine the thresholding curve
     threshold = zeros(Complex{Float64}, max_iter, num_src)
     for i = 1 : num_src  # loop over boat

         # aline traces for current source
         if  norm(dither[:,i]) > 0.0

             ds = apply_time_shift(d, -dither[:,i], dt) # to aline the signal, need to shift in opposite direction

             # apply 2D FT
             Ds = fft(ds)
         end
     end



end





n1 = 256; n2 = 2
dt = 0.002;

# the dither time has to be positive for time domain implementation
# It's better to set the delay time as the integer times of sampling interval
dither = [0.04, 0.08];
d  = zeros(n1, n2);
s1 = ricker(20, dt); l1 = length(s1);
s2 = ricker(30, dt); l2 = length(s2);

d[1:l1,1] .= s1;
d[1:l2,2] .= s2;

# apply time delay in frequency domain (pay attention to periodic property of DFT)
d1 = apply_time_shift(d, dither, dt);

# apply time delay in time domain(assume)
d2 = zeros(eltype(d), n1, n2);
m1 = floor(Int64, dither[1]/dt);
m2 = floor(Int64, dither[2]/dt);

d2[1:end,1] = vcat(zeros(m1), d[1:end-m1,1]);
d2[1:end,2] = vcat(zeros(m2), d[1:end-m2,2]);
norm(d1-d2) / norm(d2)
