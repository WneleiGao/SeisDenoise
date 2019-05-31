function fxy_eigen(d::Array{Tv,3}, rk::Ti; dt=00.2, flow=2.0, fhigh=65.0) where{Tv<:AbstractFloat, Ti<:Integer}

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

        (u, s, v) = svd(H)
        H = u[:,1:rk] * diagm(0=>s[1:rk]) * (v[:,1:rk])'

        dp[i,:,:] .= H
        dp[n1-i+2,:,:] .= conj(H)
    end

    dp = real(ifft(dp, 1))
    return dp
end

"""
   divide the cube into patches along spatial direction and apply mssa
to each patches independently.
"""
function local_fxy_eigen(cube::Array{Tv,3}, rk::Ti, work_dir::String;
                         dt=0.004, flow=2.0, fhigh=60.0,
                         x1_wl = 30, x1_wo = 10,
                         x2_wl = 30, x2_wo = 10) where {Tv<:AbstractFloat, Ti<:Int64}

     # divide the cube into patches
     vec_dir = spatial_patch(work_dir, cube, x1_wl, x1_wo, x2_wl, x2_wo);

     # pack arguments into a Dictionary
     num_patches = length(vec_dir)
     params = Vector{Dict}(undef, num_patches)
     for i = 1 : num_patches
         params[i] = Dict(:path=>vec_dir[i], :rk=>rk, :dt=>dt, :flow=>flow, :fhigh=>fhigh)
     end

     # wrap_fxy_prediction
     function wrap_fxy_eigen(params::Dict)

         # read one patch
         (d, nt, n1, n2, x1l, x1u, x2l, x2u, x1_wo, x2_wo, code) = read_one_spatial_patch(params[:path], return_flag=true)

         # denoising
         s = fxy_eigen(d, params[:rk]; dt=params[:dt], flow=params[:flow], fhigh=params[:fhigh])

         # write the result back
         fid = open(params[:path], "r+")
         pos = sizeof(Int64) * 10
         seek(fid, pos)

         if code == 1
            write(fid, convert(Vector{Float64}, vec(s)))
         elseif code == 2
            write(fid, convert(Vector{Float32}, vec(s)))
         end
         close(fid)

         return nothing
     end

     # apply fxy_prediction to each patch
     pmap(wrap_fxy_eigen, params)

     # taper the boundary of each patch
     par_spatial_taper(vec_dir)

     # merge patches
     s = spatial_unpatch(vec_dir)

     # remove the patches
     pmap(rm, vec_dir)

     return s
end
