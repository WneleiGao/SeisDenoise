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
function mssa(d::Array{Tv,3}, rk::Ti; dt=0.002, flow=2.0, fhigh=65.0) where{Tv<:AbstractFloat, Ti<:Integer}

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


"""
   divide the cube into patches along spatial direction and apply mssa
to each patches independently.
"""
function local_mssa(cube::Array{Tv,3}, rk::Ti, work_dir::String;
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
     function wrap_mssa(params::Dict)

         # read one patch
         (d, nt, n1, n2, x1l, x1u, x2l, x2u, x1_wo, x2_wo, code) = read_one_spatial_patch(params[:path], return_flag=true)

         # denoising
         s = mssa(d, params[:rk]; dt=params[:dt], flow=params[:flow], fhigh=params[:fhigh])

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
     pmap(wrap_mssa, params)

     # taper the boundary of each patch
     par_spatial_taper(vec_dir)

     # merge patches
     s = spatial_unpatch(vec_dir)

     # remove the patches
     pmap(rm, vec_dir)

     return s
end
