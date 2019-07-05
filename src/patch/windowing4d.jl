"""
   get the directory of all patches, return a vector of String
"""
function get_patches_dir(path::String, it_nw::Ti, x1_nw::Ti, x2_nw::Ti, x3_nw::Ti) where {Ti<:Int64}

    dir = Vector{String}(undef, it_nw * x1_nw * x2_nw * x3_nw)

    for i3 = 1 : x3_nw
        idx3 = (i3-1) * x2_nw * x1_nw * it_nw

        for i2 = 1 : x2_nw
            idx2 = idx3 + (i2-1) * x1_nw * it_nw

            for i1 = 1 : x1_nw
                idx1 = idx2 + (i1-1) * it_nw

                for it = 1 : it_nw
                    idx= idx1 + it

                    file_name = join(["patch" "_" it "_" i1 "_" i2 "_" i3 ".bin"])
                    dir[idx]  = joinpath(path, file_name)
                end
            end
        end
    end

    return dir
end

"""
   d is a 3D array, it_wl is the window length along time axis, it_wo is the
overlapping of window.
"""
function patch4d(path::String, d::Array{Tv,4},
               it_wl::Ti, it_wo::Ti,
               x1_wl::Ti, x1_wo::Ti,
               x2_wl::Ti, x2_wo::Ti,
               x3_wl::Ti, x3_wo::Ti) where {Tv<:AbstractFloat, Ti<:Int64}

    # size of data
    (nt, n1, n2, n3) = size(d)

    # shift in each direction
    it_ws = it_wl - it_wo
    x1_ws = x1_wl - x1_wo
    x2_ws = x2_wl - x2_wo
    x3_ws = x3_wl - x3_wo

    # make sure the window size not too big
    it_wl = it_wl > nt ? nt : it_wl
    x1_wl = x1_wl > n1 ? n1 : x1_wl
    x2_wl = x2_wl > n2 ? n2 : x2_wl
    x3_wl = x3_wl > n3 ? n3 : x3_wl

    # number of windows in each direction
    it_nw = floor(Int64,(nt-it_wl)/it_ws)+1; it_nw = (it_nw-1)*it_ws+it_wl < nt ? it_nw+1 : it_nw
    x1_nw = floor(Int64,(n1-x1_wl)/x1_ws)+1; x1_nw = (x1_nw-1)*x1_ws+x1_wl < n1 ? x1_nw+1 : x1_nw
    x2_nw = floor(Int64,(n2-x2_wl)/x2_ws)+1; x2_nw = (x2_nw-1)*x2_ws+x2_wl < n2 ? x2_nw+1 : x2_nw
    x3_nw = floor(Int64,(n3-x3_wl)/x3_ws)+1; x3_nw = (x3_nw-1)*x3_ws+x3_wl < n3 ? x3_nw+1 : x3_nw

    # third spatial direction
    for i3 = 1 : x3_nw
        i3l = (i3-1)*x3_ws + 1
        i3u =  i3l  +x3_wl - 1
        if i3u > n3
           i3u = n3
        end

        # second spatial direction
        for i2 = 1 : x2_nw
            i2l = (i2-1)*x2_ws + 1
            i2u =  i2l  +x2_wl - 1
            if i2u > n2
               i2u = n2
            end

            # first spatial direction
            for i1 = 1 : x1_nw
                i1l = (i1-1)*x1_ws + 1
                i1u = i1l  + x1_wl - 1
                if i1u > n1
                   i1u = n1
                end

                # time direction
                for it = 1 : it_nw
                    itl = (it-1)*it_ws + 1
                    itu = itl   +it_wl - 1
                    if itu > nt
                       itu = nt
                    end

                    # path to one binary file
                    file_name = join(["patch" "_" it "_" i1 "_" i2 "_" i3 ".bin"])
                    pout      = joinpath(path, file_name)

                    # dimensions of the patch
                    fid = open(pout, "w")
                    write(fid, nt, n1, n2, n3)
                    write(fid, itl, itu, i1l, i1u, i2l, i2u, i3l, i3u)
                    write(fid, it_wo, x1_wo, x2_wo, x3_wo)

                    # data format code
                    if eltype(d) == Float64
                       write(fid, 1)
                    elseif eltype(d) == Float32
                       write(fid, 2)
                    end

                    # data save as Float32
                    tmp = d[itl:itu, i1l:i1u, i2l:i2u, i3l:i3u]
                    write(fid, vec(tmp))

                    close(fid)
                end
            end
        end
    end

    return get_patches_dir(path, it_nw, x1_nw, x2_nw, x3_nw)
end

"""
   read one patch of a cube
"""
function read_one_patch4d(pin::String; return_flag=false)

    fid = open(pin, "r")

    # size of original cube
    nt  = read(fid, Int64)
    n1  = read(fid, Int64)
    n2  = read(fid, Int64)
    n3  = read(fid, Int64)

    # size of current patch
    itl = read(fid, Int64); itu = read(fid, Int64); nt_w = itu - itl + 1;
    x1l = read(fid, Int64); x1u = read(fid, Int64); n1_w = x1u - x1l + 1;
    x2l = read(fid, Int64); x2u = read(fid, Int64); n2_w = x2u - x2l + 1;
    x3l = read(fid, Int64); x3u = read(fid, Int64); n3_w = x3u - x3l + 1;

    # window overlap
    it_wo = read(fid, Int64)
    x1_wo = read(fid, Int64)
    x2_wo = read(fid, Int64)
    x3_wo = read(fid, Int64)

    # data format code
    code = read(fid, Int64)

    # read the data
    if code == 1
       d = zeros(Float64, nt_w * n1_w * n2_w * n3_w)

    elseif code == 2
       d = zeros(Float32, nt_w * n1_w * n2_w * n3_w)
    end

    read!(fid, d)
    d = reshape(d,  nt_w, n1_w, n2_w, n3_w)
    close(fid)

    # return
    if return_flag
       return d, nt, n1, n2, n3, itl, itu, x1l, x1u, x2l, x2u, x3l, x3u, it_wo, x1_wo, x2_wo, x3_wo, code
    else
       return d
    end
end

"""
   taper the boundary part of one patch
"""
function taper_one_patch4d(pin::String) where {Ti<:Int64}

    # read one patch
    (d, nt, n1, n2, n3, itl, itu, x1l, x1u, x2l, x2u, x3l, x3u, it_wo, x1_wo, x2_wo, x3_wo, code) = read_one_patch4d(pin, return_flag=true)

    # window size of one patch
    nt_w = itu-itl+1;
    n1_w = x1u-x1l+1;
    n2_w = x2u-x2l+1;
    n3_w = x3u-x3l+1;

    tapti = itl >  1  ? it_wo : 0
    taptf = itu <  nt ? it_wo : 0

    tap1i = x1l >  1  ? x1_wo : 0
    tap1f = x1u <  n1 ? x1_wo : 0

    tap2i = x2l >  1  ? x2_wo : 0
    tap2f = x2u <  n2 ? x2_wo : 0

    tap3i = x3l >  1  ? x3_wo : 0
    tap3f = x3u <  n3 ? x3_wo : 0

    tt = 1.; t1 = 1.; t2 = 1.; t3 = 1.;

    # taper in the second spatial direction
    for i3 = 1 : n3_w
        if i3 >= 1 && i3 <= tap3i
           t3 = 1. - cos(pi/2*((i3-1)/tap3i))
        end
        if i3 > tap3i && i3 <= n3_w-tap3f
           t3 = 1.
        end
        if i3 > n3_w-tap3f && i3 <= n3_w
           t3 = cos(pi/2*(i3-1-n3_w+tap3f)/tap3f)
        end

        for i2 = 1 : n2_w
            if i2 >= 1 && i2 <= tap2i
               t2 = 1. - cos(pi/2*((i2-1)/tap2i))
            end
            if i2 > tap2i && i2 <= n2_w-tap2f
               t2 = 1.
            end
            if i2 > n2_w-tap2f && i2 <= n2_w
               t2 = cos(pi/2*(i2-1-n2_w+tap2f)/tap2f)
            end

            # taper in the first spatial direction
            for i1 = 1 : n1_w
                if i1 >= 1 && i1 <= tap1i
                   t1 = 1. - cos(pi/2*((i1-1)/tap1i))
                end
                if i1 > tap1i && i1 <= n1_w-tap1f
                   t1 = 1.
                end
                if i1 > n1_w-tap1f && i1 <= n1_w
                   t1 = cos(pi/2*(i1-1-n1_w+tap1f)/tap1f)
                end

                # taper in the time axis
                for it = 1 : nt_w
                    if it >= 1 && it <= tapti
                       tt = 1. - cos(pi/2*((it-1)/tapti))
                    end
                    if it > tapti && it <= nt_w-taptf
                       tt = 1.
                    end
                    if it > nt_w-taptf && it <= nt_w
                       tt = cos(pi/2*(it-1-nt_w+taptf)/taptf)
                    end

                    d[it,i1,i2,i3] = d[it,i1,i2,i3] * tt * t1 * t2 * t3

                end
            end
        end
    end

    fid = open(pin, "r+")
    pos = sizeof(Int64) * 17
    seek(fid, pos)

    write(fid, vec(d))
    close(fid)
    return nothing
end

"""
   tapering the overlapping boundary of each patch
"""
function par_taper4d(path::Vector{String})

    pmap(taper_one_patch4d, path)
    return nothing
end

"""
   sum over the patches to get back the original cube
"""
function unpatch4d(pin::Vector{String})

    (d1, nt, n1, n2, n3, itl, itu, x1l, x1u, x2l, x2u, x3l, x3u, it_wo, x1_wo, x2_wo, x3_wo, code)= read_one_patch4d(pin[1]; return_flag=true);

    if code == 1
       d = zeros(Float64, nt, n1, n2, n3)
    elseif code == 2
       d = zeros(Float32, nt, n1, n2, n3)
    end

    d[itl:itu, x1l:x1u, x2l:x2u, x3l:x3u] .= d1[:,:,:,:]

    for i = 2 : length(pin)
        (d1, nt, n1, n2, n3, itl, itu, x1l, x1u, x2l, x2u, x3l, x3u, it_wo, x1_wo, x2_wo, x3_wo, code)= read_one_patch4d(pin[i], return_flag=true);
        d[itl:itu, x1l:x1u, x2l:x2u, x3l:x3u] .= d[itl:itu, x1l:x1u, x2l:x2u, x3l:x3u] .+ d1[:,:,:,:]
    end

    return d
end

# # generate a 4D cube
# path = joinpath(homedir(), "Desktop/randCP");
# d0   = randn(Float64, 55, 56, 57, 58);
#
# # divided into patches
# it_wl=17; it_wo=4;
# x1_wl=16; x1_wo=5;
# x2_wl=21; x2_wo=5;
# x3_wl=31; x3_wo=5;
# vec_dir = patch4d(path, d0, it_wl, it_wo, x1_wl, x1_wo, x2_wl, x2_wo, x3_wl, x3_wo);
#
# # taper the boundary of each patch
# par_taper4d(vec_dir);
#
# # merge patches
# d1 = unpatch4d(vec_dir);
#
# # test the difference
# norm(d1-d0) / norm(d0)
