"""
   struct for tucker tensor
"""
struct Tucker{Tv<:Number, Ti<:Int64}
    order :: Ti
    rk    :: Vector{Ti}
    dims  :: Vector{Ti}
    core  :: Array{Tv}
    u     :: Vector{Matrix{Tv}}
end

"""
   constructor of tucker tensor with provided core tensor and factor matrices
"""
function Tucker(core::Array{Tv}, u::Vector{Matrix{Tv}}) where {Tv<:Number}

    # check dimensions
    order = ndims(core)
    rk    = collect(size(core))

    order == length(u) || error("core dimension does not match number of factors")
    dims  =  zeros(Int64, order)

    # make sure the size of u
    for i = 1 : order
        rk[i] == size(u[i], 2) || error("core dimension does not match factors")
        dims[i]= size(u[i], 1)
    end

    return Tucker(order, rk, dims, copy(core), deepcopy(u))
end

"""
   copy a tucker tensor
"""
function copy_tucker(t::Tucker)
    return Tucker(t.order, copy(t.rk), copy(t.dims), copy(t.core), deepcopy(t.u))
end

"""
   convert low-rank tucker tensor to dense tensor
"""
function tucker2tensor(t::Tucker)

    # create a core tensor
    c = Tensor(t.core)

    # tensor matrix multiplication
    n = collect(1 : t.order)
    return ttm(t.u, c, n)
end

"""
   compute the n^{th} mode factor matrix by keeping the first rk eigenvectors
"""
function nvecs(t::Tensor, n::Ti, rk::Ti) where {Ti<:Int64}

    tn = matricization(t, n)

    # shrink the size of matrix
    y  = tn * tn'
    (u, s, v) = svd(y)

    # return the first rk left eigenvectors
    return u[:,1:rk]
end

"""
   tucker decomposition by alternative minimization
"""
function tucker_als(t::Tensor{Tv,Ti}, rk::Union{Vector{Ti},Ti};
                   tol=1e-6, max_iter=10, init_flag="eigvec") where {Tv<:Number, Ti<:Int64}

    order = t.order
    dims  = t.dims

    # Frobenous norm of original tensor
    normt = norm(t.d)

    # identical rank for all the dimensions
    if typeof(rk) == Int64
       rk = rk * ones(Int64, order)
    end
    order == length(rk) || error("dimension of t does not match the length of rank")

    # initialization factors
    u = Vector{Matrix{Tv}}(undef, order)

    # initialize as random matrix
    if init_flag == "random"
       for i = 1 : order
           u[i] = randn(Tv, dims[i], rk[i])
       end

    # initialize as eigenvectors
    elseif init_flag == "eigvec"
       for i = 1 : order
           u[i] = nvecs(t, i, rk[i])
       end

    else
       error("unsupported initialization method")
    end

    fit = 0.0; fitold = 0.0;
    core= zeros(Tv, rk...)

    # iterations
    for iter = 1 : max_iter

        # save the fitting
        fitold = fit

        # save the last
        utilde = zeros()

        for i = 1 : order

            # update factor matrix
            v = setdiff((1:order), i)
            utilde = ttm(u, t, v; transpose_flag=true)

            # update ith factor matrix
            u[i]   = nvecs(utilde, i, rk[i])
        end

        # update core tensor
        core = ttm(u[order], utilde, order; transpose_flag=true)

        # estimate residue
        normresidual = sqrt(abs(normt^2 - norm(core.d)^2))
        fit = 1 - normresidual / normt
        fitchange = fit - fitold

        @printf("iter %4d %.16f %.16f\n", iter, fit, fitchange)

        # if fitchange < tol
        #    break
        # end
    end

    return Tucker(order, dims, rk, core.d, u)
end

# c  = randn(3,4,5);
# rk = collect(size(c));
# u  = Vector{Matrix{Float64}}(undef, 3)
# dims = [21, 22, 23]
# for i = 1 : 3
#     a = randn(dims[i], dims[i])
#     (p, s, v) = svd(a)
#     u[i] = p[:,1:rk[i]]
# end
# t = Tucker(c, u);
# d = tucker2tensor(t);
#
# dn = copy_tensor(d);
# dn.d .= dn.d + randn(dims...)*0.05;
# ratio = norm(d.d-dn.d) / norm(d.d)
#
# r = tucker_als(dn, rk; max_iter=10, init_flag="eigvec");
# d1= tucker2tensor(r);
# norm(d1.d-d.d) / norm(d.d)

# """
#    Wrap tucker decomposition, the input is multi-dimensional array
# """
# function WrapTuckerAls(dp::Array{Tv}, rk::Int64) where Tv <: AbstractFloat
#     N    = ndims(dp)
#
#     X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
#     T    = tuckerAls(X, rk, tol=tol)
#     dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
#     return dp
#     return nothing
# end
#
# """
#    Wrap tucker decomposition, the input is a tuple
# """
# function WrapTuckerAls(par::Tuple{String, Vector{Int64}, Float64})
#     path = par[1]
#     println("$path")
#     rk   = par[2]
#     tol  = par[3]
#     dp   = getOnePatch(path)
#     N    = ndims(dp)
#     X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
#     T    = tuckerAls(X, rk, tol=tol)
#     dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
#     fid  = open(path, "r+"); pos = sizeof(Int32)*9;
#     seek(fid, pos); write(fid, dp); close(fid);
#     return nothing
# end
