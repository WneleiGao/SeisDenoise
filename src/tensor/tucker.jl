"""
   struct for tucker tensor
"""
struct Tucker{Tv<:Number, Ti<:Int64}
     N    :: Ti
     rk   :: Vector{Ti}
     I    :: Vector{Ti}
     core :: Array{Tv}
     U    :: Vector{Matrix{Tv}}
end

"""
   constructor of tucker tensor with provided core tensor and factor matrices
"""
function Tucker(core::Array{Tv}, U::Vector{Matrix{Tv}}) where {Tv<:Number}

    # check dimensions
    N  = ndims(core)
    rk = collect(size(core))

    N == length(U) || error("core dimension does not match number of factors")
    I = zeros(Int64,N)

    for i = 1 : N
        rk[i] == size(U[i], 2) || error("core dimension does not match factors")
        I[i]  =  size(U[i], 1)
    end

    return Tucker(N, rk, I, copy(core), deepcopy(U))
end

"""
   copy a tucker tensor
"""
function copy_tucker(T::Tucker)
    return Tucker(T.N, copy(T.rk), copy(T.I), copy(T.core), deepcopy(T.U))
end

"""
   convert low-rank tucker tensor to dense tensor
"""
function tucker2tensor(t::Tucker)

    C = Tensor(t.core)
    for i = 1 : t.N
        C = ttm(T.U[i], C, i)
    end
    return C
end

"""
   compute the n^{th} mode factor matrix by keeping the first rk eigenvectors
"""
function nvecs(X::Tensor, n::Ti, rk::Ti) where {Ti<:Int64}

    Xn = matricization(X, n)
    Y  = Xn * Xn'
    (u, s, v) = svd(Y)
    return u[:,1:rk]
end

"""
   tucker decomposition by alternative minimization
"""
function tucker_als(X::Tensor{Tv,Ti}, rk::Union{Vector{Ti},Ti};
                   tol=1e-6, Niter=10, init_flag="eigvec") where {Tv<:Number, Ti<:Int64}

    N = X.N
    I = X.I
    fitchangetol = tol
    normX = norm(X.D)

    # one rank for all the dimensions
    if typeof(rk) <: Int64
       rk = rk*ones(Int64, N)
    end
    N == length(rk) || error("dimension of D does not match length of rank")

    # initialization factors
    U = Vector{Matrix{Tv}}(undef, N)

    # initialize as random matrix
    if init_flag == "rand"
       for i = 1 : N
           U[i] = randn(Tv, I[i], rk[i])
       end

    # initialize as eigenvectors
    elseif init_flag == "eigvec"
       for i = 1 : N
           U[i] = nvecs(X, i, rk[i])
       end

    else
       error("unsupported initialization method")
    end

    fit = 0.0; fitold = 0.0;
    core = zeros(Tv, rk...)

    # iterations
    for iter = 1 : Niter
        fitold = fit
        Utilde = zeros()

        for n = 1 : N

            # update factor matrix
            sets   = setdiff((1:N), n)
            Utilde = ttm(U, X, set; transpose_flag=true)

            # factor matrix obtained by svd
            U[n]   = nvecs(Utilde, n, rk[n])
        end

        # update core tensor
        core = ttm(U[N], Utilde, N; transpose_flag=true)

        # estimate residue
        normresidual = sqrt(abs(normX^2 - norm(core.D)^2))
        fit = 1 - normresidual/normX
        fitchange = fit - fitold

        println("iter $iter, fit $fit, fitchange $fitchange")

        # if fitchange < tol
        #    break
        # end
    end
    return Tucker(N, I, rk, core.D, U)
end

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
