"""
   struct for tucker tensor
"""
struct Tucker{Tv<:Number}
     N    :: Int64
     rk   :: Vector{Int64}
     I    :: Vector{Int64}
     cory :: Array{Tv}
     U    :: Vector{Matrix{Tv}}
end

"""
   copy a tucker tensor
"""
function copy_tucker(T::Tucker)
    return Tucker(T.N, copy(T.rk), copy(T.I), copy(T.core), deepcopy(U))
end

"""
   create a tucker tensor with provided core tensor and factors in each dimension
"""
function init_tucker(core::Array{Tv}, U::Vector{Matrix{Tv}}) where {Tv<:Number}

    # check dimensions
    N = ndims(core)
    rk = collect(size(core))

    N == length(U) || error("core dimension does not match factors")
    I = zeros(Int64,N)
    for i = 1 : N
        rk[i] == size(U[i], 2) || error("core dimension does not match factors")
        I[i] = size(U[i],1)
    end

    return Tucker(N, rk, I, copy(core), deepcopy(U))
end

"""
   convert low-rank tensor to dense tensor
"""
function tucker2tensor(t::Tucker)

    C = init_tensor(t.core)
    for i = 1 : t.N
        C = ttm(C, T.U[i], i)
    end
    return C
end

"""
   compute the n^{th} mode factor matrix
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
function tucker_als(X::Tensor, rk::Union{Vector{Ti},Ti};
                   tol=1e-4, Niter=50, init_flag="rand") where {Ti<:Int64}

    data_format = eltype(X.d)
    N = X.N
    I = X.I
    fitchangetol = tol
    normX = norm(X.d)

    if typeof(rk) <: Integer
       rk = rk*ones(Int64, N)
    end
    N == length(rk) || error("dimension of D does not match length of rank")

    # initialization factors
    U = Vector{Matrix{data_format}}(undef, N)
    if init == "rand"
       for i = 1 : N
           U[i] = randn(data_format, I[i], rk[i])
       end
    elseif init == "eigvec"
       for i = 1 : N
           U[i] = nvecs(X, i, rk[i])
       end
    else
       error("unsupported initialization method")
    end

    fit = 0.0; fitold = 0.0;
    core = zeros(etp, rk...)
    for iter = 1 : maxIter
        fitold = fit
        Utilde = zeros()
        for n = 1 : N
            Utilde = ttm(X, U, -n, tflag=true)
            U[n]   = nvecs(Utilde, n, rk[n])
        end
        core = ttm(Utilde, U[N], N, tflag=true)
        normresidual = sqrt(abs(normX^2 - vecnorm(core.D)^2))
        fit = 1 - normresidual/normX
        fitchange = fit - fitold

        if fitchange < tol
           break
        end
        # println("iter $iter, fit $fit, fitchange $fitchange")
    end
    return tucker(N, I, rk, core.D, U)
end

"""
   Wrap tucker decomposition, the input is multi-dimensional array
"""
function WrapTuckerAls(dp::Array{Tv}, rk::Int64) where Tv <: AbstractFloat
    N    = ndims(dp)

    X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
    T    = tuckerAls(X, rk, tol=tol)
    dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
    return dp
    return nothing
end

"""
   Wrap tucker decomposition, the input is a tuple
"""
function WrapTuckerAls(par::Tuple{String, Vector{Int64}, Float64})
    path = par[1]
    println("$path")
    rk   = par[2]
    tol  = par[3]
    dp   = getOnePatch(path)
    N    = ndims(dp)
    X    = tensor(N, collect(size(dp)), convert(Array{Float64,N}, dp))
    T    = tuckerAls(X, rk, tol=tol)
    dp   = tucker2Tensor(T); dp = convert(Array{Float32,1}, vec(dp.D))
    fid  = open(path, "r+"); pos = sizeof(Int32)*9;
    seek(fid, pos); write(fid, dp); close(fid);
    return nothing
end
