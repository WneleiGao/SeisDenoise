"""
   tensor structure
"""
struct Tensor{Tv<:Number, Ti<:Int64}
     N :: Ti
     I :: Vector{Ti}
     D :: Array{Tv}
end

"""
   create a tensor object when the input is a multi-dimensional array
"""
function Tensor(D::Array{Tv}) where {Tv <: Number}
    return Tensor(ndims(D), collect(size(D)), copy(D))
end

"""
   copy a tucker tensor
"""
function copy_tensor(T::Tensor)
    return Tensor(T.N, copy(T.I), copy(T.D))
end

"""
   unfolding tensor along n dimension
"""
function matricization(X::Tensor, n::Int64)

    N = X.N
    I = X.I

    if n == 1 # first dimension
       Xn = copy(reshape(X.D, I[1], prod(I[2 : N])))
    else
       p = collect(1:N); p[n] = 1; p[1] = n;  # permute first and nth dimension
       Xn = permutedims(X.D, p)
       Xn = reshape(Xn, I[n], round(Int64, prod(I)/I[n]))
    end

    return Xn
end

"""
   the input is multi-dimensional array
"""
function matricization(X::Array{Tv}, n::Int64) where Tv<:Number

    I  = collect(size(X))
    N  = ndims(X)

    if n<1 || n > N
       error("n is out of bound")
    end

    if n == 1 # first dimension
       Xn = copy(reshape(X, I[1], prod(I[2:N])))
    else
       p = collect(1:N); p[n] = 1; p[1] = n;  # permute first and nth dimension
       Xn = permutedims(X, p)
       Xn = reshape(Xn, I[n], round(Int64, prod(I)/I[n]))
    end
    return Xn
end

"""
folding an matrix back to a tensor
"""
function unmatricization(D::Matrix{Tv}, I::Vector{Ti}, n::Ti; return_flag="tensor") where {Tv <: Number, Ti<:Int64}

    if n == 1
       D = reshape(D, I...)
    else
       I1 = copy(I); I1[1] = I[n]; I1[n] = I[1];
       D = reshape(D, I1...)

       p = collect(1:length(I)); p[n] = 1; p[1] = n;
       D = permutedims(D, p)
    end

    if return_flag == "tensor"
       return Tensor(D)
    elseif return_flag == "array"
       return D
    end
end

"""
   Khatri-Rao product a2 ⊙ a1
"""
function kr_product(a2::Vector{Tv}, a1::Vector{Tv}) where {Tv<:Number}

    m1 = length(a1)
    m2 = length(a2)
    m  = m1*m2

    r  = Vector{Tv}(undef, m)

    for i2 = 1 : m2
        idx = (i2-1) * m1

        for i1 = 1 : m1
            j = idx + i1
            r[j] = a2[i2] * a1[i1]
        end
    end
    return r
end

"""
   Khatri-Rao product for two matrices
"""
function kr_product(a2::Matrix{Tv}, a1::Matrix{Tv}) where {Tv<:Number}

    # check dimension
    size(a1,2) == size(a2,2) || throw(DimensionMismatch("size(a1,2) !=  size(a2,2)"))

    m1 = size(a1,1)
    m2 = size(a2,1)

    m  = m1*m2
    n  = size(a1,2)

    r = Matrix{Tv}(undef, m, n)
    for j = 1 : n

        for i2 = 1 : m2
            idx= (i2-1)*m1

            for i1 = 1 : m1
                k = idx + i1
                r[k,j] = a2[i2,j] * a1[i1,j]
            end
        end
    end
    return r
end

"""
  recursive Khatri-Rao product A[N]⊙...⊙A[1]
"""
function recursive_krp(A::Vector{Vector{Tv}}) where {Tv<:Number}

    # number of vectors
    N = length(A)

    if N == 1
       return A[1]
    elseif N == 2
       return kr_product(A[2], A[1])
    else
       r = kr_product(A[2], A[1])
       for m = 3 : N
           r = kr_product(A[m], r)
       end
       return r
    end
end

"""
  recursive Khatri-Rao product A[N]⊙...⊙A[1]
"""
function recursive_krp(A::Vector{Matrix{Tv}}) where {Tv<:Number}

    # number of factors
    N = length(A)

    if N == 1
       return A[1]

    elseif N == 2
       return kr_product(A[2], A[1])

    else
       r = kr_product(A[2], A[1])
       for m = 3 : N
           r = kr_product(A[m], r)
       end
       return r

    end
end

"""
tensor times one matrix
"""
function ttm(U::Matrix{Tv}, X::Tensor{Tv,Ti}, n::Ti;
         transpose_flag=false) where {Tv <: Number, Ti<:Int64}

    # check dimension
    (n1, n2) = size(U)
    Xn = matricization(X, n)

    # copy the dimensions of tensor
    I = copy(X.I)

    if transpose_flag
       D = U' * Xn
       I[n] = n2

    else
       D = U  * Xn
       I[n] = n1
    end

    D = unmatricization(D, I, n)
    return D
end

"""
   tensor times a series of matrices
"""
function ttm(U::Vector{Matrix{Tv}}, X::Tensor{Tv,Ti}, n::Vector{Ti}; transpose_flag=false) where {Tv<:Number, Ti<:Int64}

    # check dimensions
    length(U)  == X.N || error("number of factor matrix does not match")
    maximum(n) <= X.N || error("too large dimension")
    minimum(n) >= 1   || error("too small dimension")

    D = copy_tensor(X)
    for i in n
        D = ttm(U[i], D, i; transpose_flag=transpose_flag)
    end
    return D
end
