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

# """
#    compute tensor times a series of vectors
# """
# function ttv(X::Tensor{Tv,Ti}, v::Vector{Vector{Tv}}, dims::Vector{Ti}) {Tv<:Number, Ti<:Int64}
#
#     for i = 1 : length(dims)
#         1 <= dims[i] <= X.N || throw(DimensionMismatch())
#         X.I[dims[i]] == length(v[i]) || throw(DimensionMismatch())
#     end
#
#     remdims = setdiff(collect(1 : X.N), dims)
#     pdims = vcat(remdims, dims)
#
#     newI  = X.I[pdims]
#     des = zeros(Tv, newI...)
#     if length(dims) > 1
#        permutedims!(des, X.D, pdims)
#     end
#
#     for i = length(dims) : -1 : 1
#         des = reshape(des, prod(newI[1:n-1]), newI[n])
#         des = des * v[i]
#         n = n - 1
#     end
#     return des
# end

"""
   tensor times one matrix, which M × Xn
"""
function ttm(M::Matrix{Tv}, T::Tensor{Tv,Ti}, n::Ti) where {Tv<:Number, Ti<:Int64}

    # check dimension
    (n1, n2) = size(M)
    n2 == T.I[n] || throw(DimensionMismatch())

    # the dimension of result tensor
    I = copy(T.I)
    I[n] = n1

    # matricization
    Tn = matricization(T, n)
    D  = M * Tn

    return unmatricization(D, I, n)
end

# """
#    tensor times a series of matrices
# """
# function ttm(U::Vector{Matrix{Tv}}, X::Tensor, n::Vector{Ti}) where {Tv<:Number, Ti<:Int64}
#     N = X.N
#     length(U) == N || error("number of factor matrix does not match")
#     if n < 0
#        n = setdiff(collect(1:N), -n)
#     end
#     maximum(n) <= N || error("too large dimension")
#     minimum(n) >= 1 || error("too small dimension")
#     D = tensor(X)
#     for i in n
#         D = ttm(D, U[i], i, tflag=tflag)
#     end
#     return D
# end
