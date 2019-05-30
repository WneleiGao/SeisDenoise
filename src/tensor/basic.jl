"""
   tensor structure
"""
struct Tensor{Tv<:Number, Ti<:Int64}
     order :: Ti
     dims  :: Vector{Ti}
     d     :: Array{Tv}
end

"""
   create a tensor object when the input is a multi-dimensional array
"""
function Tensor(d::Array{Tv}) where {Tv <: Number}
    return Tensor(ndims(d), collect(size(d)), copy(d))
end

"""
   copy a tucker tensor
"""
function copy_tensor(t::Tensor)
    return Tensor(t.order, copy(t.dims), copy(t.d))
end

"""
   unfolding tensor along n dimension
"""
function matricization(t::Tensor, n::Int64)

    if n == 1 # first dimension
       nr = prod(t.dims[2:t.order])
       tn = copy(reshape(t.d, t.dims[1], nr))

    else
       p = collect(1 : t.order)
       p[1] = n
       p[n] = 1
       tn = permutedims(t.d, p)

       nr = round(Int64, prod(t.dims) / t.dims[n])
       tn = reshape(tn, t.dims[n], nr)
    end

    return tn
end

"""
   folding an matrix back to a tensor, dims is the dimension of the original tensor
"""
function unmatricization(d::Matrix{Tv}, dims::Vector{Ti}, n::Ti) where {Tv <: Number, Ti<:Int64}

    if n == 1
       d = reshape(d, dims...)

    else
       tmp = copy(dims)
       tmp[1] = dims[n]
       tmp[n] = dims[1]
       d = reshape(d, tmp...)

       p = collect(1:length(dims))
       p[1] = n
       p[n] = 1
       d = permutedims(d, p)
    end

    return Tensor(d)
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
function ttm(u::Matrix{Tv}, t::Tensor{Tv,Ti}, n::Ti;
         transpose_flag=false) where {Tv <: Number, Ti<:Int64}

    (k, l) = size(u)

    # check dimension
    if transpose_flag
       k == t.dims[n] || error("size mismatch")
    else
       l == t.dims[n] || error("size mismatch")
    end

    tn = matricization(t, n)

    # copy the dimensions of tensor
    tmp = copy(t.dims)

    if transpose_flag
       d = u' * tn
       tmp[n] = l

    else
       d = u  * tn
       tmp[n] = k
    end

    return unmatricization(d, tmp, n)
end

"""
   tensor times a series of matrices, designed specially for HOSVD
"""
function ttm(u::Vector{Matrix{Tv}}, t::Tensor{Tv,Ti}, n::Vector{Ti}; transpose_flag=false) where {Tv<:Number, Ti<:Int64}

    # check dimensions
    length(u)  == t.order || error("number of factor matrix does not match")
    maximum(n) <= t.order || error("too large dimension")
    minimum(n) >= 1       || error("too small dimension")

    d = copy_tensor(t)
    for i in n
        d = ttm(u[i], d, i; transpose_flag=transpose_flag)
    end

    return d
end

# tensor
# d = rand(11, 12, 13);
# t = Tensor(d);
#
# u = Vector{Matrix{Float64}}(undef, 3)
# for i = 1 : 3
#     a = randn(100, 100);
#     (p, s, v) = svd(a);
#     u[i] = p[:,1 : t.dims[i]]
# end
# r = ttm(u, t, [3,1,2])
# c = ttm(u, r, [2,3,1]; transpose_flag=true)
