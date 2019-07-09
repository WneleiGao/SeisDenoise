@testset "multi-level Hankel matrix or its adjoint times a vector" begin

  for (Tv, etol) in [(Float32,1.0e-5), (Float64,1.0e-14), (Complex{Float32},1.0e-5), (Complex{Float64},1.0e-14)]

      # 1D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector

      d  = randn(Tv, N1)
      v  = randn(Tv, K1)

      T  = build_toeplitz_matrix(d)
      H  = build_hankel_matrix(d)

      r1 = T * vec(v)
      v_hat = reverse_order(v)
      r2    = H * vec(v_hat)
      @test norm(r1-r2)/norm(r1) < etol # test the difference


      # 2D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector

      d  = randn(Tv, N1, N2)
      v  = randn(Tv, K1, K2)

      T  = build_toeplitz_matrix(d)
      H  = build_hankel_matrix(d)

      r1 = T * vec(v)
      v_hat = reverse_order(v)
      r2    = H * vec(v_hat)
      @test norm(r1-r2)/norm(r1) < etol # test the difference


      # 3D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 10; K3 = L3-1; N3 = L3 + K3 -1 # length of vector

      d  = randn(Tv, N1, N2, N3)
      v  = randn(Tv, K1, K2, K3)

      T  = build_toeplitz_matrix(d)
      H  = build_hankel_matrix(d)

      r1 = T * vec(v)
      v_hat = reverse_order(v)
      r2    = H * vec(v_hat)
      @test norm(r1-r2)/norm(r1) < etol # test the difference


      # 4D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector

      d  = randn(Tv, N1, N2, N3, N4)
      v  = randn(Tv, K1, K2, K3, K4)

      T  = build_toeplitz_matrix(d)
      H  = build_hankel_matrix(d)

      r1 = T * vec(v)
      v_hat = reverse_order(v)
      r2    = H * vec(v_hat)
      @test norm(r1-r2)/norm(r1) < etol # test the difference


      # 5D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector
      L5 = 4; K5 = L5-1; N5 = L5 + K5 -1 # length of vector

      d  = randn(Tv, N1, N2, N3, N4, N5)
      v  = randn(Tv, K1, K2, K3, K4, K5)

      T  = build_toeplitz_matrix(d)
      H  = build_hankel_matrix(d)

      r1 = T * vec(v)
      v_hat = reverse_order(v)
      r2    = H * vec(v_hat)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

  end
end


"""
   slow method to count the repeated copy time when building Hankel matrix from a multi-dimensional array
"""
function count_copy_times_slow(dims::Union{Ti,Vector{Ti}}) where {Ti<:Int64}

    # determine the dimensions of an array
    N = length(dims)

    # level 1
    if N == 1

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1


       # build hankel
       for j1 = 1 : K1
           for i1 = 1 : L1
               n1 = i1+j1-1

               count_num[n1] += 1
           end
       end


    # level 2
    elseif N == 2

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1

       # build Hankel
       for j2 = 1 : K2
           c2 = (j2-1)*K1
           for i2 = 1 : L2
               r2 = (i2-1)*L1
               n2 = i2 + j2 - 1

               # first layer
               for j1 = 1 : K1
                   c1 = c2 + j1
                   for i1 = 1 : L1
                       r1 = r2 + i1
                       n1 = i1 + j1 - 1

                       count_num[n1,n2] += 1
                   end
               end #first
           end
       end #second


    # level 3
    elseif N == 3

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1

       # build Hankel
       for j3 = 1 : K3
           c3 = (j3-1)*K2*K1
           for i3 = 1 : L3
               r3 = (i3-1)*L2*L1
               n3 = i3 + j3 - 1

               # second layer
               for j2 = 1 : K2
                   c2 = (j2-1)*K1
                   for i2 = 1 : L2
                       r2 = (i2-1)*L1
                       n2 = i2 + j2 - 1

                       # first layer
                       for j1 = 1 : K1
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L1
                               r1 = r3 + r2 + i1
                               n1 = i1 + j1 - 1

                               count_num[n1,n2,n3] += 1
                           end
                       end #first
                   end
               end #second
           end
       end #third


    # level 4
    elseif N == 4

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1

       # build Hankel
       for j4 = 1 : K4
           c4 = (j4-1)*K3*K2*K1
           for i4 = 1 : L4
               r4 = (i4-1)*L3*L2*L1
               n4 = i4 + j4 - 1

               # third layer
               for j3 = 1 : K3
                   c3 = (j3-1)*K2*K1
                   for i3 = 1 : L3
                       r3 = (i3-1)*L2*L1
                       n3 = i3 + j3 - 1

                       # second layer
                       for j2 = 1 : K2
                           c2 = (j2-1)*K1
                           for i2 = 1 : L2
                               r2 = (i2-1)*L1
                               n2 = i2 + j2 - 1

                               # first layer
                               for j1 = 1 : K1
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L1
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = i1 + j1 - 1

                                       count_num[n1,n2,n3,n4] += 1
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    # level 5
    elseif N == 5

       # allocate memory to record the number
       count_num = zeros(Int64, dims[1], dims[2], dims[3], dims[4], dims[5])

       # size of hankel matrix
       L1 = floor(Int64,dims[1]/2)+1; K1 = dims[1]-L1+1
       L2 = floor(Int64,dims[2]/2)+1; K2 = dims[2]-L2+1
       L3 = floor(Int64,dims[3]/2)+1; K3 = dims[3]-L3+1
       L4 = floor(Int64,dims[4]/2)+1; K4 = dims[4]-L4+1
       L5 = floor(Int64,dims[5]/2)+1; K5 = dims[5]-L5+1

       # build Hankel
       for j5 = 1 : K5
           c5 = (j5-1)*K4*K3*K2*K1
           for i5 = 1 : L5
               r5 = (i5-1)*L4*L3*L2*L1
               n5 = i5 + j5 - 1

               for j4 = 1 : K4
                   c4 = (j4-1)*K3*K2*K1
                   for i4 = 1 : L4
                       r4 = (i4-1)*L3*L2*L1
                       n4 = i4 + j4 - 1

                       # third layer
                       for j3 = 1 : K3
                           c3 = (j3-1)*K2*K1
                           for i3 = 1 : L3
                               r3 = (i3-1)*L2*L1
                               n3 = i3 + j3 - 1

                               # second layer
                               for j2 = 1 : K2
                                   c2 = (j2-1)*K1
                                   for i2 = 1 : L2
                                       r2 = (i2-1)*L1
                                       n2 = i2 + j2 - 1

                                       # first layer
                                       for j1 = 1 : K1
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L1
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = i1 + j1 - 1

                                               count_num[n1,n2,n3,n4,n5] += 1
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up to 5 dimensions")
    end

    return count_num
end

@testset "count the copy times when build Hankel matrix" begin

    etol = 1e-15

    # 1D case
    dims = 22
    count_num1 = count_copy_times(dims)
    count_num2 = count_copy_times_slow(dims)
    @test norm(count_num1-count_num2) / norm(count_num2) < etol

    # 2D case
    dims = [22,21]
    count_num1 = count_copy_times(dims)
    count_num2 = count_copy_times_slow(dims)
    @test norm(count_num1-count_num2) / norm(count_num2) < etol

    # 3D case
    dims = [10,9,8]
    count_num1 = count_copy_times(dims)
    count_num2 = count_copy_times_slow(dims)
    @test norm(count_num1-count_num2) / norm(count_num2) < etol

    # 4D case
    dims = [10,9,8,7]
    count_num1 = count_copy_times(dims)
    count_num2 = count_copy_times_slow(dims)
    @test norm(count_num1-count_num2) / norm(count_num2) < etol

    # 5D case
    dims = [10,9,8,7,6]
    count_num1 = count_copy_times(dims)
    count_num2 = count_copy_times_slow(dims)
    @test norm(count_num1-count_num2) / norm(count_num2) < etol
end


"""
   slow anti-diagonal summation via building multi-level hankel matrix first.
"""
function anti_diagonal_summation_slow(u::Vector{Tv}, v::Vector{Tv},
                                      L::Union{Ti,Vector{Ti}},
                                      K::Union{Ti,Vector{Ti}}) where {Tv <: Number, Ti<:Int64}

    # order of multi-level hankel matrix
    order =  length(L)
    order == length(K) || error(DimensionMismatch("length of L mismatch length of K"))

    # outer product u * v'
    H = u * v'

    if order == 1

       # check the size of input
       L[1] == length(u) || error(DimensionMismatch("check the length of u"))
       K[1] == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1

       # allocate memory for the result
       d  = zeros(Tv,N1)

       # first layer
       for j1 = 1 : K[1]
           for i1 = 1 : L[1]
               n1 = i1 + j1 - 1

               d[n1] += H[i1,j1]
           end
       end #first


    elseif order == 2

       # check the size of input
       prod(L[1:2]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:2]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2)

       # second layer
       for j2 = 1 : K[2]
           c2 = (j2-1)*K[1]
           for i2 = 1 : L[2]
               r2 = (i2-1)*L[1]
               n2 = i2 + j2 - 1

               # first layer
               for j1 = 1 : K[1]
                   c1 = c2 + j1
                   for i1 = 1 : L[1]
                       r1 = r2 + i1
                       n1 = i1 + j1 - 1

                       d[n1,n2] += H[r1,c1]
                   end
               end #first
           end
       end #second


    elseif order == 3

       # check the size of input
       prod(L[1:3]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:3]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3)

       # anti-diagonal summation
       for j3 = 1 : K[3]
           c3 = (j3-1)*K[2]*K[1]
           for i3 = 1 : L[3]
               r3 = (i3-1)*L[2]*L[1]
               n3 = i3 + j3 - 1

               # second layer
               for j2 = 1 : K[2]
                   c2 = (j2-1)*K[1]
                   for i2 = 1 : L[2]
                       r2 = (i2-1)*L[1]
                       n2 = i2 + j2 - 1

                       # first layer
                       for j1 = 1 : K[1]
                           c1 = c3 + c2 + j1
                           for i1 = 1 : L[1]
                               r1 = r3 + r2 + i1
                               n1 = i1 + j1 - 1

                               d[n1,n2,n3] += H[r1,c1]
                           end
                       end #first
                   end
               end #second
           end
       end #thirds


    elseif order == 4

       # check the size of input
       prod(L[1:4]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:4]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4)

       # anti-diagonal summation
       for j4 = 1 : K[4]
           c4 = (j4-1)*K[3]*K[2]*K[1]
           for i4 = 1 : L[4]
               r4 = (i4-1)*L[3]*L[2]*L[1]
               n4 = i4 + j4 - 1

               # third layer
               for j3 = 1 : K[3]
                   c3 = (j3-1)*K[2]*K[1]
                   for i3 = 1 : L[3]
                       r3 = (i3-1)*L[2]*L[1]
                       n3 = i3 + j3 - 1

                       # second layer
                       for j2 = 1 : K[2]
                           c2 = (j2-1)*K[1]
                           for i2 = 1 : L[2]
                               r2 = (i2-1)*L[1]
                               n2 = i2 + j2 - 1

                               # first layer
                               for j1 = 1 : K[1]
                                   c1 = c4 + c3 + c2 + j1
                                   for i1 = 1 : L[1]
                                       r1 = r4 + r3 + r2 + i1
                                       n1 = i1 + j1 - 1

                                       d[n1,n2,n3,n4] += H[r1,c1]
                                   end
                               end #first
                           end
                       end #second
                   end
               end #third
           end
       end #fourth


    elseif order == 5

       # check the size of input
       prod(L[1:5]) == length(u) || error(DimensionMismatch("check the length of u"))
       prod(K[1:5]) == length(v) || error(DimensionMismatch("check the length of v"))

       # the size of array
       N1 = L[1]+K[1]-1
       N2 = L[2]+K[2]-1
       N3 = L[3]+K[3]-1
       N4 = L[4]+K[4]-1
       N5 = L[5]+K[5]-1

       # allocate memory for the result
       d  = zeros(Tv, N1, N2, N3, N4, N5)

       # anti-diagonal summation
       for j5 = 1 : K[5]
           c5 = (j5-1)*K[4]*K[3]*K[2]*K[1]
           for i5 = 1 : L[5]
               r5 = (i5-1)*L[4]*L[3]*L[2]*L[1]
               n5 = i5 + j5 - 1

               for j4 = 1 : K[4]
                   c4 = (j4-1)*K[3]*K[2]*K[1]
                   for i4 = 1 : L[4]
                       r4 = (i4-1)*L[3]*L[2]*L[1]
                       n4 = i4 + j4 - 1

                       # third layer
                       for j3 = 1 : K[3]
                           c3 = (j3-1)*K[2]*K[1]
                           for i3 = 1 : L[3]
                               r3 = (i3-1)*L[2]*L[1]
                               n3 = i3 + j3 - 1

                               # second layer
                               for j2 = 1 : K[2]
                                   c2 = (j2-1)*K[1]
                                   for i2 = 1 : L[2]
                                       r2 = (i2-1)*L[1]
                                       n2 = i2 + j2 - 1

                                       # first layer
                                       for j1 = 1 : K[1]
                                           c1 = c5 + c4 + c3 + c2 + j1
                                           for i1 = 1 : L[1]
                                               r1 = r5 + r4 + r3 + r2 + i1
                                               n1 = i1 + j1 - 1

                                               d[n1,n2,n3,n4,n5] += H[r1,c1]
                                           end
                                       end #first
                                   end
                               end #second
                           end
                       end #third
                   end
               end #fourth
           end
       end # fiveth


    else
       error("only support up-to fiveth order")
    end

    return d
end


@testset "anti-diagonal summation of Hankel matrix" begin

  for (Tv, etol) in [(Float32,1.0e-5), (Float64,1.0e-14), (Complex{Float32},1.0e-5), (Complex{Float64},1.0e-14)]

      # 1D case
      L = [12]
      K = L .- 1
      u = randn(Tv, prod(L[1:1]))
      v = randn(Tv, prod(K[1:1]))
      d1 = anti_diagonal_summation(u, v, L, K);
      d2 = anti_diagonal_summation_slow(u, v, L, K);
      @test norm(d1-d2) / norm(d1) < etol

      # 2D case
      L = [12,11]
      K = L .- 1
      u = randn(Tv, prod(L[1:2]))
      v = randn(Tv, prod(K[1:2]))
      d1 = anti_diagonal_summation(u, v, L, K);
      d2 = anti_diagonal_summation_slow(u, v, L, K);
      @test norm(d1-d2) / norm(d1) < etol

      # 3D case
      L = [12,11,10]
      K = L .- 1
      u = randn(Tv, prod(L[1:3]))
      v = randn(Tv, prod(K[1:3]))
      d1 = anti_diagonal_summation(u, v, L, K);
      d2 = anti_diagonal_summation_slow(u, v, L, K);
      @test norm(d1-d2) / norm(d1) < etol

      # 4D case
      L = [9,8,7,6]
      K = L .- 1
      u = randn(Tv, prod(L[1:4]))
      v = randn(Tv, prod(K[1:4]))
      d1 = anti_diagonal_summation(u, v, L, K);
      d2 = anti_diagonal_summation_slow(u, v, L, K);
      @test norm(d1-d2) / norm(d1) < etol

      # 5D case
      L = [8,7,6,5,4]
      K = L .- 1
      u = randn(Tv, prod(L[1:5]))
      v = randn(Tv, prod(K[1:5]))
      d1 = anti_diagonal_summation(u, v, L, K);
      d2 = anti_diagonal_summation_slow(u, v, L, K);
      @test norm(d1-d2) / norm(d1) < etol
  end
end
