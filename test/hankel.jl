@testset "multi-level Hankel matrix or its adjoint times a vector" begin

  for (Tv, etol) in [(Float32,1.0e-5), (Float64,1.0e-14), (Complex{Float32},1.0e-5), (Complex{Float64},1.0e-14)]

      # 1D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      d  = randn(Tv, N1)                # randnom vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1)                # randnom vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="forward")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1)                # randnom vector
      r1 = H' * vec(v)                  # matrix times vector
      r2 = hankel_multiplication(d, v; flag="adjoint"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="adjoint")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference


      # 2D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      d  = randn(Tv, N1, N2)            # randnom vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2)            # randnom vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="forward")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2)            # randnom vector
      r1 = H' * vec(v)                  # matrix times vector
      r2 = hankel_multiplication(d, v; flag="adjoint"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="adjoint")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference


      # 3D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 10; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      d  = randn(Tv, N1, N2, N3)        # randnom vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3)        # randnom vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="forward")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3)        # randnom vector
      r1 = H' * vec(v)                  # matrix times vector
      r2 = hankel_multiplication(d, v; flag="adjoint"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="adjoint")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference


      # 4D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector
      d  = randn(Tv, N1, N2, N3, N4)     # randnom vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3, K4)    # randnom vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="forward")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3, L4)    # randnom vector
      r1 = H' * vec(v)                  # matrix times vector
      r2 = hankel_multiplication(d, v; flag="adjoint"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="adjoint")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # 5D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector
      L5 = 4; K5 = L5-1; N5 = L5 + K5 -1 # length of vector
      d  = randn(Tv, N1, N2, N3, N4, N5) # randnom vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3, K4, K5)# randnom vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="forward")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3, L4, L5)# randnom vector
      r1 = H' * vec(v)                  # matrix times vector
      r2 = hankel_multiplication(d, v; flag="adjoint"); r2 = vec(r2);  # efficient hankel multiplication
      r3 = hankel_multiplication(d, vec(v); flag="adjoint")
      @test norm(r1-r2)/norm(r1) < etol # test the difference
      @test norm(r1-r3)/norm(r3) < etol # test the difference
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
