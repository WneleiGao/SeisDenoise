@testset "multi-level Hankel matrix or its adjoint times a vector" begin

  for (Tv, etol) in [(Float32,1.0e-5), (Float64,1.0e-14), (Complex{Float32},1.0e-5), (Complex{Float64},1.0e-14)]

      # 1D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      d  = randn(Tv, N1)                # random vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1)                # random vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward")  # efficient hankel multiplication
      r2 = vec(r2)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1)                # random vector
      r3 = H' * vec(v)                  # matrix times vector
      r4 = hankel_multiplication(d, v; flag="adjoint")  # efficient hankel multiplication
      r4 = vec(r4)
      @test norm(r3-r4)/norm(r4) < etol # test the difference


      # 2D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      d  = randn(Tv, N1, N2)            # random vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2)            # random vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward")  # efficient hankel multiplication
      r2 = vec(r2)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2)            # random vector
      r3 = H' * vec(v)                  # matrix times vector
      r4 = hankel_multiplication(d, v; flag="adjoint")  # efficient hankel multiplication
      r4 = vec(r4)
      @test norm(r3-r4)/norm(r4) < etol # test the difference


      # 3D case
      L1 = 12; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 11; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 10; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      d  = randn(Tv, N1, N2, N3)        # random vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3)        # random vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward")  # efficient hankel multiplication
      r2 = vec(r2)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3)        # random vector
      r3 = H' * vec(v)                  # matrix times vector
      r4 = hankel_multiplication(d, v; flag="adjoint")  # efficient hankel multiplication
      r4 = vec(r4)
      @test norm(r3-r4)/norm(r4) < etol # test the difference


      # 4D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector
      d  = randn(Tv, N1, N2, N3, N4)     # random vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3, K4)    # random vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward")  # efficient hankel multiplication
      r2 = vec(r2)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3, L4)    # random vector
      r3 = H' * vec(v)                  # matrix times vector
      r4 = hankel_multiplication(d, v; flag="adjoint")  # efficient hankel multiplication
      r4 = vec(r4)
      @test norm(r3-r4)/norm(r4) < etol # test the difference


      # 5D case
      L1 = 6; K1 = L1-1; N1 = L1 + K1 -1 # length of vector
      L2 = 5; K2 = L2-1; N2 = L2 + K2 -1 # length of vector
      L3 = 4; K3 = L3-1; N3 = L3 + K3 -1 # length of vector
      L4 = 4; K4 = L4-1; N4 = L4 + K4 -1 # length of vector
      L5 = 4; K5 = L5-1; N5 = L5 + K5 -1 # length of vector
      d  = randn(Tv, N1, N2, N3, N4, N5) # random vector

      H  = build_hankel_matrix(d)       # build hankel matrix

      # forward
      v  = randn(Tv, K1, K2, K3, K4, K5)# random vector
      r1 = H * vec(v)                   # matrix times vector
      r2 = hankel_multiplication(d, v; flag="forward")  # efficient hankel multiplication
      r2 = vec(r2)
      @test norm(r1-r2)/norm(r1) < etol # test the difference

      # adjoint
      v  = randn(Tv, L1, L2, L3, L4, L5)# random vector
      r3 = H' * vec(v)                  # matrix times vector
      r4 = hankel_multiplication(d, v; flag="adjoint")  # efficient hankel multiplication
      r4 = vec(r4)
      @test norm(r3-r4)/norm(r4) < etol # test the difference
  end

end
