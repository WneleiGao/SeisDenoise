@testset "multi-level teoplitz matrix or its adjoint times a vector" begin

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
