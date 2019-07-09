@testset "circulant matrix times a vector" begin

  for (Tv, etol) in [(Float32,1.0e-5), (Float64,1.0e-14), (Complex{Float32},1.0e-5), (Complex{Float64},1.0e-14)]

      # 1D case
      N1 = 123            # length of vector
      d  = randn(Tv, N1)  # random vector
      v  = randn(Tv, N1)  # random vector

      C  = build_circulant_matrix(d)  # build circulant matrix

      r1 = C * vec(v)                 # matrix times vector
      r2 = ifft(fft(d) .* fft(v))     # fast method
      r2 = vec(r2)

      @test norm(r1-r2)/norm(r1) < etol  # test the difference


      # 2D case
      N1 = 13
      N2 = 14

      d  = randn(Tv, N1, N2)
      v  = randn(Tv, N1, N2)

      C  = build_circulant_matrix(d)

      r1 = C * vec(v)
      r2 = ifft(fft(d) .* fft(v))
      r2 = vec(r2)

      @test norm(r1-r2)/norm(r1) < etol


      # 3D case
      N1 = 6
      N2 = 7
      N3 = 8

      d  = randn(Tv, N1, N2, N3);
      v  = randn(Tv, N1, N2, N3);

      C  = build_circulant_matrix(d);

      r1 = C * vec(v)
      r2 = ifft(fft(d) .* fft(v))
      r2 = vec(r2)

      @test norm(r1-r2)/norm(r1) < etol


      # 4D case
      N1 = 3
      N2 = 4
      N3 = 5
      N4 = 6

      d  = randn(Tv, N1, N2, N3, N4)
      v  = randn(Tv, N1, N2, N3, N4)

      C  = build_circulant_matrix(d)

      r1 = C * vec(v)
      r2 = ifft(fft(d) .* fft(v))
      r2 = vec(r2)

      @test norm(r1-r2)/norm(r1) < etol


      # 5D case
      N1 = 3
      N2 = 4
      N3 = 5
      N4 = 4
      N5 = 3

      d  = randn(Tv, N1, N2, N3, N4, N5)
      v  = randn(Tv, N1, N2, N3, N4, N5)

      C  = build_circulant_matrix(d)

      r1 = C * vec(v)
      r2 = ifft(fft(d) .* fft(v))
      r2 = vec(r2)
      norm(r1-r2) / norm(r1)

      @test norm(r1-r2)/norm(r1) < etol
  end

end
