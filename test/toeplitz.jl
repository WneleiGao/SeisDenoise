# using FFTW, LinearAlgebra
# #==============================================================================#
# #                                   1D case
# #==============================================================================#
# Tv = Complex{Float64};
# # Tv = Float32;
#
# L1 = 123; K1 = L1-1;
# # L1 = 5; K1 = L1;
# N1 = L1 + K1 - 1;
#
#
# d  = randn(Tv, N1);
# v  = randn(Tv, K1);
#
# H  = build_hankel_matrix(d);
# T  = build_toeplitz_matrix(d);
#
# r1 = H * vec(v);
# v_hat = reverse_order(v);
# r2 = T * vec(v_hat);
# norm(r1-r2) / norm(r1)
#
#
# # forward
# r1 = H * vec(v);
# r2 = hankel_multiplication(d, v; flag="forward");
# r2 = vec(r2);
# norm(r1-r2) / norm(r1)
#
# # adjoint
# v = randn(Tv, L1);
# r3 = H' * vec(v);
# r4 = hankel_multiplication(d, v; flag="adjoint");
# r4 = vec(r4);
# norm(r3-r4) / norm(r3)
#
#
#
# #==============================================================================#
# #                                   2D case
# #==============================================================================#
# Tv = Complex{Float64};
# # Tv = Float64;
#
# L1 = 5; K1 = L1-1;
# L2 = 6; K2 = L2-1;
#
# # L1 = 5; K1 = L1;
# # L2 = 6; K2 = L2;
#
# N1 = L1 + K1 - 1;
# N2 = L2 + K2 - 1;
#
# d  = randn(Tv, N1, N2);
# v  = randn(Tv, K1, K2);
#
# H  = build_hankel_matrix(d);
# T  = build_toeplitz_matrix(d);
#
# r1 = H * vec(v);
# v_hat = reverse_order(v);
# r2 = T * vec(v_hat);
# norm(r1-r2) / norm(r1)
#
#
# # forward
# r1 = H * vec(v);
# r2 = hankel_multiplication(d, v; flag="forward");
# r2 = vec(r2);
# norm(r1-r2) / norm(r1)
#
# # adjoint
# v = randn(Tv, L1, L2);
# r3 = H' * vec(v);
# r4 = hankel_multiplication(d, v; flag="adjoint");
# r4 = vec(r4);
# norm(r3-r4) / norm(r3)
#
#
# #==============================================================================#
# #                                   3D case
# #==============================================================================#
# # Tv = Complex{Float64};
# Tv = Float64;
#
# L1 = 5; K1 = L1-1;
# L2 = 6; K2 = L2-1;
# L3 = 7; K3 = L3-1;
#
# # L1 = 5; K1 = L1;
# # L2 = 6; K2 = L2;
# # L3 = 7; K3 = L3;
#
# N1 = L1 + K1 - 1;
# N2 = L2 + K2 - 1;
# N3 = L3 + K3 - 1;
#
# d  = randn(Tv, N1, N2, N3);
# v  = randn(Tv, K1, K2, K3);
#
# H  = build_hankel_matrix(d);
# T  = build_toeplitz_matrix(d);
#
# r1 = H * vec(v);
# v_hat = reverse_order(v);
# r2 = T * vec(v_hat);
# norm(r1-r2) / norm(r1)
#
#
# # forward
# r1 = H * vec(v);
# r2 = hankel_multiplication(d, v; flag="forward");
# r2 = vec(r2);
# norm(r1-r2) / norm(r1)
#
# # adjoint
# v = randn(Tv, L1, L2, L3);
# r3 = H' * vec(v);
# r4 = hankel_multiplication(d, v; flag="adjoint");
# r4 = vec(r4);
# norm(r3-r4) / norm(r3)
#
# # data format
# Tv = Complex{Float64};
#
#
# #==============================================================================#
# #                                   4D case
# #==============================================================================#
#
# # data format
# Tv = Complex{Float64};
# # Tv = Float64;
#
# # dimensions
# L1 = 5; K1 = L1-1;
# L2 = 6; K2 = L2-1;
# L3 = 7; K3 = L3-1;
# L4 = 8; K4 = L4-1;
#
# # L1 = 5; K1 = L1;
# # L2 = 6; K2 = L2;
# # L3 = 7; K3 = L3;
# # L4 = 8; K4 = L4;
#
# N1 = L1 + K1 - 1;
# N2 = L2 + K2 - 1;
# N3 = L3 + K3 - 1;
# N4 = L4 + K4 - 1;
#
#
# d  = randn(Tv, N1, N2, N3, N4);
# v  = randn(Tv, K1, K2, K3, K4);
#
# H  = build_hankel_matrix(d);
# T  = build_toeplitz_matrix(d);
#
# r1 = H * vec(v);
# v_hat = reverse_order(v);
# r2 = T * vec(v_hat);
# norm(r1-r2) / norm(r1)
#
# # forward
# r1 = H * vec(v);
# r2 = hankel_multiplication(d, v; flag="forward");
# r2 = vec(r2);
# norm(r1-r2) / norm(r1)
#
# # adjoint
# v = randn(Tv, L1, L2, L3, L4);
# r3 = H' * vec(v);
# r4 = hankel_multiplication(d, v; flag="adjoint");
# r4 = vec(r4);
# norm(r3-r4) / norm(r3)
#
#
# #==============================================================================#
# #                                   5D case
# #==============================================================================#
# # data format
# # Tv = Complex{Float64};
# Tv = Float64;
#
# dimensions
# L1 = 4; K1 = L1-1;
# L2 = 4; K2 = L2-1;
# L3 = 4; K3 = L3-1;
# L4 = 4; K4 = L4-1;
# L5 = 4; K5 = L5-1;
#
# # L1 = 4; K1 = L1;
# # L2 = 4; K2 = L2;
# # L3 = 4; K3 = L3;
# # L4 = 4; K4 = L4;
# # L5 = 4; K5 = L5;
#
# N1 = L1 + K1 - 1;
# N2 = L2 + K2 - 1;
# N3 = L3 + K3 - 1;
# N4 = L4 + K4 - 1;
# N5 = L5 + K5 - 1;
#
# d  = randn(Tv, N1, N2, N3, N4, N5);
# v  = randn(Tv, K1, K2, K3, K4, K5);
#
# H  = build_hankel_matrix(d);
# T  = build_toeplitz_matrix(d);
#
# r1 = H * vec(v);
# v_hat = reverse_order(v);
# r2 = T * vec(v_hat);
# norm(r1-r2) / norm(r1)
