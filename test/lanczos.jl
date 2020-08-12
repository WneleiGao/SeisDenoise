# ======================================================================
#             test lanczos with partial bidiagonalization
# ======================================================================
using PyPlot, AcousticWave

m = 100; n = 67;
k = n;
A = randn(m, n);
(U, B, K) = lanbpro(A, k, r0=randn(m));

m = 212; n = 100;
path = "/Users/wenlei/Documents/AcousticHessianApprox/test/data/helio.bin"
fid  = open(path, "r"); A = read(fid, Float64, m*n); close(fid);
A = reshape(A, m, n);

k  = n

# the svd computed by lacpack
(ls, sv_true, rs) = svd(A);
A = ls * diagm(sv_true) * rs;
(ls, sv, rs) = svd(A)

# the relative error
relerr = zeros(k)
for i = 1 : k
    relerr[i] = abs(sv[i]-sv_true[i]) / sv_true[i]
end

# computed by matrix
r0 = rand(m) - 0.5*ones(m)
(U, B, V) = lanbpro(A, k, r0=r0)

# computed by operator
function Amul(x, iflag; A=[])

    if iflag == 1
       y = A * x
    elseif iflag == 2
       y = A' * x
    else
       error("no third option")
    end
    return y
end

params = Dict(:A=>A)
(U1, B1, V1) = lanbpro(Amul, k; m=m, n=n, params=params, r0=r0)

vecnorm(U1-U)
vecnorm(B1-B)
vecnorm(V1-V)


(ls1, sv1, rs1) = svd(full(B));
relerr1 = zeros(k)
for i = 1 : k
    relerr1[i] = abs(sv1[i]-sv_true[i]) / sv_true[i]
end
plot(1:k, log10.(abs.(sv_true)), label="true")
plot(1:k, log10.(abs.(sv     )), label="lacpack")
plot(1:k, log10.(abs.(sv1    )), label="lanbpro")

plot(1:k, log10.(relerr), label="lacpack relerr")
plot(1:k, log10.(relerr1), label="lanbpro relerr")
legend();

# check the profile of loss orthogonality
delta = sqrt(eps()/k)
eta   = (eps()^(3/4)) / sqrt(k)

# plot for vector V
j = 19
plot(1:j, log10(delta)*ones(j))
plot(1:j, log10(eta)*ones(j))
plot(1:j, log10.(abs.(optout[:iptV_his][  1:j,j])))
plot(1:j, log10.(abs.(optout[:iptV_true][ 1:j,j])))

# plot for vector U
figure();
plot(1:j, log10(delta)*ones(j))
plot(1:j, log10(eta)*ones(j))
plot(1:j, log10.(abs.(optout[:iptU_his][  1:j,j])))
plot(1:j, log10.(abs.(optout[:iptU_true][ 1:j,j])))


figure();
plot(1:k, log10(delta)*ones(k))
plot(1:k, log10(eta)*ones(k))
plot(1:k, log10.(abs.(optout[:iptUmax][1:k])))
plot(1:k, log10.(abs.(optout[:iptVmax][1:k])))


j = 6; figure();
plot(1:j, log10(delta)*ones(j))
plot(1:j, log10(eta)*ones(j))
plot(1:j, log10.(abs.(iptU_his[1:j,j])))
plot(1:j, log10.(abs.(iptU_true[1:j,j])))
plot(1:j, log10.(abs.(iptU_after[1:j,j])))


# ======================================================================
#               basic version of lanczos bidiagonalization
# ======================================================================
m = 375; n = 300;
A = randn(m, n)
k = 6
(U, B, V, R) = lanczos(A, k);

# A * V_K = U_k * B_k + U[:.k+1] * beta[k+1] * e_k'
residue = A * V - U * B
residue = residue[:,end]
norm(residue - R) / norm(residue)

plot(abs.(U'*R))

indexSet = collect(1:k)
(re, rnorm, num) = reorthogonalization(U, R, norm(R), indexSet, iflag=2)

plot(residue - re)
#
# # A' * U_k = V_k * B_k'
# residue = A' * U - V * B'
# vecnorm(residue)
#
# # inner product
# iptU = U' * U
# iptV = V' * V
# imshow(log10.(abs.(iptV-diagm(ones(k)))))

# ======================================================================
#                      test Reorthogonalization
# ======================================================================
m = 1777; n = 23;
Q = rand(m, n);
for i = 1 : n
    Q[:,i] = Q[:,i] / norm(Q[:,i])
end
(Q, s, v) = svd(Q)

r = rand(m)
IndexSet = collect(1:n)

(re , rnorm , num) = reorthogonalization(Q, r, norm(r), IndexSet, iflag=1)
(re1, rnorm1, num) = reorthogonalization(Q, r, norm(r), IndexSet, iflag=2)

plot(log10.(abs.(Q[:,IndexSet]'*re )))
plot(log10.(abs.(Q[:,IndexSet]'*re1)))


# test select vectors to orthogonalize
N = 20
ipt = rand(N); x = collect(1:N)
j   = 15
delta = 0.8
eta = 0.3
LL  = 0
strategy = 0
extra = 0
index = selectVectors(ipt, j, delta, eta, LL, 0, 0)

plot(x, ipt, "o"); plot(delta*ones(N+1)); plot(eta*ones(N+1), "-");
plot(j*ones(N), collect(linspace(0, 1, N)))
plot(LL*ones(N), collect(linspace(0, 1, N)))
scatter(index, ipt[index], marker="x", s=60, c="r")
