# ==============================================================================
#                        real data example
# ==============================================================================
using Registration

i3 = 24; i4 = 5; dt =0.004;
pin = "/Users/wenlei/Desktop/ArcisData/SmallCube/PP_wave/patches/patch"
pathin  = join([pin "_" "$i3" "_" "$i4" ".bin"]); fid = open(pathin, "r")
nt  = read(fid, Int32); ntrace = read(fid, Int32);
mark= read(fid, Int32, ntrace); mark = convert(Vector{Int64}, mark);
dim = read(fid, Int32, 4); dim = convert(Vector{Int64}, dim);
dobs= read(fid, Float32, nt*ntrace); dobs=reshape(dobs,Int64(nt),Int64(ntrace));
close(fid)
dori = zeros(Float32, nt, dim[1], dim[2], dim[3], dim[4]);
for i = 1 : length(mark)
    idx = ind2sub((dim[1],dim[2],dim[3],dim[4]), mark[i])
    dori[:,idx[1], idx[2], idx[3], idx[4]] = dobs[:,i]
end
dori = reshape(dori, Int64(nt), dim[1], dim[2], dim[3], dim[4]);
dori = permutedims(dori, (1,3,2,4,5));

itl=376; nt=165; itu=itl+nt-1;
dobs = convert(Array{Float64,5}, dori)[itl:itu, :,:,:,:];
(nt, N1, N2, N3, N4) = size(dobs);
(Amp, f) = freqSpectrum(reshape(dobs, nt, N1*N2*N3*N4), dt)

Fd = fft(dobs, 1)
f0 = 32.0; iw=floor(Int64, f0*nt*dt)+1;

Indicator = ones(N1, N2, N3, N4)
for i4 = 1 : N4
    for i3 = 1 : N3
        for i2 = 1 : N2
            for i1 = 1 : N1
                if norm(dobs[:,i1,i2,i3,i4]) == 0.0
                   Indicator[i1,i2,i3,i4] = 0.0
                end
            end
        end
    end
end

# find the row and column index of missing trace embeded in hankel matrix
L1 = floor(Int64, N1/2) + 1; K1 = N1 - L1 + 1;
L2 = floor(Int64, N2/2) + 1; K2 = N2 - L2 + 1;
L3 = floor(Int64, N3/2) + 1; K3 = N3 - L3 + 1;
L4 = floor(Int64, N4/2) + 1; K4 = N4 - L4 + 1;
# save the index the trace is observed
rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[];
# fourth layer
for j4 = 1 : K4
    col4 = (j4-1)*K3*K2*K1
    for i4 = 1 : L4
        row4 = (i4-1)*L3*L2*L1
        m4 = i4 + j4 - 1
        # third layer
        for j3 = 1 : K3
            col3 = (j3-1)*K2*K1
            for i3 = 1 : L3
                row3 = (i3-1)*L2*L1
                m3 = i3 + j3 - 1
                # second layer
                for j2 = 1 : K2
                    col2 = (j2-1)*K1
                    for i2 = 1 : L2
                        row2 = (i2-1)*L1
                        m2 = i2 + j2 - 1
                        # first layer
                        for j1 = 1 : K1
                            col1 = col4 + col3 + col2 + j1
                            for i1 = 1 : L1
                                row1 = row4 + row3 + row2 + i1
                                m1 = i1 + j1 - 1
                                if Indicator[m1, m2, m3, m4] != 0.0
                                   push!(rowIdx, row1)
                                   push!(colIdx, col1)
                                   push!(dc, Fd[iw, m1, m2, m3, m4])
                                end
                            end
                        end #first
                    end
                end #second
            end
        end #third
    end
end #fourth

mu = 1.0 / 100.0 * vecnorm(dc); nrow = L1*L2*L3*L4; ncol = K1*K2*K3*K4
rk = 50; maxIter = 50;
(count, sigma) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);

plot(collect(1:rk-1), abs.(diff(sort(abs.(sigma[1:rk]), rev=true))), label="mu=0.01")
plot(collect(1:rk), sort(abs.(sigma[1:rk]), rev=true), label="mu=0.01")


# read the reconstructed one
pout = "/Users/wenlei/Desktop/ArcisData/SmallCube/PP_wave/recon/recon"
pathout  = join([pout "_" "$i3" "_" "$i4" ".bin"]); fid = open(pathout, "r")
nt = read(fid, Int32); n1 = read(fid,Int32); n2 = read(fid,Int32); n3 = read(fid,Int32); n4 = read(fid,Int32)
drec = read(fid, Float32, nt*n1*n2*n3*n4); drec = reshape(drec, Int64(nt), Int64(n1), Int64(n2), Int64(n3), Int64(n4)); close(fid);
drec = permutedims(drec, (1,3,2,4,5));

# ploting assistant variables
sp = 2; noff = size(drec,2); naz = size(drec,3);
ntrace = noff*naz + sp*(naz-1)

itl=376; nt=165; itu=itl+nt-1
i4=8; i5=8;
tori = zeros(Float32,nt,ntrace)
trec = zeros(Float32,nt,ntrace)
for i = 1 : naz
    il = (i-1)*(sp+noff)+1; iu = il+noff-1;
    tori[:,il:iu] = dori[itl:itu,:,i,i4,i5]
    trec[:,il:iu] = drec[itl:itu,:,i,i4,i5]
end

tmp = "/Users/wenlei/Desktop/ori.bin" ; tmpfid = open(tmp, "w"); write(tmpfid, vec(tori)); close(tmpfid);
tmp = "/Users/wenlei/Desktop/rec.bin" ; tmpfid = open(tmp, "w"); write(tmpfid, vec(trec)); close(tmpfid);
pswigb < /Users/wenlei/Desktop/ori.bin xbox=0.5 ybox=0.5 labelsize=14 n1=165 d1=0.004 f1=1.5 d2=0.05263 d2num=1 f2=1 f2num=1 hbox=6 wbox=8 label1="Time (s)" label2="Azimuth Index" perc=96 interp=1 > /Users/wenlei/Desktop/ori1.eps
pswigb < /Users/wenlei/Desktop/rec.bin xbox=0.5 ybox=0.5 labelsize=14 n1=165 d1=0.004 f1=1.5 d2=0.05263 d2num=1 f2=1 f2num=1 hbox=6 wbox=8 label1="Time (s)" label2="Azimuth Index" perc=96 interp=1 > /Users/wenlei/Desktop/rec1.eps


# ==============================================================================
#                          test for 2D MSSA
# ==============================================================================
using PyPlot, Seismic
(d0, ext)=SeisLinearEvents(nt= 256, dt=0.008, nx1=32, dx1=30,
                         tau=[0.47, 0.83, 1.23],
                         p1=[0.00045, -0.00039, 0.0],
                         p2=[0.0, 0.0,  0.],
                         p3=[0., 0., 0.],
                         p4=[0., 0., 0.],
                         amp=[1., 0.7, -0.5])

SR = 0.5; Indicator = ones(n1);
for i = 1 : n1
    if rand() > SR
       d[:,i] = 0.0
       Indicator[i] = 0.0
    end
end

Fd = fft(d, 1)
f0=20.; dt=0.008; iw = floor(Int64, f0*dt*nt) + 1;
L = floor(Int64, n1/2) + 1; K = n1+1-L;

rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[];
for i2 = 1 : K
    for i1 = 1 : L
        idx = i1+i2-1
        if Indicator[idx] == 1.0
           println("$idx")
           push!(rowIdx, i1)
           push!(colIdx, i2)
           push!(dc, Fd[iw, idx])
        end
    end
end

H = zeros(L, K)
for i2 = 1 : K
    for i1 = 1 : L
        idx = i1+i2-1
        if Indicator[idx] == 1.0
           H[i1,i2] = 1.
        end
    end
end




rk=K; maxIter = 100; nrow = L; ncol = K;
mu = 1.0;  (count, w1 ) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 10.0; (count, w10) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 20.0; (count, w20) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 30.0; (count, w30) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 40.0; (count, w40) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 50.0; (count, w50) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 60.0; (count, w60) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 70.0; (count, w70) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 80.0; (count, w80) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 90.0; (count, w90) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu =100.0; (count, w100)= L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);

# plot(collect(1:10), abs.(diff((abs.(w[1:11])))));
plot(collect(1:rk), sort(abs.(w1[1:rk]), rev=true), label="mu=1")
plot(collect(1:rk), sort(abs.(w10[1:rk]), rev=true), label="mu=10")
plot(collect(1:rk), sort(abs.(w20[1:rk]), rev=true), label="mu=20")
plot(collect(1:rk), sort(abs.(w30[1:rk]), rev=true), label="mu=30")
plot(collect(1:rk), sort(abs.(w40[1:rk]), rev=true), label="mu=40")
plot(collect(1:rk), sort(abs.(w50[1:rk]), rev=true), label="mu=50")
plot(collect(1:rk), sort(abs.(w60[1:rk]), rev=true), label="mu=60")
plot(collect(1:rk), sort(abs.(w70[1:rk]), rev=true), label="mu=70")
plot(collect(1:rk), sort(abs.(w80[1:rk]), rev=true), label="mu=80")
plot(collect(1:rk), sort(abs.(w90[1:rk]), rev=true), label="mu=90")
plot(collect(1:rk), sort(abs.(w100[1:rk]), rev=true), label="mu=100")

plot(collect(1:rk-1), abs.(diff(sort(abs.(w1[1:rk]), rev=true))), label="mu=1")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w10[1:rk]), rev=true))), label="mu=10")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w20[1:rk]), rev=true))), label="mu=20")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w30[1:rk]), rev=true))), label="mu=30")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w40[1:rk]), rev=true))), label="mu=40")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w50[1:rk]), rev=true))), label="mu=50")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w60[1:rk]), rev=true))), label="mu=60")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w70[1:rk]), rev=true))), label="mu=70")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w80[1:rk]), rev=true))), label="mu=80")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w90[1:rk]), rev=true))), label="mu=90")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w100[1:rk]), rev=true))), label="mu=100")


# ==============================================================================
#                          test super parameter for 3D MSSA
# ==============================================================================
using PyPlot, Seismic
(d0, ext)=SeisLinearEvents(nt= 250, dt=0.008, nx1=32, nx2=32, dx1=30, dx2=30,
                         tau=[0.4, 0.8, 1.2, 1.6],
                         p1=[ 0.0003,  0.00057, -0.0004, 0.],
                         p2=[-0.0002, -0.0004,  0.0001, 0.],
                         p3=[0., 0., -0., -0.],
                         p4=[0., 0., -0., -0.],
                         amp=[1., 1.2, -0.9, -0.9])

snr = 1.0; d = SeisAddNoise(d0, snr, L=5);
(nt, N1, N2) = size(d)

Fd = fft(d, 1);
f0 = 20.0; dt = 0.008; iw = floor(Int64, f0*dt*nt) + 1;

SR = 0.3; Indicator = ones(N1, N2)
for i2 = 1 : N2
    for i1 = 1 : N1
        if rand() > SR
           Indicator[i1,i2] = 0.
        end
    end
end

L1 = floor(Int64, N1/2) + 1; K1 = N1 - L1 + 1;
L2 = floor(Int64, N2/2) + 1; K2 = N2 - L2 + 1;

rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[]
for j2 = 1 : K2
    col2 = (j2-1)*K1
    for i2 = 1 : L2
        row2 = (i2-1)*L1
        m2 = i2 + j2 - 1
        # first layer
        for j1 = 1 : K1
            col1 = col2 + j1
            for i1 = 1 : L1
                row1 = row2 + i1
                m1 = i1 + j1 - 1
                if Indicator[m1, m2] == 1.0
                   push!(rowIdx, row1)
                   push!(colIdx, col1)
                   push!(dc, Fd[iw, m1, m2])
                end
            end
        end #first
    end
end #second

rk=30; maxIter = 50; nrow = L1 * L2; ncol = K1 * K2;
denum = collect(1000:-100:100); append!(denum, 50); append!(denum, 10);

W = Vector{Vector}(length(denum))
for i = 1 : length(denum)
    mu = norm(dc) / denum[i];
    (count, W[i]) = L1MC_improvedInitial(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
end

# plot(collect(1:10), abs.(diff((abs.(w[1:11])))));

for i = 1 : 12
    mu = denum[i]
    plot(collect(1:rk), sort(abs.(W[i][1:rk]), rev=true), label="mu=$mu"); legend();
    path = join(["/Users/wenlei/Desktop/fig" "$i" ".pdf"]);
    savefig(path); close();
end

for i = 1 : 12
    mu = 1/denum[i]
    plot(collect(1:rk-1), abs.(diff(sort(abs.(W[i][1:rk]), rev=true))), label="SR=$mu"); legend();
    path = join(["/Users/wenlei/Desktop/fig" "$i" ".pdf"]);
    savefig(path); close();
end

tmp="/Users/wenlei/Desktop/d0.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d0))); close(tmpfid);
tmp="/Users/wenlei/Desktop/dn.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d ))); close(tmpfid);
tmp="/Users/wenlei/Desktop/dobs.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d ))); close(tmpfid);

pscube < /Users/wenlei/Desktop/d0.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/d0.eps
pscube < /Users/wenlei/Desktop/dn.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/dn.eps
pscube < /Users/wenlei/Desktop/dobs.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=32 d2num=8 f2num=8 label2="X" n3=32 d3num=8 f3num=8 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/dobs.eps


# ==============================================================================
#               test for 3D MSSA with different sampling rate
# ==============================================================================
using PyPlot, Seismic
(d0, ext)=SeisLinearEvents(nt= 250, dt=0.008, nx1=32, nx2=32, dx1=30, dx2=30,
                           tau=[0.4, 0.8, 1.2, 1.6],
                           p1=[ 0.0003,  0.00057, -0.0004, 0.],
                           p2=[-0.0002, -0.0004,  0.0001, 0.],
                           p3=[0., 0., -0., -0.],
                           p4=[0., 0., -0., -0.],
                           amp=[1., -0.9, 1.2, -0.9])

snr = 1.0; d = SeisAddNoise(d0, snr, L=5);
(nt, N1, N2) = size(d);

Fd = fft(d, 1); f0=20.0; dt= 0.008;
iw = floor(Int64, f0*dt*nt) + 1;

L1 = floor(Int64, N1/2) + 1; K1 = N1 - L1 + 1;
L2 = floor(Int64, N2/2) + 1; K2 = N2 - L2 + 1;
nrow = L1 * L2; ncol = K1 * K2;

rk=30; maxIter = 50;
# test different sampling rate
SampleRate = collect(0.05:0.05:0.9)
W = Vector{Vector}(length(SampleRate))
for s = 1 : length(SampleRate)
    SR = SampleRate[s]
    Indicator = ones(N1, N2)
    for i2 = 1 : N2
        for i1 = 1 : N1
            if rand() > SR
               Indicator[i1,i2] = 0.0;
            end
        end
    end
    rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[];
    # fourth layer
    for j2 = 1 : K2
        col2 = (j2-1)*K1
        for i2 = 1 : L2
            row2 = (i2-1)*L1
            m2 = i2 + j2 - 1
            # first layer
            for j1 = 1 : K1
                col1 = col2 + j1
                for i1 = 1 : L1
                    row1 = row2 + i1
                    m1 = i1 + j1 - 1
                    if Indicator[m1, m2] == 1.0
                       push!(rowIdx, row1)
                       push!(colIdx, col1)
                       push!(dc, Fd[iw, m1, m2])
                    end
                end
            end #first
        end
    end #second
    mu = 1.0 / 50.0 * vecnorm(dc)
    (count, W[s]) = L1MC_improvedInitial(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
    # (count, W[s]) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
    println("finish $s")
end

for i = 1 : 12
    sr = i * 5
    plot(collect(1:rk), sort(abs.(W[i][1:rk]), rev=true), label="SR=$sr"); legend();
    path = join(["/Users/wenlei/Desktop/fig" "$i" ".pdf"]);
    savefig(path); close();
end

for i = 1 : 12
    sr = i * 5
    plot(collect(1:rk-1), abs.(diff(sort(abs.(W[i][1:rk]), rev=true))), label="SR=$sr"); legend();
    path = join(["/Users/wenlei/Desktop/fig" "$i" ".pdf"]);
    savefig(path); close();
end

# ==============================================================================
#              5D examples
# ==============================================================================
# generate 5D synthetic seismic data
(d0, ext) = SeisLinearEvents(dt=0.008, nt=256,
                     dx1=30.0, nx1=12,
                     dx2=30.0, nx2=12,
                     dx3=30.0, nx3=12,
                     dx4=30.0, nx4=12,
                     tau=[0.32, 0.85, 1.13, 1.53, 1.8],
                     p1 =[ 0.0003,  0.0004, -0.0004, 0.,  0.0001],
                     p2 =[-0.0002, -0.0003,  0.0002, 0., -0.0002],
                     p3 =[ 0.0002,  0.0003, -0.0002, 0.,  0.0002],
                     p4 =[-0.0004, -0.0002,  0.0003, 0., -0.0003],
                     amp = [1, -1, 0.7, 0.5, 0.89], f0=20.0)
(nt, N1, N2, N3, N4) = size(d0)
snr = 1.0; d = SeisAddNoise(d0, snr, L=5);

# randomly sampling the data
SR = 0.3; dobs = copy(d); Indicator = ones(N1, N2, N3, N4)
for i4 = 1 : N4
    for i3 = 1 : N3
        for i2 = 1 : N2
            for i1 = 1 : N1
                if rand() > SR
                   dobs[:,i1,i2,i3,i4] = 0.0;
                   Indicator[i1,i2,i3,i4] = 0;
                end
            end
        end
    end
end

# transform data to frequency domain
Fd = fft(dobs, 1); f0 = 20.0; dt = 0.008;
iw = floor(Int64, f0*nt*dt) + 1;

# find the row and column index of missing trace embeded in hankel matrix
L1 = floor(Int64, N1/2) + 1; K1 = N1 - L1 + 1;
L2 = floor(Int64, N2/2) + 1; K2 = N2 - L2 + 1;
L3 = floor(Int64, N3/2) + 1; K3 = N3 - L3 + 1;
L4 = floor(Int64, N4/2) + 1; K4 = N4 - L4 + 1;
# save the index the trace is observed
rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[];
# fourth layer
for j4 = 1 : K4
    col4 = (j4-1)*K3*K2*K1
    for i4 = 1 : L4
        row4 = (i4-1)*L3*L2*L1
        m4 = i4 + j4 - 1
        # third layer
        for j3 = 1 : K3
            col3 = (j3-1)*K2*K1
            for i3 = 1 : L3
                row3 = (i3-1)*L2*L1
                m3 = i3 + j3 - 1
                # second layer
                for j2 = 1 : K2
                    col2 = (j2-1)*K1
                    for i2 = 1 : L2
                        row2 = (i2-1)*L1
                        m2 = i2 + j2 - 1
                        # first layer
                        for j1 = 1 : K1
                            col1 = col4 + col3 + col2 + j1
                            for i1 = 1 : L1
                                row1 = row4 + row3 + row2 + i1
                                m1 = i1 + j1 - 1
                                if Indicator[m1, m2, m3, m4] != 0.0
                                   push!(rowIdx, row1)
                                   push!(colIdx, col1)
                                   push!(dc, Fd[iw, m1, m2, m3, m4])
                                end
                            end
                        end #first
                    end
                end #second
            end
        end #third
    end
end #fourth
# length(rowIdx) / (L1*L2*L3*L4*K1*K2*K3*K4)

rk=30; maxIter = 25; nrow = L1*L2*L3*L4; ncol = K1*K2*K3*K4;
mu = 1.0;  (count, w1 ) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 10.0; (count, w10) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 20.0; (count, w20) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 30.0; (count, w30) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 40.0; (count, w40) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 50.0; (count, w50) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 60.0; (count, w60) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 70.0; (count, w70) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 80.0; (count, w80) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu = 90.0; (count, w90) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
mu =100.0; (count, w100)= L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);

# plot(collect(1:10), abs.(diff((abs.(w[1:11])))));
plot(collect(1:rk), sort(abs.(w1[1:rk]), rev=true), label="mu=1")
plot(collect(1:rk), sort(abs.(w10[1:rk]), rev=true), label="mu=10")
plot(collect(1:rk), sort(abs.(w20[1:rk]), rev=true), label="mu=20")
plot(collect(1:rk), sort(abs.(w30[1:rk]), rev=true), label="mu=30")
plot(collect(1:rk), sort(abs.(w40[1:rk]), rev=true), label="mu=40")
plot(collect(1:rk), sort(abs.(w50[1:rk]), rev=true), label="mu=50")
plot(collect(1:rk), sort(abs.(w60[1:rk]), rev=true), label="mu=60")
plot(collect(1:rk), sort(abs.(w70[1:rk]), rev=true), label="mu=70")
plot(collect(1:rk), sort(abs.(w80[1:rk]), rev=true), label="mu=80")
plot(collect(1:rk), sort(abs.(w90[1:rk]), rev=true), label="mu=90")
plot(collect(1:rk), sort(abs.(w100[1:rk]), rev=true), label="mu=100"); legend();

plot(collect(1:rk-1), abs.(diff(sort(abs.(w1[1:rk]), rev=true))), label="mu=1")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w10[1:rk]), rev=true))), label="mu=10")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w20[1:rk]), rev=true))), label="mu=20")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w30[1:rk]), rev=true))), label="mu=30")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w40[1:rk]), rev=true))), label="mu=40")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w50[1:rk]), rev=true))), label="mu=50")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w60[1:rk]), rev=true))), label="mu=60")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w70[1:rk]), rev=true))), label="mu=70")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w80[1:rk]), rev=true))), label="mu=80")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w90[1:rk]), rev=true))), label="mu=90")
plot(collect(1:rk-1), abs.(diff(sort(abs.(w100[1:rk]), rev=true))), label="mu=100"); legend();

tmp="/Users/wenlei/Desktop/d0.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d0[:,:,:,12,12]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/d.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(d[:,:,:,12,12]))); close(tmpfid);
tmp="/Users/wenlei/Desktop/dobs.bin"; tmpfid=open(tmp,"w"); write(tmpfid,convert(Vector{Float32}, vec(dobs[:,:,:,12,12]))); close(tmpfid);

pscube < /Users/wenlei/Desktop/d0.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=256 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=12 d2num=3 f2num=3 label2="X" n3=12 d3num=3 f3num=3 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/d0.eps
pscube < /Users/wenlei/Desktop/d.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=256 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=12 d2num=3 f2num=3 label2="X" n3=12 d3num=3 f3num=3 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/d.eps
pscube < /Users/wenlei/Desktop/dobs.bin size1=3.0 size2=2.5 size3=1.2 labelsize=12 n1=256 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=12 d2num=3 f2num=3 label2="X" n3=12 d3num=3 f3num=3 label3="Y" xbox=0.5 ybox=0.5 perc=0.96 > /Users/wenlei/Desktop/dobs.eps


# ==============================================================================
#              5D examples, test different sampling rate
# ==============================================================================
# generate 5D synthetic seismic data
(d0, ext) = SeisLinearEvents(dt=0.008, nt=256,
                     dx1=30.0, nx1=12,
                     dx2=30.0, nx2=12,
                     dx3=30.0, nx3=12,
                     dx4=30.0, nx4=12,
                     tau=[0.32, 0.85, 1.13, 1.53, 1.8],
                     p1 =[ 0.0003,  0.0004, -0.0004, 0.,  0.0001],
                     p2 =[-0.0002, -0.0003,  0.0002, 0., -0.0002],
                     p3 =[ 0.0002,  0.0003, -0.0002, 0.,  0.0002],
                     p4 =[-0.0004, -0.0002,  0.0003, 0., -0.0003],
                     amp = [1, -1, 0.7, 0.5, 0.89], f0=20.0)
(nt, N1, N2, N3, N4) = size(d0)
snr = 1.0; dn = SeisAddNoise(d0, snr, L=5);
Fd = fft(dn, 1); f0 = 20.0; dt = 0.008;
iw = floor(Int64, f0*nt*dt) + 1;


# test different sampling rate
SampleRate = collect(0.05:0.05:0.9)
W = Vector{Vector}(length(SampleRate))
for s = 1 : length(SampleRate)
    SR = SampleRate[s]
    Indicator = ones(N1, N2, N3, N4)
    for i4 = 1 : N4
        for i3 = 1 : N3
            for i2 = 1 : N2
                for i1 = 1 : N1
                    if rand() > SR
                       Indicator[i1,i2,i3,i4] = 0.0;
                    end
                end
            end
        end
    end
    rowIdx = Int64[]; colIdx = Int64[]; dc = Complex128[];
    # fourth layer
    for j4 = 1 : K4
        col4 = (j4-1)*K3*K2*K1
        for i4 = 1 : L4
            row4 = (i4-1)*L3*L2*L1
            m4 = i4 + j4 - 1
            # third layer
            for j3 = 1 : K3
                col3 = (j3-1)*K2*K1
                for i3 = 1 : L3
                    row3 = (i3-1)*L2*L1
                    m3 = i3 + j3 - 1
                    # second layer
                    for j2 = 1 : K2
                        col2 = (j2-1)*K1
                        for i2 = 1 : L2
                            row2 = (i2-1)*L1
                            m2 = i2 + j2 - 1
                            # first layer
                            for j1 = 1 : K1
                                col1 = col4 + col3 + col2 + j1
                                for i1 = 1 : L1
                                    row1 = row4 + row3 + row2 + i1
                                    m1 = i1 + j1 - 1
                                    if Indicator[m1, m2, m3, m4] == 1.0
                                       push!(rowIdx, row1)
                                       push!(colIdx, col1)
                                       push!(dc, Fd[iw, m1, m2, m3, m4])
                                    end
                                end
                            end #first
                        end
                    end #second
                end
            end #third
        end
    end #fourth
    mu = 1.0 / 50.0 * vecnorm(dc)
    (count, W[s]) = L1MC(dc, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
    println("finish $s")
end

for i = 1 : length(SampleRate)
    sr = i * 5
    plot(collect(1:rk), sort(abs.(W[i][1:rk]), rev=true), label="SR=$sr")
end
legend(); savefig("/Users/wenlei/Desktop/fig6.pdf");

for i = 1 : 18
    sr = i * 5
    plot(collect(1:rk-1), abs.(diff(sort(abs.(W[i][1:rk]), rev=true))), label="SR=$sr")
end
legend(); savefig("/Users/wenlei/Desktop/fig7.pdf");


# ==============================================================================
#                         brute force to determine turning point (works)
# ==============================================================================
w = sort(abs.(w10[1:30]), rev=true);
function DetectTurningPoint(w::Vector{Float64})
    n = length(w)
    cost = zeros(n-2)
    for k = 2 : n-2
          b  = w[1:k]
          A  = hcat(collect(1:k), ones(k))
          lc = (A'*A) \ (A'*b)

          # fit the right line
          b  = w[k+1:n]
          A  = hcat(collect(k+1:n), ones(n-k))
          rc = (A'*A) \ (A'*b)
          for i = 1 : k
              cost[k] = cost[k] + (w[i]-lc[1]*i - lc[2])^2
          end
          for i = k+1 : n
              cost[k] = cost[k] + (w[i]-rc[1]*i - rc[2])^2
          end
     end
end

# ==============================================================================
#                          test for real data
# ==============================================================================
M = 377; N = 412; R= 27;
A = randn(M, R); B = randn(N, R);
C = A * B'; C = SeisAddNoise(C, 1.5, L=5);
SR = 0.35;
rowIdx = Int64[];
colIdx = Int64[];
d      = Float64[];
for j = 1 : N
    for i = 1 : M
        if rand() < SR
           push!(rowIdx, i)
           push!(colIdx, j)
           push!(d, C[i,j])
        end
    end
end

rk = 100; maxIter = 100; mu = 50.0;
(count, w) = L1MC(d, rowIdx, colIdx, M, N, rk, mu, maxIter);
plot(abs.(diff(abs.(w))));
plot(abs.(w))

# ==============================================================================
#                          test for complex data
# ==============================================================================
M = 377; N = 412; R= 27;
A = randn(M,R)+im*randn(M,R);
B = randn(N,R)+im*randn(N,R);
C = A * B';
SR = 0.35;
rowIdx = Int64[];
colIdx = Int64[];
d      = Complex128[];
for j = 1 : N
    for i = 1 : M
        if rand() < SR
           push!(rowIdx, i)
           push!(colIdx, j)
           push!(d, C[i,j])
        end
    end
end
rk=100; mu = 50.0; maxIter = 100;

# ==============================================================================
#                          test for complex data
# ==============================================================================
using PyPlot, Seismic
# generate clean data with 3 linear events
f0 = 20.0; dt = 0.004; tmax= 1.2; t=collect(0.0:dt:tmax); nt = length(t);
tau = [0.1  ,  0.4   , 0.55    ];
p   = [0.001, -0.0004, 0.000001]; nx =60; dx = 10.0; ph=zeros(length(p));
amp = [1.0  , -1.0   , 1.0     ];
(dobs, extent) = SeisLinearEvents(dt=dt, nt=nt, dx1=10, nx1=nx, tau=tau, p1=p, p2=ph, p3=ph, p4=ph, f0=f0, amp=amp)
snr = 1.5; dobs = SeisAddNoise(dobs, snr, L=5);

# add band-liminited gaussian noise to data, snr is defined as the ratio of maximum amplitude
(nt, nx) = size(dobs)
nf = nextpow2(nt);
nf = nf >= 512 ? nf : 512;
if nf > nt
   df = vcat(convert(Matrix{Complex{Float64}},dobs), zeros(Complex{Float64},nf-nt,nx))
else
   df = convert(Matrix{Complex{Float64}},dobs)
end
fft!(df, 1)  #in-place fft
ncol = floor(Int64, nx/2) + 1
nrow = nx - ncol + 1

iw = 40; SR=0.6
(rowIdx, colIdx, dobs) = Hankel(nrow, ncol, df, iw, SR)
rk=30; mu = 30.0; maxIter = 100;
(count, w) = L1MC(dobs, rowIdx, colIdx, nrow, ncol, rk, mu, maxIter);
plot(abs.(diff((abs.(w)))));
plot(abs.(w))


# ==============================================================================
#                          Float64 cases
# ==============================================================================
function FormHankel2D(df::Array{Complex{Tv}, 3}, iw::Ti, SR::Tv,
                      nrow2::Ti, ncol2::Ti, nrow3::Ti, ncol3::Ti) where {Tv<:Float64, Ti<:Int64}
     nx1 = nrow2 + ncol2 - 1
     nx2 = nrow3 + ncol3 - 1
     idx1 = Int64[]
     idx2 = Int64[]
     rowIdx = Int64[]
     colIdx = Int64[]
     dobs   = Complex128[]
     for j = 1 : nx2
         for i = 1 : nx1
             if rand() < SR
                push!(idx1, i)
                push!(idx2, j)
             end
         end
     end
     for j = 1 : ncol3
         for i = 1 : nrow3
             tmp2 = i+j-1
             rnum = (i-1)*nrow2
             cnum = (j-1)*ncol2
             for k = 1 : ncol2
                 for l = 1 : nrow2
                     tmp1 = k + l - 1
                     for m = 1 : length(idx1)
                         if tmp1 == idx1[m] && tmp2 == idx2[m]
                            irh  = rnum + l
                            irc  = cnum + k
                            push!(rowIdx, irh)
                            push!(colIdx, irc)
                            push!(dobs, df[iw, tmp1, tmp2])
                         end
                     end
                 end
             end
         end
     end
     return rowIdx, colIdx, dobs
end


function Hankel(nrow::Ti, ncol::Ti, df::Matrix{Complex{Tv}}, iw::Ti, SR::Tv) where {Tv<:Float64, Ti<:Int64}
    nx = nrow + ncol - 1
    H = zeros(Complex{Float64}, nrow, ncol)
    idx = Int64[];
    for i = 1 : nx
        if rand() < SR
           push!(idx, i)
        end
    end
    rowIdx = Int64[]
    colIdx = Int64[]
    dt     = Complex128[]
    for i2 = 1 : ncol
        for i1 = 1 : nrow
            tmp= i1 + i2 - 1
            for k = 1 : length(idx)
                if tmp == idx[k]
                   push!(rowIdx, i1)
                   push!(colIdx, i2)
                   push!(dt, df[iw, tmp])
                end
            end
        end
    end
    return rowIdx, colIdx, dt
end

function L1MC(d::Vector{Tv}, rowIdx::Vector{Ti}, colIdx::Vector{Ti}, M::Ti, N::Ti,
              rk::Ti, mu::Tv, maxIter::Ti) where {Tv<:Float64, Ti<:Int64}
    # initialization
    w = randn(rk)
    U = randn(M, rk)
    V = randn(N, rk)
    for i = 1 : rk
        U[:,i] = U[:,i] / norm(U[:,i])
        V[:,i] = V[:,i] / norm(V[:,i])
    end
    X = zeros(M, N)
    Xr= zeros(M, N)
    Z = zeros(M, N)
    for i = 1 : length(d)
        X[rowIdx[i], colIdx[i]] = d[i]
    end
    # iterations
    for iter = 1 : maxIter
        Xr = copy(X)
        for i = 1 : rk
            if w[i] != 0
               # update column factors
               U[:,i] = (Xr * V[:,i]) / w[i]
               U[:,i] = U[:,i] / norm(U[:,i])
               # update row    factors
               V[:,i] = (Xr'* U[:,i]) / w[i]
               V[:,i] = V[:,i] / norm(V[:,i])
               # update weight by soft shrinking
               w[i]   = sum(Xr .* (U[:,i] * V[:,i]'))
               if w[i] > mu
                  w[i] = w[i] - mu
               elseif w[i] < -mu
                  w[i] = w[i] + mu
               else
                  w[i] = 0.0
               end
               Xr = Xr - (w[i]*U[:,i]*V[:,i]')
            end
        end
        # update X
        Z = X - Xr; res = 0.0
        X = copy(Z)
        for i = 1 : length(d)
            X[rowIdx[i], colIdx[i]] = d[i]
            res = res + (Z[rowIdx[i], colIdx[i]]-d[i])^2
        end
        res = sqrt(res) / norm(d)
        println("$iter relerr: $res")
    end
    # final trim of rank
    threshold = 1e-3 * length(d) / (M*N) * sum(abs.(w))
    count = 0;
    for i = 1 : rk
        if abs(w[i]) > threshold
           count = count + 1
        else
           w[i] = 0.0
        end
    end
    return count, w
end

# ==============================================================================
#                          Complex128 cases
# ==============================================================================
function L1MC(d::Vector{Tv}, rowIdx::Vector{Ti}, colIdx::Vector{Ti}, M::Ti, N::Ti,
              rk::Ti, mu::Float64, maxIter::Ti) where {Tv<:Complex128, Ti<:Int64}
    # initialization
    w = randn(  rk) + im*randn(  rk)
    U = randn(M,rk) + im*randn(M,rk)
    V = randn(N,rk) + im*randn(N,rk)
    # normalize factor vectors
    for i = 1 : rk
        U[:,i] = U[:,i] / norm(U[:,i])
        V[:,i] = V[:,i] / norm(V[:,i])
    end
    X = zeros(Complex128, M, N)
    Xr= zeros(Complex128, M, N)
    Z = zeros(Complex128, M, N)
    # insert observations
    for i = 1 : length(d)
        X[rowIdx[i], colIdx[i]] = d[i]
    end
    # iterations
    for iter = 1 : maxIter
        Xr = copy(X)
        for i = 1 : rk
            if norm(w[i]) > 0.0
               # update column factors
               U[:,i] = conj(w[i]) / (norm(w[i]))^2 * (Xr * V[:,i])
               U[:,i] = U[:,i] / norm(U[:,i])
               # update row    factors
               V[:,i] =      w[i]  / (norm(w[i]))^2 * (Xr'* U[:,i])
               V[:,i] = V[:,i] / norm(V[:,i])
               # update weight by soft shrinking
               w[i]   = sum(Xr .* conj(U[:,i] * V[:,i]'))
               if norm(w[i]) > mu
                  w[i] = (norm(w[i])-mu)/norm(w[i]) * w[i]
               else
                  w[i] = zero(Complex128)
               end
               Xr = Xr - (w[i]*U[:,i]*V[:,i]')
            end
        end
        # update X
        Z = X - Xr; res = 0.0
        X = copy(Z)
        for i = 1 : length(d)
            X[rowIdx[i], colIdx[i]] = d[i]
            res = res + (norm(Z[rowIdx[i], colIdx[i]]-d[i]))^2
        end
        res = sqrt(res) / norm(d)
        println("$iter, relative error: $res")
    end
    # final trim of rank
    threshold = 1e-3 * length(d) / (M*N) * sum(abs.(w))
    count = 0;
    for i = 1 : rk
        if abs(w[i]) > threshold
           count = count + 1
        else
           w[i] = zero(Complex128)
        end
    end
    return count, w
end


function L1MC_improvedInitial(d::Vector{Tv}, rowIdx::Vector{Ti}, colIdx::Vector{Ti}, M::Ti, N::Ti,
              rk::Ti, mu::Float64, maxIter::Ti) where {Tv<:Complex128, Ti<:Int64}

    X = zeros(Complex128, M, N)
    Xr= zeros(Complex128, M, N)
    Z = zeros(Complex128, M, N)
    # insert observations
    for i = 1 : length(d)
        X[rowIdx[i], colIdx[i]] = d[i]
    end
    (U, w, V) = svd(X)
    U = U[:,1:rk]
    w = w[1:rk]; w = convert(Vector{Complex128}, w)
    V = V[:,1:rk]
    # iterations
    for iter = 1 : maxIter
        Xr = copy(X)
        for i = 1 : rk
            if norm(w[i]) > 0.0
               # update column factors
               U[:,i] = conj(w[i]) / (norm(w[i]))^2 * (Xr * V[:,i])
               U[:,i] = U[:,i] / norm(U[:,i])
               # update row    factors
               V[:,i] =      w[i]  / (norm(w[i]))^2 * (Xr'* U[:,i])
               V[:,i] = V[:,i] / norm(V[:,i])
               # update weight by soft shrinking
               w[i]   = sum(Xr .* conj(U[:,i] * V[:,i]'))
               if norm(w[i]) > mu
                  w[i] = (norm(w[i])-mu)/norm(w[i]) * w[i]
               else
                  w[i] = zero(Complex128)
               end
               Xr = Xr - (w[i]*U[:,i]*V[:,i]')
            end
        end
        # update X
        Z = X - Xr; res = 0.0
        X = copy(Z)
        for i = 1 : length(d)
            X[rowIdx[i], colIdx[i]] = d[i]
            res = res + (norm(Z[rowIdx[i], colIdx[i]]-d[i]))^2
        end
        res = sqrt(res) / norm(d)
        println("$iter, relative error: $res")
    end
    # final trim of rank
    threshold = 1e-3 * length(d) / (M*N) * sum(abs.(w))
    count = 0;
    for i = 1 : rk
        if abs(w[i]) > threshold
           count = count + 1
        else
           w[i] = zero(Complex128)
        end
    end
    return count, w
end


# ==============================================================================
#                     learn the properties of circulant matrix
# ==============================================================================
n = 128; c = rand(Complex128, n);
# c = convert(Vector{Complex{Float64}}, collect(1:n))
C = zeros(Complex128, n, n)
for k = 1 : n
    for i = k : n
        C[i,k] = c[i+1-k]
    end
    for i = 1 : k-1
        C[i,k] = c[n-k+i+1]
    end
end

# the eigenvector of circulant matrix is the basis vector of Fourier transform
# the corresponding eigenvalue is the fourier coefficients of the vector

# build fourier transform matrix
Fn = zeros(Complex128, n, n)
for i2 = 1 : n
    for i1 = 1 : n
        Fn[i1,i2] = exp(im*2*pi*(i2-1)*(i1-1)/n) / n
    end
end

coef = fft(c)
T  = Fn * diagm(coef)

T1 = C * Fn



# the multiplication of circulant matrix with a vector is
# C * v = ifft * (Fn * c) \cdot (Fn * v)
v = rand(Complex128, n)
b = C * v
b1= ifft( fft(c) .*  fft(v))
(vecnorm(b1-b)) / (vecnorm(b))

# ==============================================================================
#                    test the properties of Toeplitze matrix
# ==============================================================================
L = 556; K = 392;
N = L + K - 1;

t = rand(N); T = zeros(eltype(t), L, K);
# form toeplitze matrix
for j = 1 : K
    istart = K - j + 1
    for i = 1 : L
        T[i,j] = t[istart+i-1]
    end
end

v = rand(eltype(t), K);
b = T * v;
c = zeros(eltype(t), N); c[1:L] = t[K:end]; c[L+1:end] = t[1:K-1];
v1 = vcat(v, zeros(eltype(t), L-1));
b1= ifft(fft(c) .* fft(v1))[1:L];
(vecnorm(b1-b)) / (vecnorm(b))

# ==============================================================================
#                    test the properties of Hankel matrix (1D) case
# ==============================================================================
N = 1234;
K = floor(Int64, N/2) + 1; L = N - K + 1;
d = rand(Complex128, N)
H = zeros(eltype(d), L, K)
for j = 1 : K
    for i = 1 : L
        H[i,j] = d[i+j-1]
    end
end

 # compute H \times v
v = rand(eltype(d), K);
b = H * v
c = vcat(d[K:end], d[1:K-1]);
vhat = vcat(reverse(v), zeros(eltype(d), L-1))
b1 = ifft(fft(c) .* fft(vhat))[1:L]
(vecnorm(b1-b)) / (vecnorm(b))

# test compute the conjugate transpose of H times a vector
v = rand(L);
r = H' * v
c = conj(vcat(d[K:end], d[1:K-1])); # very important to use the conjugate of the vector
w = vcat(zeros(eltype(d),K-1), reverse(v))
r1 = ifft(fft(c) .* fft(w))[L:N]
(vecnorm(r1-r)) / (vecnorm(r))

# an equvialent algorithm to compute the conjugate transpose of H times a vector
c = conj(vcat(d[L:N], d[1:L-1]))
w = vcat(reverse(v), zeros(eltype(d), K-1))
r2 = ifft(fft(c) .* fft(w))[1:K]
(vecnorm(r2-r)) / (vecnorm(r))


# ==============================================================================
# test the multiplication of Hankel matrix with a vector and its transpose (2D) case
# ==============================================================================
(data, ext) = SeisParabEvents(nx1=51, nx2=51,p1=[0.3,-0.8,0.7], p2=[-0.4, 0.6, -0.7])

space = 30.0; N=20;
x1 = collect(0.0:space:space*N)
tmp = randn(length(x1))*10;
x1 = x1 + tmp

x2 = collect(0.0:space:space*N)
tmp = randn(length(x1))*5;
x2 = x2 + tmp
ds = SeisParab3D(x1, x2, apx1=0.*ones(3), apx2=0.*ones(3))

# d = SeisHyp3D();
# idx1 = sort(shuffle(collect(1:51))[1:17])
# idx2 = sort(shuffle(collect(1:51))[1:17])
# ds = data[:,idx1, idx2]

df = fft(ds,1);
idx = 40;

d = df[idx,:,:]
(N1, N2) = size(d)

K1 = floor(Int64, N1/2) + 1; L1 = N1 - K1 + 1;
K2 = floor(Int64, N2/2) + 1; L2 = N2 - K2 + 1;
H = zeros(eltype(d), L1*L2, K1*K2)

for n = 1 : K2
    for m = 1 : L2
        row = (m-1)*L1
        col = (n-1)*K1
        i2  = m + n -1
        for j = 1 : K1
            for i = 1 : L1
                i1 = i + j - 1
                H[row+i, col+j] = d[i1, i2]
            end
        end
    end
end
(U, S, V) = svd(H)


N1 = 33; N2 = 66; d = rand(Complex128, N1, N2);
K1 = floor(Int64, N1/2) + 1; L1 = N1 - K1 + 1;
K2 = floor(Int64, N2/2) + 1; L2 = N2 - K2 + 1;
H = zeros(eltype(d), L1*L2, K1*K2)

for n = 1 : K2
    for m = 1 : L2
        row = (m-1)*L1
        col = (n-1)*K1
        i2  = m + n -1
        for j = 1 : K1
            for i = 1 : L1
                i1 = i + j - 1
                H[row+i, col+j] = d[i1, i2]
            end
        end
    end
end

# direct multiplication
x = rand(eltype(d), K1*K2)
r = H * x

# fast multiplication
C = zeros(eltype(d), N1, N2)
for i2 = 1 : N2
    j2 = i2 <= L2 ? N2-L2+i2 : i2-L2
    for i1 = 1 : N1
        j1 = i1 <= L1 ? N1-L1+i1 : i1-L1
        C[i1,i2] = d[j1,j2]
    end
end
Xm = reshape(x, K1, K2)
V  = zeros(eltype(d), N1, N2)
for i2 = 1 : K2
    for i1 = 1 : K1
        V[i1,i2] = Xm[K1+1-i1, K2+1-i2]
    end
end
r1 = ifft(fft(C) .* fft(V));
r1 = r1[1:L1,1:L2]
r1 = vec(r1)
vecnorm(r1-r) / vecnorm(r)  # it works
# ==============================================================================
# test the Hermitian(H) * x
x = rand(eltype(d), L1*L2)
r = H' * x

# fast multiplication
C = zeros(eltype(d), N1, N2)
for i2 = 1 : N2
    j2 = i2 <= K2 ? N2-K2+i2 : i2-K2
    for i1 = 1 : N1
        j1 = i1 <= K1 ? N1-K1+i1 : i1-K1
        C[i1,i2] = conj(d[j1,j2])
    end
end
Xm = reshape(x, L1, L2)
V  = zeros(eltype(d), N1, N2)
for i2 = 1 : L2
    for i1 = 1 : L1
        V[i1,i2] = Xm[L1+1-i1, L2+1-i2]
    end
end
r1 = ifft(fft(C) .* fft(V));
r1 = r1[1:K1,1:K2]
r1 = vec(r1)
vecnorm(r1-r) / vecnorm(r)  # it works


# ==============================================================================
# test the multiplication of Hankel matrix with a vector and its transpose (4D) case
# ==============================================================================
N1 = 11; N2 = 12; N3 = 13; N4 = 14;
d = rand(Complex128, N1, N2, N3, N4);
K1 = floor(Int64, N1/2) + 1; L1 = N1 - K1 + 1;
K2 = floor(Int64, N2/2) + 1; L2 = N2 - K2 + 1;
K3 = floor(Int64, N3/2) + 1; L3 = N3 - K3 + 1;
K4 = floor(Int64, N4/2) + 1; L4 = N4 - K4 + 1;
H = zeros(eltype(d), L1*L2*L3*L4, K1*K2*K3*K4)

# fourth layer
for j4 = 1 : K4
    col4 = (j4-1)*K3*K2*K1
    for i4 = 1 : L4
        row4 = (i4-1)*L3*L2*L1
        m4 = i4 + j4 - 1

        # third layer
        for j3 = 1 : K3
            col3 = (j3-1)*K2*K1
            for i3 = 1 : L3
                row3 = (i3-1)*L2*L1
                m3 = i3 + j3 - 1

                # second layer
                for j2 = 1 : K2
                    col2 = (j2-1)*K1
                    for i2 = 1 : L2
                        row2 = (i2-1)*L1
                        m2 = i2 + j2 - 1

                        # first layer
                        for j1 = 1 : K1
                            col1 = col4 + col3 + col2 + j1
                            for i1 = 1 : L1
                                row1 = row4 + row3 + row2 + i1
                                m1 = i1 + j1 - 1
                                H[row1, col1] = d[m1,m2,m3,m4]
                            end
                        end #first

                    end
                end #second
            end
        end #third
    end
end #fourth

# direct multiplication
x = rand(eltype(d), K1*K2*K3*K4)
r = H * x

# fast multiplication
C = zeros(eltype(d), N1, N2, N3, N4)
for i4 = 1 : N4
    j4 = i4 <= L4 ? N4-L4+i4 : i4-L4
    for i3 = 1 : N3
        j3 = i3 <= L3 ? N3-L3+i3 : i3-L3
        for i2 = 1 : N2
            j2 = i2 <= L2 ? N2-L2+i2 : i2-L2
            for i1 = 1 : N1
                j1 = i1 <= L1 ? N1-L1+i1 : i1-L1
                C[i1,i2,i3,i4] = d[j1,j2,j3,j4]
            end
        end
    end
end
Xm = reshape(x, K1, K2, K3, K4)
V  = zeros(eltype(d), N1, N2, N3, N4)
for i4 = 1 : K4
    for i3 = 1 : K3
        for i2 = 1 : K2
             for i1 = 1 : K1
                 V[i1,i2,i3,i4] = Xm[K1+1-i1, K2+1-i2, K3+1-i3, K4+1-i4]
             end
        end
    end
end
r1 = ifft(fft(C) .* fft(V));
r1 = r1[1:L1, 1:L2, 1:L3, 1:L4]
r1 = vec(r1)
vecnorm(r1-r) / vecnorm(r)  # it works
# ==============================================================================

# test the Hermitian(H) * x
x = rand(eltype(d), L1*L2*L3*L4)
r = H' * x

# fast multiplication
C = zeros(eltype(d), N1, N2, N3, N4)
for i4 = 1 : N4
    j4 = i4 <= K4 ? N4-K4+i4 : i4-K4
    for i3 = 1 : N3
        j3 = i3 <= K3 ? N3-K3+i3 : i3-K3
        for i2 = 1 : N2
            j2 = i2 <= K2 ? N2-K2+i2 : i2-K2
            for i1 = 1 : N1
                j1 = i1 <= K1 ? N1-K1+i1 : i1-K1
                C[i1,i2,i3,i4] = conj(d[j1,j2,j3,j4])
            end
        end
    end
end

Xm = reshape(x, L1, L2, L3, L4)
V  = zeros(eltype(d), N1, N2, N3, N4)
for i4 = 1 : L4
    for i3 = 1 : L3
        for i2 = 1 : L2
            for i1 = 1 : L1
                V[i1,i2,i3,i4] = Xm[L1+1-i1, L2+1-i2, L3+1-i3, L4+1-i4]
            end
        end
    end
end
r1 = ifft(fft(C) .* fft(V));
r1 = r1[1:K1, 1:K2, 1:K3, 1:K4]
r1 = vec(r1)
vecnorm(r1-r) / vecnorm(r)  # it works

# ==============================================================================
# test anti-diagonal average (1D) case
# ==============================================================================
N = 1234;
K = floor(Int64, N/2) + 1; L = N - K + 1;
u = rand(Complex128, L); v = rand(Complex128, K); sigma = rand(Complex128;

Y = sigma * u * v';
r = zeros(Complex128, N)
for i = 1 : L
    for j = 1 : i
        r[i] = r[i] + Y[j, i-j+1]
    end
    r[i] = r[i] / i
end
for i = K : N
    for j = i-K+1 : L
        r[i] = r[i] + Y[j, i-j+1]
    end
    r[i] = r[i] / (N-i+1)
end


uhat = zeros(Complex128, N); uhat[1:L] = u[1:L];
vhat = zeros(Complex128, N); vhat[1:K] = conj(v[1:K]);

r1  = ifft(fft(uhat) .* fft(vhat));
weight = ones(N)
if iseven(N)
   weight = vcat(collect(1:L), collect(L:-1:1))
else
   weight = vcat(collect(1:L), collect(L-1:-1:1))
end
r1 = (r1 ./weight) * sigma


# ==============================================================================
# test anti-diagonal average (2D) case
# ==============================================================================
N1 = 34; N2 = 12;
K1 = floor(Int64, N1/2) + 1; L1 = N1 - K1 + 1;
K2 = floor(Int64, N2/2) + 1; L2 = N2 - K2 + 1;

sigma = rand(Complex128); u = rand(Complex128, L1*L2); v = rand(Complex128, K1*K2);
H = sigma * u * v';
d = zeros(Complex128, N1, N2); countnum = zeros(Int64, N1, N2);
for j2 = 1 : K2
    col2 = (j2-1) * K1
    for i2 = 1 : L2
        row2 = (i2-1) * L1
        idx2 = i2 + j2 - 1
        for j1 = 1 : K1
            col1 = col2 + j1
            for i1 = 1 : L1
                row1 = row2 + i1
                idx1 = i1 + j1 - 1
                countnum[idx1, idx2] = countnum[idx1, idx2] + 1
                d[idx1, idx2] = d[idx1, idx2] + H[row1, col1]
            end
        end
    end
end
for i2 = 1 : N2
    for i1 = 1 : N1
        d[i1,i2] = d[i1,i2] / countnum[i1,i2]
    end
end

uhat = zeros(Complex128, N1, N2);
vhat = zeros(Complex128, N1, N2);
uhat[1:L1, 1:L2] = reshape(u, L1, L2);
vhat[1:K1, 1:K2] = reshape(conj(v), K1, K2);
d1 = ifft(fft(uhat).*fft(vhat));
for i2 = 1 : N2
    for i1 = 1 : N1
        d1[i1,i2] = d1[i1,i2] / countnum[i1,i2] * sigma
    end
end


# ==============================================================================
# test anti-diagonal average (4D) case
# ==============================================================================
N1 = 17; N2 = 12; N3 = 16; N4 = 18;
K1 = floor(Int64, N1/2) + 1; L1 = N1 - K1 + 1;
K2 = floor(Int64, N2/2) + 1; L2 = N2 - K2 + 1;
K3 = floor(Int64, N3/2) + 1; L3 = N3 - K3 + 1;
K4 = floor(Int64, N4/2) + 1; L4 = N4 - K4 + 1;

sigma = rand(Complex128); u = rand(Complex128, L1*L2*L3*L4); v = rand(Complex128, K1*K2*K3*K4);
H = sigma * u * v';
d = zeros(Complex128, N1, N2, N3, N4); countnum = zeros(Int64, N1, N2, N3, N4);
for j4 = 1 : K4
    col4 = (j4-1)*K3*K2*K1    # column index for hankel matrix
    for i4 = 1 : L4
        row4 = (i4-1)*L3*L2*L1    # row index for hankel matrix
        idx4 = i4 + j4 - 1
        for j3 = 1 : K3
            col3 = (j3-1)*K2*K1
            for i3 = 1 : L3
                row3 = (i3-1)*L2*L1
                idx3 = i3 + j3 - 1
                for j2 = 1 : K2
                    col2 = (j2-1) * K1
                    for i2 = 1 : L2
                        row2 = (i2-1) * L1
                        idx2 = i2 + j2 - 1
                        for j1 = 1 : K1
                            col1 = col4 + col3 + col2 + j1
                            for i1 = 1 : L1
                                row1 = row4 + row3 + row2 + i1
                                idx1 = i1 + j1 - 1
                                countnum[idx1, idx2, idx3, idx4] = countnum[idx1, idx2, idx3, idx4] + 1
                                d[idx1, idx2, idx3, idx4] = d[idx1, idx2, idx3, idx4] + H[row1, col1]
                            end
                        end
                    end
                end
            end
        end
    end
end

for i4 = 1 : N4
    for i3 = 1 : N3
        for i2 = 1 : N2
            for i1 = 1 : N1
                d[i1,i2,i3,i4] = d[i1,i2,i3,i4] / countnum[i1,i2,i3,i4]
            end
        end
    end
end


uhat = zeros(Complex128, N1, N2, N3, N4);
vhat = zeros(Complex128, N1, N2, N3, N4);
uhat[1:L1, 1:L2, 1:L3, 1:L4] = reshape(u, L1, L2, L3, L4);
vhat[1:K1, 1:K2, 1:K3, 1:K4] = reshape(conj(v), K1, K2, K3, K4);
d1 = ifft(fft(uhat).*fft(vhat));
for i4 = 1 : N4
    for i3 = 1 : N3
        for i2 = 1 : N2
            for i1 = 1 : N1
                d1[i1,i2,i3,i4] = d1[i1,i2,i3,i4] / countnum[i1,i2,i3,i4] * sigma
            end
        end
    end
end

# a better way to compute the countnum
# only applicable to the L=K or L=K-1, It only need to compute once for all the frequencies
countnum1 = ones(Int64, N1, N2, N3, N4)
for i4 = 1 : N4
    c4 = i4 <= L4 ? i4 : N4-i4+1
    for i3 = 1 : N3
        c3 = i3 <= L3 ? i3 : N3-i3+1
        for i2 = 1 : N2
            c2 = i2 <= L2 ? i2 : N2-i2+1
            for i1 = 1 : N1
                c1 = i1 <= L1 ? i1 : N1-i1+1
                countnum1[i1,i2,i3,i4] = c4*c3*c2*c1
            end
        end
    end
end
vecnorm(countnum - countnum1)

# ==============================================================================
#                 test the properties of kronnecker product
# ==============================================================================
m1 = 33; n1 = 14;
m2 = 13; n2 = 17;
m3 = 21; n3 = 19;
M  = m1 * m2 * m3;
N  = n1 * n2 * n3;

A1 = rand(Complex128, m1, n1);
A2 = rand(Complex128, m2, n2);
A3 = rand(Complex128, m3, n3);

A = kron(A3, kron(A2, A1));
x = rand(Complex128, N)
y = A * x
y = reshape(y, m1, m2, m3);

x1 = reshape(x, n1, n2, n3);
y1 = zeros(Complex128, m1, n2, n3);
y2 = zeros(Complex128, m1, m2, n3);
y3 = zeros(Complex128, m1, m2, m3);

# process first dimension
for i3 = 1 : n3
    for i2 = 1 : n2
        y1[:,i2,i3] = A1 * x1[:,i2,i3]
    end
end

for i3 = 1 : n3
    for i1 = 1 : m1
        y2[i1,:,i3] = A2 * vec(y1[i1,:,i3])
    end
end

for i2 = 1 : m2
    for i1 = 1 : m1
        y3[i1,i2,:] = A3 * vec(y2[i1,i2,:])
    end
end

vecnorm(y3-y) / vecnorm(y);
# extended to high-dimensional separable transform
# ==============================================================================
#                      test the properties of kronnecker product 2D case
# ==============================================================================
m1 = 33; n1 = 14;
m2 = 13; n2 = 17;

M  = m1 * m2
N  = n1 * n2

A1 = randn(m1, n1);
A2 = randn(m2, n2);

A = kron(A2, A1);
x = randn(N)
y = A * x
y = reshape(y, m1, m2);

x1 = reshape(x, n1, n2, n3);
y1 = zeros(m1, n2);
y2 = zeros(m1, m2);
# process first dimension
    for i2 = 1 : n2
        y1[:,i2] = A1 * x1[:,i2]
    end


    for i1 = 1 : m1
        y2[i1,:,i3] = A2 * vec(y1[i1,:,i3])
    end
end

for i2 = 1 : m2
    for i1 = 1 : m1
        y3[i1,i2,:] = A3 * vec(y2[i1,i2,:])
    end
end

# ==============================================================================
#                    the comumutative properties of kronnecker product
# ==============================================================================
# create stride permutation matrix
function stride_permutation(n, p)
    m = round(Int64, n/p)
    P = zeros(n, n)
    for i2 = 1 : p
        for i1 = 1 : m
            row_idx = (i2-1)*m + i1
            col_idx = (i1-1)*p + i2
            P[row_idx, col_idx] = 1
        end
    end
    return P
end

m1 = 33; n1 = 48;
m2 = 49; n2 = 17;

A = rand(m1,n1);
B = rand(m2,n2);

M = kron(A, B);
T = kron(B, A);

P1 = stride_permutation(m1*m2, m1)
P2 = stride_permutation(n1*n2, n2)

M1 = P1 * T * P2
vecnorm(M-M1)
