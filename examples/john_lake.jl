using SeisMain, SeisPlot, SeisDenoise, SeisProcessing

# read the binary file
path = joinpath(homedir(), "Desktop/Channel3D/binary/cube.bin");
fid = open(path, "r");
nt = read(fid, Int64);
n1 = read(fid, Int64);
n2 = read(fid, Int64);
d  = zeros(Float32, nt*n1*n2);
read!(fid, d); close(fid);
d  = reshape(d, nt, n1, n2);

# time sampling rate
dt = 0.004

# add noise to the real data
dn = SeisAddNoise(d, 0.75);

# plotting the real data
i = 55;
SeisPlotTX( d[:,:,i], cmap="gray", wbox=5, hbox=6, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.5:0.5:2.5));
SeisPlotTX(dn[:,:,i], cmap="gray", wbox=5, hbox=6, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.5:0.5:2.5));

path="/Users/wenlei/Desktop/dc.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(d ))); close(fid);
path="/Users/wenlei/Desktop/nc.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(dn))); close(fid);

pscube < /Users/wenlei/Desktop/dc.bin size1=2.7 size2=2.5 size3=1.2 labelsize=13 n1=751 d1=0.004 d1num=0.7 f1num=0.7 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="X" n3=115 d3num=30 f3num=30 label3="Y" xbox=0.5 ybox=0.5 bclip=-12979 wclip=12979 > /Users/wenlei/Desktop/dc.eps
pscube < /Users/wenlei/Desktop/nc.bin size1=2.7 size2=2.5 size3=1.2 labelsize=13 n1=751 d1=0.004 d1num=0.7 f1num=0.7 label1="Time (s)" n2=180 d2num=40 f2num=40 label2="X" n3=115 d3num=30 f3num=30 label3="Y" xbox=0.5 ybox=0.5 bclip=-12979 wclip=12979 > /Users/wenlei/Desktop/nc.eps

# slices
t = 206; i = 90; j = 55;
SeisPlotTX( d[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/d1.pdf"); close();
SeisPlotTX(dn[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/n1.pdf"); close();
SeisPlotTX( d[:,j,:], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(20:20:100), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/d2.pdf"); close();
SeisPlotTX(dn[:,j,:], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(20:20:100), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/n2.pdf"); close();
SeisPlotTX( d[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/d3.pdf"); close();
SeisPlotTX(dn[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/n3.pdf"); close();

# frequency spectra
(fx, amp) = amplitude_spectra(d[:,1:100,1], dt)

# processing
# work directory
work_dir = joinpath(homedir(), "Desktop/randCP");

# mssa
s1 = local_mssa(dn, 12, work_dir; dt=dt, flow=3.0, fhigh=60.0,
                x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# fxy prediction
s2 = local_fxy_prediction(dn, 3, work_dir; dt=dt, flow=3.0, fhigh=60.0,
                          x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10,
                          max_iter=10, mu=1.0e-6, tol=1.0e-9);

# tensor als
s3 = local_tucker_als(dn, 9, work_dir; tol=1.0e-6, max_iter=10, init_flag="random",
                      it_wl = 30, it_wo = 10, x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# compare the three result
t = 206; i = 90; j = 55;

# frontal slice
df = dn - s1;
SeisPlotTX(s1[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s1.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f1.pdf"); close();

df = dn - s2;
SeisPlotTX(s2[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s2.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f2.pdf"); close();

df = dn - s3;
SeisPlotTX(s3[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s3.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f3.pdf"); close();



# lateral slice
df = dn - s1;
SeisPlotTX(dn[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/dn.pdf"); close();
SeisPlotTX(s1[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s1.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f1.pdf"); close();

df = dn - s2;
SeisPlotTX(s2[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(20:20:100), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s2.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f2.pdf"); close();

df = dn - s3;
SeisPlotTX(s3[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="Y", xticks=collect(20:20:100), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/s3.pdf"); close();
SeisPlotTX(df[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.004, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.7:0.7:2.8)); tight_layout(); savefig("/Users/wenlei/Desktop/f3.pdf"); close();





# time slice
df = dn - s1;
SeisPlotTX( d[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/d.pdf") ; close();
SeisPlotTX(dn[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/dn.pdf"); close();
SeisPlotTX(s1[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/s1.pdf"); close();
SeisPlotTX(df[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/f1.pdf"); close();

df = dn - s2;
SeisPlotTX(s2[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/s2.pdf"); close();
SeisPlotTX(df[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/f2.pdf"); close();

df = dn - s3;
SeisPlotTX(s3[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/s3.pdf"); close();
SeisPlotTX(df[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", yticks=collect(40:40:160), xticks=collect(20:20:100)); tight_layout(); savefig("/Users/wenlei/Desktop/f3.pdf"); close();







# read the data
# path = joinpath(homedir(), "Desktop/Channel3D/seis/john_lake")
# (d, h, ext) = SeisRead(path);
#
# # reshape into a cube
# nt = 751; n1 = 225; n2 = 217;
# d = reshape(d, nt, n1, n2);
#
# # cute part of the data
# ix_1l = 11; ix_1u = 190;
# ix_2l = 71; ix_2u = 185;
# c = d[:, ix_1l:ix_1u, ix_2l:ix_2u];
#
# # save it as a binary file
# (nt, n1, n2) = size(c);
#
# path = joinpath(homedir(), "Desktop/Channel3D/binary/cube.bin");
# fid = open(path, "w");
# write(fid, nt, n1, n2);
# write(fid, vec(c));
# close(fid);
