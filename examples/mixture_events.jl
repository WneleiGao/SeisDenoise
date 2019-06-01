using SeisPlot, SeisDenoise, SeisProcessing

dt = 8.0/1000.

d = get_mixture_events(nt=250, nx1=200, nx2=200, dt=dt, dx1=10, dx2=10,
                              v1=[6000. ,2000,-2700.,-2700.,6000. ],
                              v2=[14000.,4000, 4000., 4000.,14000.],
                              tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                              amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                              apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              event_type=['l', 'h', 'p', 'p', 'l' ],
                              f0=20.0, static_flag=true, L=13, sigma=0.1);

# add random noise to data
dn = SeisAddNoise(d, 0.75);

# work directory
work_dir = joinpath(homedir(), "Desktop/randCP");

# fxy eigen
s1 = local_fxy_eigen(dn, 4, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                     x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# mssa
s2 = local_mssa(dn, 18, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# fxy prediction
s3 = local_fxy_prediction(dn, 2, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                          x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10,
                          max_iter=10, mu=1.0e-6, tol=1.0e-9);

# tensor als
s4 = local_tucker_als(dn, 9, work_dir; tol=1.0e-6, max_iter=10, init_flag="random",
                      it_wl = 30, it_wo = 10, x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);


snr0 = MeasureSNR(d, dn, db=true)
snr1 = MeasureSNR(d, s1, db=true)
snr2 = MeasureSNR(d, s2, db=true)
snr3 = MeasureSNR(d, s3, db=true)
snr4 = MeasureSNR(d, s4, db=true)

path="/Users/wenlei/Desktop/dcube.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(d ))); close(fid);
path="/Users/wenlei/Desktop/ncube.bin"; fid=open(path,"w"); write(fid,convert(Vector{Float32}, vec(dn))); close(fid);

pscube < /Users/wenlei/Desktop/dcube.bin size1=2.7 size2=2.5 size3=1.2 labelsize=13 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="X" n3=200 d3num=40 f3num=40 label3="Y" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/d.eps
pscube < /Users/wenlei/Desktop/ncube.bin size1=2.7 size2=2.5 size3=1.2 labelsize=13 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=200 d2num=40 f2num=40 label2="X" n3=200 d3num=40 f3num=40 label3="Y" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/n.eps
# psimage < /Users/wenlei/Desktop/msecdiff.bin height=3.0 width=2.0 labelsize=12 n1=250 d1=0.008 d1num=0.4 f1num=0.4 label1="Time (s)" n2=60 d2num=40 f2num=40 label2="Crossline" xbox=0.5 ybox=0.5 bclip=0.96 wclip=-0.96 > /Users/wenlei/Desktop/msecdiff.eps

# plotting the result
i = 55;
SeisPlotTX( d[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/d.pdf"); close();
SeisPlotTX(dn[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/dn.pdf"); close();
SeisPlotTX(dn[:,:,i]-d[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n.pdf"); close();

SeisPlotTX(s1[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/s1.pdf"); close();
SeisPlotTX(dn[:,:,i]-s1[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n1.pdf"); close();

SeisPlotTX(s2[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/s2.pdf"); close();
SeisPlotTX(dn[:,:,i]-s2[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n2.pdf"); close();

SeisPlotTX(s3[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/s3.pdf"); close();
SeisPlotTX(dn[:,:,i]-s3[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n3.pdf"); close();

SeisPlotTX(s4[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/s4.pdf"); close();
SeisPlotTX(dn[:,:,i]-s4[:,:,i], cmap="gray", wbox=3, hbox=2.5, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n4.pdf"); close();



# plot in gray scale
tmp1 = hcat(s1[:,:,i], dn[:,:,i]-s1[:,:,i]);
tmp2 = hcat(s2[:,:,i], dn[:,:,i]-s2[:,:,i]);
tmp3 = hcat(s3[:,:,i], dn[:,:,i]-s3[:,:,i]);
tmp4 = hcat(s4[:,:,i], dn[:,:,i]-s4[:,:,i]);
tmp5 = hcat(dn[:,:,i], dn[:,:,i]- d[:,:,i]);
SeisPlotTX(tmp1, wbox=10, cmap="gray");
SeisPlotTX(tmp2, wbox=10, cmap="gray");
SeisPlotTX(tmp3, wbox=10, cmap="gray");
SeisPlotTX(tmp4, wbox=10, cmap="gray");
SeisPlotTX(tmp5, wbox=10, cmap="gray");

# wiggles
itl = 90; itu=171;
ixl = 77; ixu=137;
tmp1 = hcat(s1[itl:itu,ixl:ixu,i], dn[itl:itu,ixl:ixu,i]-s1[itl:itu,ixl:ixu,i]);
tmp2 = hcat(s2[itl:itu,ixl:ixu,i], dn[itl:itu,ixl:ixu,i]-s2[itl:itu,ixl:ixu,i]);
tmp3 = hcat(s3[itl:itu,ixl:ixu,i], dn[itl:itu,ixl:ixu,i]-s3[itl:itu,ixl:ixu,i]);
tmp4 = hcat(s4[itl:itu,ixl:ixu,i], dn[itl:itu,ixl:ixu,i]-s4[itl:itu,ixl:ixu,i]);
tmp5 = hcat( d[itl:itu,ixl:ixu,i], dn[itl:itu,ixl:ixu,i]- d[itl:itu,ixl:ixu,i]);

xcur = 2.0
SeisPlotTX(tmp1, style="wiggles", xcur=xcur, wbox=10, wiggle_trace_increment=1)
SeisPlotTX(tmp2, style="wiggles", xcur=xcur, wbox=10, wiggle_trace_increment=1)
SeisPlotTX(tmp3, style="wiggles", xcur=xcur, wbox=10, wiggle_trace_increment=1)
SeisPlotTX(tmp4, style="wiggles", xcur=xcur, wbox=10, wiggle_trace_increment=1)
SeisPlotTX(tmp5, style="wiggles", xcur=xcur, wbox=10, wiggle_trace_increment=1)

SeisPlotTX( d[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/d.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/dn.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i]-d[itl:itu,ixl:ixu,i], xcur=1.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/n.pdf"); close();

SeisPlotTX(s1[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/s1.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i]-s1[itl:itu,ixl:ixu,i], xcur=1.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/n1.pdf"); close();

SeisPlotTX(s2[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/s2.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i]-s2[itl:itu,ixl:ixu,i], xcur=1.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/n2.pdf"); close();

SeisPlotTX(s3[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/s3.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i]-s3[itl:itu,ixl:ixu,i], xcur=1.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/n3.pdf"); close();

SeisPlotTX(s4[itl:itu,ixl:ixu,i], xcur=2.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/s4.pdf"); close();
SeisPlotTX(dn[itl:itu,ixl:ixu,i]-s4[itl:itu,ixl:ixu,i], xcur=1.0, style="wiggles", wbox=3, hbox=2.5, dy=0.008, oy=0.72, ox=77, ylabel="Time (s)", xlabel="X", xticks=collect(80:10:130), yticks=collect(0.8:0.2:1.3), wiggle_trace_increment=2); tight_layout(); savefig("/Users/wenlei/Desktop/n4.pdf"); close();



t = 111;
tmp1 = hcat(s1[t,:,:], dn[t,:,:]-s1[t,:,:]);
tmp2 = hcat(s2[t,:,:], dn[t,:,:]-s2[t,:,:]);
tmp3 = hcat(s3[t,:,:], dn[t,:,:]-s3[t,:,:]);
tmp4 = hcat(s4[t,:,:], dn[t,:,:]-s4[t,:,:]);
tmp5 = hcat( d[t,:,:], dn[t,:,:]- d[t,:,:]);

# plot in gray scale
SeisPlotTX(tmp1, wbox=10, cmap="gray");
SeisPlotTX(tmp2, wbox=10, cmap="gray");
SeisPlotTX(tmp3, wbox=10, cmap="gray");
SeisPlotTX(tmp4, wbox=10, cmap="gray");
SeisPlotTX(tmp5, wbox=10, cmap="gray");

i = 55; j = 77; t = 100;
SeisPlotTX( d[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/d1.pdf"); close();
SeisPlotTX(dn[:,:,i], cmap="gray", wbox=3, hbox=4, dy=0.008, ylabel="Time (s)", xlabel="X", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n1.pdf"); close();
SeisPlotTX( d[:,j,:], cmap="gray", wbox=3, hbox=4, dy=0.008, ylabel="Time (s)", xlabel="Y", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/d2.pdf"); close();
SeisPlotTX(dn[:,j,:], cmap="gray", wbox=3, hbox=4, dy=0.008, ylabel="Time (s)", xlabel="Y", xticks=collect(40:40:160), yticks=collect(0.4:0.4:1.6)); tight_layout(); savefig("/Users/wenlei/Desktop/n2.pdf"); close();
SeisPlotTX( d[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", xticks=collect(40:40:160), yticks=collect(40:40:160)); tight_layout(); savefig("/Users/wenlei/Desktop/d3.pdf"); close();
SeisPlotTX(dn[t,:,:], cmap="gray", wbox=3, hbox=4, ylabel="Y", xlabel="X", xticks=collect(40:40:160), yticks=collect(40:40:160)); tight_layout(); savefig("/Users/wenlei/Desktop/n3.pdf"); close();


SeisPlotTX(d[:,:,i], style="wiggles", xcur=xcur, wiggle_trace_increment=3)
SeisPlotTX(dn[:,:,i], style="wiggles", xcur=xcur, wiggle_trace_increment=3)
