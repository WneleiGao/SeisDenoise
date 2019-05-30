using SeisPlot, SeisProcessing, SeisDenoise

dt = 4.0/1000.

# generate data with two linear events
d = SeisLinearEvents(;ot=0.0,dt=dt, nt=200, ox1=0.0, dx1=10.0,
                      nx1=40, ox2=0.0, dx2=10.0, nx2=40, ox3=0.0, dx3=10.0,
                      nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[.4,.5],
                      p1=[0.0001,-0.0003],p2=[0.0001,-0.0003],p3=[0.,0],p4=[0.,0.],
                      amp=[1.0,-1.0], f0=24.0);

# add random noise to data
dn = SeisAddNoise(d, .5);

# fxy prediction
s1 = fxy_prediction(dn, 2; dt=dt, flow=0.0, fhigh=60.0,
                    Niter=10, mu=1.0e-6, tol=1.0e-9);

# MSSA
s2 = mssa(dn, 2; dt=dt, flow=0.0, fhigh=60.0);

# fxy eigen
s3 = fxy_eigen(dn, 2; dt=dt, flow=0.0, fhigh=60.0);

# plotting the result
i = 10;
SeisPlotTX(d[:,:,i], style="wiggles", xcur=2)
SeisPlotTX(dn[:,:,i], style="wiggles", xcur=2)

tmp1 = hcat(s1[:,:,i], dn[:,:,10]-s1[:,:,10]);
tmp2 = hcat(s2[:,:,i], dn[:,:,10]-s2[:,:,10]);
tmp3 = hcat(s3[:,:,i], dn[:,:,10]-s3[:,:,10]);
SeisPlotTX(tmp1, style="wiggles", xcur=2, wbox=10)
SeisPlotTX(tmp2, style="wiggles", xcur=2, wbox=10)
SeisPlotTX(tmp3, style="wiggles", xcur=2, wbox=10)
