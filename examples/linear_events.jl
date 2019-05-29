# testing on
using  SeisPlot, FFTW, LinearAlgebra, SeisProcessing

dt = 4.0/1000.
d = SeisLinearEvents(;ot=0.0,dt=dt, nt=200, ox1=0.0, dx1=10.0,
                      nx1=80, ox2=0.0, dx2=10.0, nx2=83, ox3=0.0, dx3=10.0,
                      nx3=1, ox4=0.0, dx4=10.0, nx4=1, tau=[.4,.5],
                      p1=[0.0001,-0.0003],p2=[0.0001,-0.0003],p3=[0.,0],p4=[0.,0.],
                      amp=[1.0,-1.0], f0=24.0);

dn = SeisAddNoise(d, .5);
sn = fxy_prediction(dn; dt=dt, L=4, flow=0.0, fhigh=60.0,
                    Niter=10, mu=1.0e-6, tol=1.0e-9);

T  = tucker_als(Tensor(dn), 3);
st = tucker2tensor(T);

SeisPlotTX(d[:,:,10], style="wiggles", xcur=2)
SeisPlotTX(dn[:,:,10], style="wiggles", xcur=2)
SeisPlotTX(sn[:,:,10], style="wiggles", xcur=2)
SeisPlotTX(st[:,:,10], style="wiggles", xcur=2)
SeisPlotTX(dn[:,:,10] - sn[:,:,10], style="wiggles", xcur=2)
