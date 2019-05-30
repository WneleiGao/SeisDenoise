using SeisPlot, SeisDenoise

dt = 8.0/1000.

d = get_mixture_events(nt=250, nx1=200, nx2=200, dt=dt, dx1=10, dx2=10,
                              v1=[6000. ,2000,-2700.,-2700.,6000. ],
                              v2=[14000.,4000, 4000., 4000.,14000.],
                              tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                              amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                              apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              event_type=['l', 'h', 'p', 'p', 'l' ],
                              f0=20.0, static_flag=true, L=13, sigma=0.1)

# add random noise to data
dn = SeisAddNoise(d, 1.0);

# tucker decomposition
s1 = tucker_als(Tensor(dn), 10; init_flag="eigvec");
s1 = tucker2tensor(s1).d;

# MSSA
s2 = mssa(dn, 2; dt=dt, flow=0.0, fhigh=60.0);

# fxy prediction
s3 = fxy_prediction(dn, 5; dt=dt, flow=0.0, fhigh=60.0,
                    Niter=10, mu=1.0e-6, tol=1.0e-9);

# fxy eigen
s4 = fxy_eigen(dn, 2; dt=dt, flow=0.0, fhigh=60.0);

# plotting the result
i = 1;
tmp1 = hcat(s1[:,:,i], dn[:,:,10]-s1[:,:,10]);
tmp2 = hcat(s2[:,:,i], dn[:,:,10]-s2[:,:,10]);
tmp3 = hcat(s3[:,:,i], dn[:,:,10]-s3[:,:,10]);
tmp4 = hcat(s4[:,:,i], dn[:,:,10]-s4[:,:,10]);
SeisPlotTX(tmp1, style="wiggles", xcur=2, wbox=10)
SeisPlotTX(tmp2, style="wiggles", xcur=2, wbox=10)
SeisPlotTX(tmp3, style="wiggles", xcur=2, wbox=10, wiggle_trace_increment=3)
SeisPlotTX(tmp4, style="wiggles", xcur=2, wbox=10)

i = 100
SeisPlotTX(d[:,:,i], style="wiggles", xcur=2, wiggle_trace_increment=3)
SeisPlotTX(dn[:,:,i], style="wiggles", xcur=2, wiggle_trace_increment=3)
