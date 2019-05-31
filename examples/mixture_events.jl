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
dn = SeisAddNoise(d, 1.0);

# work directory
work_dir = "/Users/wenyue/Desktop/randCP"

# tucker decomposition
# s1 = tucker_als(Tensor(dn), 10; init_flag="eigvec");
# s1 = tucker2tensor(s1).d;

# fxy eigen
s1 = local_fxy_eigen(dn, 5, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                     x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# mssa
s2 = local_mssa(dn, 5, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# fxy prediction
s3 = local_fxy_prediction(dn, 2, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                          x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10,
                          max_iter=10, mu=1.0e-6, tol=1.0e-9);

# tensor als
s4 = local_tucker_als(dn, 5, work_dir; tol=1.0e-6, max_iter=10, init_flag="random",
                      it_wl = 30, it_wo = 10, x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# plotting the result
i = 10;
tmp1 = hcat(s1[:,:,i], dn[:,:,10]-s1[:,:,10]);
tmp2 = hcat(s2[:,:,i], dn[:,:,10]-s2[:,:,10]);
tmp3 = hcat(s3[:,:,i], dn[:,:,10]-s3[:,:,10]);
tmp4 = hcat(s4[:,:,i], dn[:,:,10]-s4[:,:,10]);
SeisPlotTX(tmp1, style="wiggles", xcur=1.2, wbox=10, wiggle_trace_increment=5)
SeisPlotTX(tmp2, style="wiggles", xcur=1.2, wbox=10, wiggle_trace_increment=5)
SeisPlotTX(tmp3, style="wiggles", xcur=1.2, wbox=10, wiggle_trace_increment=5)
SeisPlotTX(tmp4, style="wiggles", xcur=1.2, wbox=10, wiggle_trace_increment=5)


i = 100
SeisPlotTX(d[:,:,i], style="wiggles", xcur=2, wiggle_trace_increment=3)
SeisPlotTX(dn[:,:,i], style="wiggles", xcur=2, wiggle_trace_increment=3)
