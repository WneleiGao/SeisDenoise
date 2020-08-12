using SeisPlot, SeisDenoise, SeisProcessing

# time sample
dt = 8.0/1000.

# synthetic data
d = get_mixture_events(nt=250, nx1=120, nx2=120, dt=dt, dx1=10, dx2=10,
                              v1=[6000. ,2000,-2700.,-2700.,6000. ],
                              v2=[14000.,4000, 4000., 4000.,14000.],
                              tau=[0.1  ,0.4 ,  0.9 ,  1.2 ,1.4   ],
                              amp=[1.0  ,0.9 ,  1.0 , -1.0 ,0.6   ],
                              apx1=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              apx2=[0.0 ,0.0 ,  0.0 ,  0.0 ,0.0   ],
                              event_type=['l', 'h', 'p', 'p', 'l' ],
                              f0=20.0, static_flag=true, L=13, sigma=0.03);

# add random noise to data
dn = SeisAddNoise(d, 0.75);

# work directory
work_dir = joinpath(homedir(), "Desktop/patches");
isdir(work_dir) && rm(work_dir, recursive=true);
mkdir(work_dir);

# fxy eigen
s1 = local_fxy_eigen(dn, 4, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                     x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# mssa
s2 = local_mssa(dn, 2, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10);

# fxy prediction
s3 = local_fxy_prediction(dn, 2, work_dir; dt=dt, flow=2.0, fhigh=60.0,
                          x1_wl = 30, x1_wo = 10, x2_wl = 30, x2_wo = 10,
                          max_iter=10, mu=1.0e-6, tol=1.0e-9);

# plotting the result
i = 10;
tmp1 = hcat(s1[:,:,i], dn[:,:,i]-s1[:,:,i]);
tmp2 = hcat(s2[:,:,i], dn[:,:,i]-s2[:,:,i]);
tmp3 = hcat(s3[:,:,i], dn[:,:,i]-s3[:,:,i]);
tmp4 = hcat( d[:,:,i], dn[:,:,i]- d[:,:,i]);
SeisPlotTX(tmp1, cmap="gray", wbox=10, title="eigen")
SeisPlotTX(tmp2, cmap="gray", wbox=10, title="mssa")
SeisPlotTX(tmp3, cmap="gray", wbox=10, title="fxy_prediction")
SeisPlotTX(tmp4, cmap="gray", wbox=10, title="clean")
