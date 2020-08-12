"""
   h is the window half-length in vertical direction and w is the half-length in
horizontal direction.
"""
function scalar_median_filter(d::Matrix{Tv}, h::Ti, w::Ti) where {Tv<:Real, Ti<:Int64}

    (n1, n2) = size(d)
    d_filter = zeros(Tv, n1, n2)

    for i2 = 1 : n2
        for i1 = 1 : n1

            # vertical range
            i1_lower = i1 - h < 1  ? 1  : i1 - h
            i1_upper = i1 + h > n1 ? n1 : i1 + h

            # horizontal range
            i2_lower = i2 - w < 1  ? 1  : i2 - w
            i2_upper = i2 + w > n2 ? n2 : i2 + w

            v = median(view(d, i1_lower:i1_upper, i2_lower:i2_upper))
            d_filter[i1,i2] = v
        end
    end
    return d_filter
end

function vector_median_filter(d::Matrix{Tv}, h::Ti, w::Ti) where {Tv<:Real, Ti<:Int64}



end



dcut = d[281:449,20,77:128];
SeisPlotTX(dcut, style="wiggles", xcur=3.0);

d1 = scalar_median_filter(dcut, 2, 2);
SeisPlotTX(d1, style="wiggles", xcur=3.0);


dir_work = joinpath(homedir(), "Desktop/PGS/PGS_rsf");

path_d   = joinpath(dir_work, "data.rsf");
(hdr, d) = read_RSdata(path_d);
for i = 40 : 40 : size(d,2)
    SeisPlotTX(d[:,i,:], hbox=8, wbox=5, pclip=90, cmap="gray"); tight_layout();
end

idx_channel = 100; dt = 0.016;

d1 = scalar_median_filter(d[:,idx_channel,:], 1, 1);

i1l = 1; i1u = size(d,1); ot = (i1l-1)*dt; tmax = (i1u-1)*dt; i1_range = range(ot, length=5, stop=tmax);
i2l = 1; i2u = size(d,3); i2_range = range(i2l, length=5, stop=i2u);
SeisPlotTX(d[i1l:i1u,idx_channel,i2l:i2u], hbox=8, wbox=6, pclip=95, cmap="gray",
           dx=1, dy=dt, ox=i2l, oy=ot, xticks=i2_range, yticks=i1_range, ticksize=15,
           xlabel="Trace number", ylabel="Time (S)", labelsize=15); tight_layout();


SeisPlotTX(d1, hbox=8, wbox=6, pclip=95, cmap="gray",
           dx=1, dy=dt, ox=i2l, oy=ot, xticks=i2_range, yticks=i1_range, ticksize=15,
           xlabel="Trace number", ylabel="Time (S)", labelsize=15); tight_layout();
