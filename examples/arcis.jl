function read_header(path_header::String)

    # open the file
    fid = open(path_header, "r")

    # number samples
    nt = read(fid, Int32)

    # number of gather
    num_traces = read(fid, Int32)

    # return four vector of inline, crossline, offset, azimuth
    inline    = zeros(Int32, num_traces); read!(fid, inline   );
    crossline = zeros(Int32, num_traces); read!(fid, crossline);
    offset    = zeros(Int32, num_traces); read!(fid, offset   );
    azimuth   = zeros(Int32, num_traces); read!(fid, azimuth  );
    close(fid)

    return nt, inline, crossline, offset, azimuth
end

function read_fold_map(path::String)

   fid  = open(path, "r")
   nt = read(fid, Int32)

   # size of fold map
   num_crline = read(fid, Int32)
   num_inline = read(fid, Int32)

   # read a vector of Int32
   fold = zeros(Int32, num_crline * num_inline)
   read!(fid, fold)
   close(fid)

   # return it as a matrix
   fold = convert(Vector{Int64}, fold)
   return reshape(fold, Int64(num_crline), Int64(num_inline))

end

function estimate_fold_number(inline::Vector{Ti}, crline::Vector{Ti}) where {Ti<:Integer}

    ntraces = length(inline)

    # get number of CDP gathers
    num_gathers = 1
    for i = 2 : ntraces
        if inline[i] != inline[i-1] || crline[i] != crline[i-1]
           num_gathers = num_gathers + 1
        end
    end

    # a vector fold save the number of traces for each gather
    fold = zeros(Int64, num_gathers)
    fold[1] = 1
    k       = 1

    # loop over traces
    for i = 2 : ntraces
        if inline[i] == inline[i-1] && crline[i] == crline[i-1]
           fold[k] += 1
        else
           k = k + 1
           fold[k] += 1
        end
    end

    return fold
end

function stack_raw_data(path_fold::Ts, path_data::Ts) where {Ts<:String}

    # read the foldmap
    fold = read(path_fold)

        dout = zeros(Float32, nt, numcr, numin)
    fdata = open(pd, "r")
    esize = 4
    for j = 1 : numin
        for i = 1 : numcr
            ntrace = fold[i,j]
            tmp = read(fdata, Float32, nt*ntrace); tmp = reshape(tmp, Int64(nt), Int64(ntrace));
            dout[:,i,j] = squeeze(sum(tmp,2), 2) / fold[i,j]
        end
    end
    close(fdata);
    fid = open(pout, "w");
    write(fid, Int32(nt), Int32(numcr), Int32(numin));
    write(fid, vec(convert(Array{Float32}, dout))); close(fid)
    return dout
end

#==============================================================================#
#                                    PS wave
#==============================================================================#
# read the header information
path_ps_header = joinpath(homedir(), "Desktop/PS_wave/PSHeader.bin");
(nt, inline_ps, crline_ps, offset_ps, azimuth_ps) = read_header(path_ps_header);

# the dimension of the cube
inline_lower_ps = minimum(inline_ps)
inline_upper_ps = maximum(inline_ps)
crline_lower_ps = minimum(crline_ps)
crline_upper_ps = maximum(crline_ps)

# # get the fold number of each gather
# fold_ps = estimate_fold_number(inline_ps, crline_ps)

# read fold map
path_fold = joinpath(homedir(), "Desktop/PS_wave/foldmap.bin");
fold_ps   = read_fold_map(path_fold)

# the data file
path_ps_data = joinpath(homedir(), "Desktop/PS_wave/PSData.bin");

#==============================================================================#
#                                    PP wave
#==============================================================================#
path_pp_header = joinpath(homedir(), "Desktop/PP_wave/PPHeader.bin");
(nt, inline_pp, crline_pp, offset_pp, azimuth_pp) = read_header(path_pp_header);

inline_lower_ps = minimum(inline_pp)
inline_upper_ps = maximum(inline_pp)
crline_lower_ps = minimum(crline_pp)
crline_upper_ps = maximum(crline_pp)

# # get the fold number of each gather
# fold_pp = estimate_fold_number(inline_pp, crline_pp)

# read fold map
path_fold = joinpath(homedir(), "Desktop/PP_wave/foldmap.bin");
fold_pp   = read_fold_map(path_fold)

path_pp_data = joinpath(homedir(), "Desktop/PP_wave/PPData.bin");
