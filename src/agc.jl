"""
   compute the envelope of seismic data
"""
function envelope(d::Array{Tv}) where {Tv<:AbstractFloat}

    # the dimensions of input
    dims = size(d)

    # matrix
    if length(dims) == 2
       dout = zeros(Tv, dims[1], dims[2])
       for i = 1 : dims[2]
           tmp = hilbert(d[:,i])
           dout[:,i] .= abs.(tmp)
       end

    # cube
    elseif length(dim) == 3
       dout = zeros(Tv, dims[1], dims[2], dims[3])
       for i3 = 1 : dims[3]
           for i2 = 1 : dims[2]
               tmp = hilbert(d[:,i2,i3])
               dout[:,i2,i3] .= abs.(tmp)
           end
       end
    end

    return dout
end

"""
   automatic gaining control, L1 is the number of windows for one trace,
"""
function agc(din::Array{Tv}; L1=10, L2=10) where {Tv<:AbstractFloat}

    # get the dimensions
    dims = size(din)
    nt   = dims[1]
    num_traces = prod(dims[2:end])

    # reshape the input as matrix
    din  = reshape(din, dims[1], num_traces)
    dout = zeros(Tv, nt, num_traces)

    # create smooth filter
    hl        = floor(Int64, nt/L1)
    op_length = 2 * hl + 1
    smoother  = triang(op_length)
    smoother  = convert(Vector{Tv}, smoother)

    # envelope of the data
    env  = envelope(din)
    trip = floor(Int64, op_length/L2)

    # loop over traces
    for i = 1 : num_traces

        # smoothing the envelope
        tmp = env[:,i]
        tmp = conv(tmp, smoother)

        # truncate the result of smoothing to have same length
        env_sm = tmp[hl+1 : nt+hl]
        env_sm[1:trip]        .= env_sm[trip+1]
        env_sm[nt-trip+1: end].= env_sm[nt-trip]

        trace_in  = din[:, i]
        trace_out = trace_in ./ env_sm
        dout[:,i].= trace_out / maximum(abs.(trace_out))
    end

    return reshape(dout, dims...)
end

function balance(din::Array{Tv}) where {Tv<:AbstractFloat}

    dims = size(din)
    nt   = dims[1]
    num_traces = prod(dims[2:end])
    din  = reshape(din, nt, num_traces)

    dout  = zeros(Tv, nt, num_traces)
    trref = ones(nt)
    top   = norm(trref)

    for i = 1 : num_traces
        bottom     = norm(din[:,i])
        dout[:,i] .= din[:,i] * top / bottom
    end

    return reshape(dout, dims...)
end
