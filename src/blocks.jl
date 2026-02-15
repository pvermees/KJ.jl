function blocks(run::Vector{Sample},
                blocksize::Integer;
                absolute_buffer::AbstractFloat=2.0,
                relative_buffer::AbstractFloat=0.1)
    ns = Int(round(length(run)/blocksize))
    out = Vector{Sample}(undef,ns)
    for i in eachindex(out)
        start = (i-1)*blocksize + 1
        stop = start + blocksize - 1
        out[i] = merge_samples(run[start:stop];
                               absolute_buffer=absolute_buffer,
                               relative_buffer=relative_buffer)
    end
    return out
end
export blocks

function merge_samples(run::Vector{Sample};
                       absolute_buffer::AbstractFloat=2.0,
                       relative_buffer::AbstractFloat=0.1)
    nblocks = length(run)
    if nblocks < 2
        return run[1]
    end
    sname = ""
    datetime = run[1].datetime
    dat = copy(run[1].dat)
    t0 = 0.0
    bwin = nothing
    swin = nothing
    group = run[1].group
    maxsig = -Inf
    prev_length = size(dat,1)
    next_length = 0
    for i in 2:nblocks
        next_dat = run[i].dat
        next_dat[:,1] .+= dat[end,1] + next_dat[2,1] - next_dat[1,1]
        next_length = size(next_dat,1)
        append!(dat,next_dat)
        signals = getSignals(run[i])
        total_sig = dataframe_sum(signals)
        if total_sig > maxsig
            sname = run[i].sname
            maxsig = total_sig
            t0 = next_dat[1,1]
            end_of_signal = size(dat,1)
            start_of_signal = end_of_signal - next_length
            start_of_blank = start_of_signal - prev_length + 1
            bwin = crop_windows(dat[:,1],
                                dat[start_of_blank,1],
                                dat[start_of_signal,1];
                                absolute_buffer=absolute_buffer,
                                relative_buffer=relative_buffer)
            swin = crop_windows(dat[:,1],
                                dat[start_of_signal,1],
                                dat[end_of_signal,1];
                                absolute_buffer=absolute_buffer,
                                relative_buffer=relative_buffer)
        end
        prev_length = next_length
    end
    return Sample(sname,datetime,dat,t0,bwin,swin,group)
end

function crop_windows(t::AbstractVector,
                      tstart::Number,
                      tstop::Number;
                      absolute_buffer::AbstractFloat=2.0,
                      relative_buffer::AbstractFloat=0.1)
    duration = tstop - tstart
    if duration > 2*absolute_buffer
        window_start = tstart + absolute_buffer
        window_stop = tstop - absolute_buffer
    else
        window_start = tstart + duration*relative_buffer
        window_stop = tstop - duration*relative_buffer
    end
    istart = findfirst(t .> window_start)
    istop = findfirst(t .> window_stop)
    [(istart,istop)]
end
