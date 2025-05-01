function parseData(data::AbstractDataFrame,
                   timestamps::AbstractDataFrame)
    ICPtime = data[:,1] # "Time [Sec]"
    totsig = vec(sum(Matrix(data[:,2:end]),dims=2))
    lag = getLaserLag(totsig,ICPtime,timestamps)
    df = ICPcutter(lag,ICPtime,timestamps)
    run = Vector{Sample}(undef,size(df,1))
    for i in eachindex(run)
        run[i] = df2sample(data,df.snames[i],df.date_times[i],
                           df.ICPstart[i],df.ICPstop[i],df.ICPon[i],df.ICPoff[i])
    end
    return run
end

function getLaserLag(totsig::AbstractVector,
                     ICPtime::AbstractVector,
                     timestamps::AbstractDataFrame)
    # generated indices of laser on and off state
    ICPduration = ICPtime[end] - ICPtime[1]
    sweep = ICPduration / (length(ICPtime)-1)
    onoff = (cumsum(rle(timestamps[:,11])[2]).+1)[1:end-1] # Laser State
    LAtime = time_difference(timestamps[1,1],timestamps[onoff,1])
    start = round.(Int,LAtime[1:2:end-1]/sweep)
    stop = round.(Int,LAtime[2:2:end]/sweep)
    on_intervals = map((s, e) -> (s+1):(e), start, stop)
    off_intervals = map((s, e) -> (s+1):(e), stop[1:end-1], start[2:end])
    ion = reduce(vcat, on_intervals)
    ioff = reduce(vcat, off_intervals)
    # choose the search range
    lower = 0
    LAduration = time_difference(timestamps[onoff[1],1],
                                 timestamps[onoff[end],1])
    if LAduration > ICPduration
        @warn The laser session is longer than the ICP-MS session!
        upper = ceil(Int,ICPduration/sweep)
    else
        upper = ceil(Int,(ICPduration-LAduration)/sweep)
    end
    # maximise the signal to blank (log)ratio
    misfit = fill(0.0,upper-lower+1)
    for lag in lower:upper
        signal = totsig[lag .+ ion]
        blank = totsig[lag .+ ioff]
        misfit[lag+1] = log(sum(signal)) - log(sum(blank))
    end
    return argmax(misfit)
end

function ICPcutter(lag::AbstractFloat,
                   ICPtime::AbstractVector,
                   timestamps::AbstractDataFrame)
    sequences = findall(!ismissing,timestamps[:,2]) # "Sequence Number"
    ns = length(sequences)
    onoff = rle(timestamps[:,11])
    i_onoff = [1;cumsum(onoff[2][1:end-1]) .+ 1]
    i_on = i_onoff[onoff[1].=="On"]
    i_off = i_onoff[onoff[1].=="Off"]
    ICPstart = fill(0.0,ns)
    ICPon = fill(0.0,ns)
    ICPoff = fill(0.0,ns)
    ICPstop = fill(0.0,ns)
    for i in 2:ns-1
        icurrent, inext = sequences[i:i+1]
        i_previous_off = i_off[i_off .< icurrent][end]
        i_current_on = i_on[i_on .> icurrent][1]
        i_current_off = i_off[i_off .< inext][end]
        i_next_on = i_on[i_on .> inext][1]
        t_previous_off = time_difference(timestamps[1,1],timestamps[i_previous_off,1])
        t_current_on = time_difference(timestamps[1,1],timestamps[i_current_on,1])
        t_current_off = time_difference(timestamps[1,1],timestamps[i_current_off,1])
        t_next_on = time_difference(timestamps[1,1],timestamps[i_next_on,1])
        ICPstart[i] = lag + (t_previous_off + t_current_on)/2
        ICPon[i] = lag + t_current_on
        ICPoff[i] = lag + t_current_off
        ICPstop[i] = lag + (t_current_off + t_next_on)/2
    end
    ICPon[1] = lag
    ICPoff[1] = lag + time_difference(timestamps[1,1],timestamps[i_off[2],1])
    ICPstop[1] = ICPstart[2]
    ICPstart[end] = ICPstop[end-1]
    ICPon[end] = lag + time_difference(timestamps[1,1],timestamps[i_on[end],1])
    ICPoff[end] = lag + time_difference(timestamps[1,1],timestamps[i_off[end],1])
    ICPstop[end] = ICPtime[end]
    return DataFrame(snames=timestamps[sequences,5],
                     date_times=automatic_datetime.(timestamps[sequences,1]),
                     ICPstart=ICPstart,ICPon=ICPon,ICPoff=ICPoff,ICPstop=ICPstop)
end

function df2sample(df::AbstractDataFrame,
                   sname::AbstractString,
                   datetime::DateTime;
                   absolute_buffer::AbstractFloat=2.0,
                   relative_buffer::AbstractFloat=0.1)
    t = df[:,1]
    i0 = geti0(df[:,2:end])
    t0 = df[i0,1]
    bwin = autoWindow(t,t0;blank=true,
                      absolute_buffer=absolute_buffer,
                      relative_buffer=relative_buffer)
    swin = autoWindow(t,t0;blank=false,
                      absolute_buffer=absolute_buffer,
                      relative_buffer=relative_buffer)
    return Sample(sname,datetime,df,t0,bwin,swin,"sample")
end
function df2sample(df::AbstractDataFrame,
                   sname::AbstractString,
                   datetime::DateTime,
                   start::AbstractFloat,
                   stop::AbstractFloat,
                   on::AbstractFloat,
                   off::AbstractFloat;
                   absolute_buffer::AbstractFloat=2.0,
                   relative_buffer::AbstractFloat=0.1)
    t = df[:,1]
    selection = (t.>=start .&& t.<=stop)
    bwin = autoBwin(t[selection],on;
                    off=off,start=start,stop=stop,
                    absolute_buffer=absolute_buffer,
                    relative_buffer=relative_buffer)
    swin = autoSwin(t[selection],on;
                    off=off,start=start,stop=stop,
                    absolute_buffer=absolute_buffer,
                    relative_buffer=relative_buffer)
    
    if (on-start) > absolute_buffer
        t2 = on - absolute_buffer
    else
        t2 = on - (on - start)*(1 - relative_buffer)
    end
    i1 = 1
    i2 = findall(t[selection] .< t2)[end]
    bwin = [(i1,i2)]
    if (off-on) > 2*absolute_buffer
        t1 = on + absolute_buffer
        t2 = off - absolute_buffer
    else
        t1 = on + (off - on)*(1 - relative_buffer)
        t2 = off - (off - on)*(1 - relative_buffer)
    end
    i1 = findall(t[selection] .< t1)[end]
    i2 = findall(t[selection] .< t2)[end]
    swin = [(i1,i2)]
    dat = df[selection,:]
    dat[:,1] .= df[selection,1] .- start # change time
    return Sample(sname,datetime,dat,0.0,bwin,swin,"sample")
end
