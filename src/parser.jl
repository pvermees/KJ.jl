function parser_main(ICPdata::AbstractDataFrame,
                     timestamps::AbstractDataFrame)
    lag, sweep = parser_lag_sweep(ICPdata,timestamps)
    LA_start_indices = findall(.!ismissing.(timestamps[:,3])) # SubPoint
    LA_start_times = time_difference(timestamps[1,1],
                                     timestamps[LA_start_indices,1])
    ICP_start_indices = parser_LAtime2ICPindex(LA_start_times,lag,sweep)
    run = Vector{Sample}(undef,length(LA_start_indices))
    for i in eachindex(run)
        run[i] = parser_i2samp(i,ICPdata,timestamps,
                               LA_start_indices,ICP_start_indices)
    end
    return run
end

function parser_i2samp(i::Integer,
                       ICPdata::AbstractDataFrame,
                       timestamps::AbstractDataFrame,
                       LA_start_indices::AbstractVector,
                       ICP_start_indices::AbstractVector)
    ns = length(LA_start_indices)
    if i == 1
        ICP_start_index = 1
        ICP_stop_index = ICP_start_indices[i+1]-1
    elseif i == ns
        ICP_start_index = ICP_start_indices[i]
        ICP_stop_index = size(ICPdata,1)
    else
        ICP_start_index = ICP_start_indices[i]
        ICP_stop_index = ICP_start_indices[i+1]-1
    end
    LA_start_index = LA_start_indices[i]
    if i == ns
        LA_stop_index = size(timestamps,1)
    else
        LA_stop_index = LA_start_indices[i+1]-1
    end
    sample = parser_df2sample(ICPdata[ICP_start_index:ICP_stop_index,:],
                              timestamps[LA_start_index:LA_stop_index,:])
    return sample
end

function parser_df2sample(selected_ICP_data::AbstractDataFrame,
                          selected_timestamps::AbstractDataFrame)
    sname = selected_timestamps[1,5] # Comment
    datetime = automatic_datetime(selected_timestamps[1,1])
    vertix_number = selected_timestamps[:,4]
    i_vertices = findall(.!ismissing.(vertix_number)) # Vertix Number
    if isempty(i_vertices)
        samp = io_df2sample(selected_ICP_data,sname,datetime)
    else
        swin, bwin = parser_read_vertices(i_vertices,
                                          selected_ICP_data,
                                          selected_timestamps)
    end
    return samp
end

function parser_read_vertices(i_vertices::AbstractVector,
                              selected_ICP_data::AbstractDataFrame,
                              selected_timestamps::AbstractDataFrame)
end

function parser_lag_sweep(ICPdata::AbstractDataFrame,
                          timestamps::AbstractDataFrame)
    ICPtime = ICPdata[:,1]
    ICPsig = ICPdata[:,2:end]
    totsig = vec(sum(Matrix(ICPsig),dims=2))
    sweep = parser_sweep(ICPtime)
    ion, ioff = parser_interpolated_on_off_indices(timestamps,totsig,sweep)
    lower, upper = parser_lag_search_range(timestamps,totsig,ICPtime)
    misfit = parser_lag_misfit(lower,upper,ion,ioff,totsig)
    return argmax(misfit), sweep
end

function parser_lag_misfit(lower::Integer,
                           upper::Integer,
                           ion::AbstractVector,
                           ioff::AbstractVector,
                           totsig::AbstractVector)
    misfit = fill(0.0,upper-lower+1)
    for lag in lower:upper
        signal = totsig[lag .+ ion]
        blank = totsig[lag .+ ioff]
        misfit[lag+1] = log(sum(signal)) - log(sum(blank))
    end
    return misfit
end

function parser_on_off_indices(timestamps,totsig,sweep)
    return (cumsum(rle(timestamps[:,11])[2]).+1)[1:end-1] # Laser State
end

function parser_interpolated_on_off_indices(timestamps,totsig,sweep)
    onoff = parser_on_off_indices(timestamps,totsig,sweep)
    LAtime = time_difference(timestamps[1,1],timestamps[onoff,1])
    start = round.(Int,LAtime[1:2:end-1]/sweep)
    stop = round.(Int,LAtime[2:2:end]/sweep)
    on_intervals = map((s, e) -> (s+1):(e), start, stop)
    off_intervals = map((s, e) -> (s+1):(e), stop[1:end-1], start[2:end])
    ion = reduce(vcat, on_intervals)
    ioff = reduce(vcat, off_intervals)
    return ion, ioff
end

function parser_LAtime2ICPindex(LAtime::AbstractVector,
                                lag::Integer,
                                sweep::Number)
    return @. lag + round(Int,LAtime/sweep)
end

function parser_ICPduration(ICPtime::AbstractVector)
    return ICPtime[end] - ICPtime[1]
end

function parser_sweep(ICPtime::AbstractVector)
    return parser_ICPduration(ICPtime) / (length(ICPtime)-1)
end

function parser_lag_search_range(timestamps::AbstractDataFrame,
                                 totsig::AbstractVector,
                                 ICPtime::AbstractVector)
    sweep = parser_sweep(ICPtime)
    onoff = parser_on_off_indices(timestamps,totsig,sweep)
    LAduration = time_difference(timestamps[1,1],
                                 timestamps[end,1])
    ICPduration = parser_ICPduration(ICPtime)
    if LAduration > ICPduration
        @warn The laser session is longer than the ICP-MS session!
        upper = floor(Int,ICPduration/sweep)
    else
        upper = floor(Int,(ICPduration-LAduration)/sweep)
    end
    return 0, upper
end

function parser_df2sample(df::AbstractDataFrame,
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
