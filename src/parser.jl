function parser_main(ICP_data::AbstractDataFrame,
                     timestamps::AbstractDataFrame)
    lag, sweep = parser_lag_sweep(ICP_data,timestamps)
    LA_start_of_spot_indices = findall(.!ismissing.(timestamps[:,3])) # SubPoint
    LA_times = time_difference.(timestamps[1,1],timestamps[:,1])
    LA_start_of_spot_times = LA_times[LA_start_of_spot_indices]
    ICP_start_of_spot_indices = parser_LAtime2ICPindex(LA_start_of_spot_times,lag,sweep)
    nspot = length(LA_start_of_spot_indices)
    run = Vector{Sample}(undef,nspot)
    for i in eachindex(run)
        LA_start_of_spot_index = LA_start_of_spot_indices[i]
        ICP_start_of_spot_index = ICP_start_of_spot_indices[i]
        if i == 1
            ICP_start_of_blank_index = 1
        else
            LA_start_of_blank_time = LA_times[LA_start_of_spot_index-1]
            ICP_start_of_blank_index = parser_LAtime2ICPindex(LA_start_of_blank_time,
                                                              lag,sweep)
        end
        if i == nspot
            LA_end_of_spot_index = length(LA_times)
        else
            LA_end_of_spot_index = LA_start_of_spot_indices[i+1] - 1
        end
        LA_end_of_spot_time = LA_times[LA_end_of_spot_index]
        ICP_end_of_spot_index = parser_LAtime2ICPindex(LA_end_of_spot_time,
                                                       lag,sweep)
        selected_ICP_data = ICP_data[ICP_start_of_blank_index:ICP_end_of_spot_index,:]
        selected_timestamps = timestamps[LA_start_of_spot_index:LA_end_of_spot_index,:]
        run[i] = parser_df2sample(selected_ICP_data,
                                  selected_timestamps)
    end
    return run
end

function parser_df2sample(selected_ICP_data::AbstractDataFrame,
                          selected_timestamps::AbstractDataFrame;
                          absolute_buffer::AbstractFloat=2.0,
                          relative_buffer::AbstractFloat=0.1)
    sname = selected_timestamps[1,5] # Comment
    datetime = automatic_datetime(selected_timestamps[1,1])
    LA_on_off_indices = parser_on_off_indices(selected_timestamps)
    LA_on_off_times = time_difference.(selected_timestamps[1,1],
                                       selected_timestamps[LA_on_off_indices,1])
    x = selected_timestamps[LA_on_off_indices,6] # X(um)
    y = selected_timestamps[LA_on_off_indices,7] # Y(um)
    dat = selected_ICP_data
    ICP_times = dat[:,1] .= dat[:,1] .- dat[1,1]
    nsweep = length(ICP_times)
    sweep = ICP_times[end]/nsweep
    lag = round(Int,nsweep*(1-LA_on_off_times[end]/ICP_times[end]))
    ICP_on_off_indices = parser_LAtime2ICPindex(LA_on_off_times,
                                               lag,sweep)
    i0 = ICP_on_off_indices[1]
    t0 = ICP_times[i0]
    bwin = parser_bwin(i0,ICP_times)
    swin = parser_swin(ICP_on_off_indices,x,y,ICP_times)
    return Sample(sname,datetime,dat,t0,bwin,swin,"sample")
end

function parser_bwin(i0::Integer,
                     t::AbstractVector;
                     absolute_buffer::AbstractFloat=2.0,
                     relative_buffer::AbstractFloat=0.1)
    if t[i0] > absolute_buffer
        i_buffer = round(Int,i0*absolute_buffer/t[i0])
    else
        i_buffer = round(Int,i0*relative_buffer)
    end
    i0_minus_buffer = i0 - i_buffer
    return [(1,i0_minus_buffer)]
end

function parser_swin(i_on_off::AbstractVector,
                     x::AbstractVector,
                     y::AbstractVector,
                     t::AbstractVector;
                     absolute_buffer::AbstractFloat=2.0,
                     relative_buffer::AbstractFloat=0.1)
    nwin = round(Int,length(x)/2)
    nt = length(t)
    out = Vector{Tuple}(undef,nwin)
    for i in 1:nwin
        j = 2*i - 1
        i_start = i_on_off[j]
        i_stop = i_on_off[j+1]
        window_width_i = i_stop - i_start
        window_width_t = t[i_stop] - t[i_start]
        if window_width_t > 2*absolute_buffer
            i_buffer = round(Int,window_width_i*absolute_buffer/window_width_t)
        else
            i_buffer = round(Int,window_width_i*relative_buffer)
        end
        i_start_plus_buffer = i_start + i_buffer
        i_stop_minus_buffer = i_stop - i_buffer
        x_buffer = (x[j+1] - x[j])*i_buffer/window_width_i
        x_start_plus_buffer = x[j] + x_buffer
        x_stop_minus_buffer = x[j+1] - x_buffer
        y_buffer = (y[j+1] - y[j])*i_buffer/window_width_i
        y_start_plus_buffer = y[j] + y_buffer
        y_stop_minus_buffer = y[j+1] - y_buffer
        out[i] = (i_start_plus_buffer, i_stop_minus_buffer,
                  x_start_plus_buffer, x_stop_minus_buffer,
                  y_start_plus_buffer, y_stop_minus_buffer)
    end
    return out
end

function parser_lag_sweep(ICPdata::AbstractDataFrame,
                          timestamps::AbstractDataFrame)
    ICPtime = ICPdata[:,1]
    ICPsig = ICPdata[:,2:end]
    totsig = vec(sum(Matrix(ICPsig),dims=2))
    sweep = parser_sweep(ICPtime)
    ion, ioff = parser_interpolated_on_off_indices(timestamps,sweep)
    lower, upper = parser_lag_search_range(timestamps,totsig,ICPtime)
    misfit = parser_lag_misfit(lower,upper,ion,ioff,totsig)
    lag = argmax(misfit)
    return lag, sweep
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

function parser_on_off_indices(timestamps)
    return (cumsum(rle(timestamps[:,11])[2]).+1)[1:end-1] # Laser State
end

function parser_interpolated_on_off_indices(timestamps,sweep)
    onoff = parser_on_off_indices(timestamps)
    LAtime = time_difference.(timestamps[1,1],timestamps[onoff,1])
    start = ceil.(Int,LAtime[1:2:end-1]/sweep)
    stop = floor.(Int,LAtime[2:2:end]/sweep)
    on_intervals = map((s, e) -> (s+1):(e), start, stop)
    off_intervals = map((s, e) -> (s+1):(e), stop[1:end-1], start[2:end])
    ion = reduce(vcat, on_intervals)
    ioff = reduce(vcat, off_intervals)
    return ion, ioff
end

function parser_LAtime2ICPindex(LAtime::AbstractVector,
                                lag::Integer,
                                sweep::Number)
    return parser_LAtime2ICPindex.(LAtime,lag,sweep)
end
function parser_LAtime2ICPindex(LAtime::Number,
                                lag::Integer,
                                sweep::Number)
    return lag + round(Int,LAtime/sweep)
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
    onoff = parser_on_off_indices(timestamps)
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
