"""
    setBwin!(run::Vector{Sample}, bwin::AbstractVector; seconds=false)
    setBwin!(samp::Sample, bwin::AbstractVector; seconds=false)
    setBwin!(run::Vector{Sample})
    setBwin!(samp::Sample)

Set the blank window (bwin) for one or more samples.

If `bwin` is provided, it is directly set (as indices or time in seconds if `seconds=true`).
If no `bwin` is provided, an automatic blank window is calculated.

# Arguments
- `run`/`samp`: Vector of samples or single sample
- `bwin`: Window specification (optional)
- `seconds`: If true, interpret window as time in seconds; otherwise as indices
"""
function setBwin!(run::Vector{Sample},
                  bwin::AbstractVector;
                  seconds::Bool=false)
    for i in eachindex(run)
        setBwin!(run[i],bwin;seconds=seconds)
    end
end

function setBwin!(samp::Sample,
                  bwin::AbstractVector;
                  seconds::Bool=false)
    samp.bwin = seconds ? time2window(samp,bwin) : bwin
end

function setBwin!(run::Vector{Sample})
    for i in eachindex(run)
        setBwin!(run[i])
    end
end

function setBwin!(samp::Sample)
    bwin = autoWindow(samp,blank=true)
    setBwin!(samp,bwin)
end
export setBwin!

"""
    setSwin!(run::Vector{Sample}, swin::AbstractVector; seconds=false)
    setSwin!(samp::Sample, swin::AbstractVector; seconds=false)
    setSwin!(run::Vector{Sample})
    setSwin!(samp::Sample)

Set the signal window (swin) for one or more samples.

If `swin` is provided, it is directly set (as indices or time in seconds if `seconds=true`).
If no `swin` is provided, an automatic signal window is calculated.

# Arguments
- `run`/`samp`: Vector of samples or single sample
- `swin`: Window specification (optional)
- `seconds`: If true, interpret window as time in seconds; otherwise as indices
"""
function setSwin!(run::Vector{Sample},
                  swin::AbstractVector;
                  seconds::Bool=false)
    for i in eachindex(run)
        setSwin!(run[i],swin;seconds=seconds)
    end
end

function setSwin!(samp::Sample,
                  swin::AbstractVector;
                  seconds::Bool=false)
    samp.swin = seconds ? time2window(samp,swin) : swin
end

function setSwin!(run::Vector{Sample})
    for samp in run
        setSwin!(samp)
    end
end

function setSwin!(samp::Sample)
    swin = autoWindow(samp,blank=false)
    setSwin!(samp,swin)
end
export setSwin!

"""
    shift_windows!(run::Vector{Sample}, shift::Number=0.0)

Shift both blank and signal windows for all samples in a run by a given time offset.

# Arguments
- `run`: Vector of samples
- `shift`: Time shift in seconds (default: 0.0)
"""
function shift_windows!(run::Vector{Sample},
                        shift::Number=0.0)
    for samp in run
        samp.t0 += shift
        di = t2i(samp,shift)
        nt = size(samp.dat,1)
        bwin = []
        for win in samp.bwin
            start = max(win[1]+di,1)
            stop = min(win[2]+di,nt)
            push!(bwin,(start,stop))
        end
        setBwin!(samp,bwin)
        swin = []
        for win in samp.swin
            start = max(win[1]+di,1)
            stop = min(win[2]+di,nt)
            push!(swin,(start,stop))
        end
        setSwin!(samp,swin)
    end
end
export shift_windows!

"""
    bwinData(samp::Sample; add_xy=false)

Extract data from the blank window of a sample.

# Arguments
- `samp`: Sample to extract data from
- `add_xy`: If true, include x,y coordinates if available

# Returns
- DataFrame containing the windowed data
"""
function bwinData(samp::Sample;
                  add_xy::Bool=false)
    windows = samp.bwin
    selection, x, y = windows2selection(windows;add_xy=add_xy)
    selected_dat =  samp.dat[selection,:]
    return selected_dat
end
export bwinData

"""
    swinData(samp::Sample; add_xy=false)

Extract data from the signal window of a sample.

The returned data includes a `T` column with relative time in minutes from t0.
If x,y coordinates are available and `add_xy=true`, they are also included.

# Arguments
- `samp`: Sample to extract data from
- `add_xy`: If true, include x,y coordinates if available

# Returns
- DataFrame containing the windowed data with time column `T` in minutes
"""
function swinData(samp::Sample;
                  add_xy::Bool=false)
    windows = samp.swin
    selection, x, y = windows2selection(windows;add_xy=add_xy)
    selected_dat =  samp.dat[selection,:]
    selected_dat.T = (selected_dat[:,1] .- samp.t0)./60 # in minutes
    if !(isnothing(x) || isnothing(y))
        selected_dat.x = x
        selected_dat.y = y
    end
    return selected_dat
end
export swinData

function windows2selection(windows::AbstractVector;
                           add_xy::Bool=false)
    selection = Integer[]
    add_xy = add_xy && length(windows[1])>2
    if add_xy
        x = Float64[]
        y = Float64[]
    else
        x = y = nothing
    end
    for w in windows
        append!(selection, w[1]:w[2])
        if add_xy
            nsweeps = w[2]-w[1]+1
            append!(x,range(w[3],w[4];length=nsweeps))
            append!(y,range(w[5],w[6];length=nsweeps))
        end
    end
    return selection, x, y
end

function string2windows(run::Vector{Sample},text::AbstractString,single::Bool)
    return string2windows(run[1],text,single)
end
function string2windows(samp::Sample,text::AbstractString,single::Bool)
    if single
        parts = split(text,',')
        stime = [parse(Float64,parts[1])]
        ftime = [parse(Float64,parts[2])]
        nw = 1
    else
        parts = split(text,['(',')',','])
        stime = parse.(Float64,parts[2:4:end])
        ftime = parse.(Float64,parts[3:4:end])
        nw = Int(round(size(parts,1)/4))
    end
    windows = Vector{Tuple}(undef,nw)
    t = samp.dat[:,1]
    nt = size(t,1)
    maxt = t[end]
    for i in 1:nw
        if stime[i]>t[end]
            stime[i] = t[end-1]
            print("Warning: start point out of bounds and truncated to ")
            print(string(stime[i]) * " seconds.")
        end
        if ftime[i]>t[end]
            ftime[i] = t[end]
            print("Warning: end point out of bounds and truncated to ")
            print(string(maxt) * " seconds.")
        end
        windows[i] = time2window(samp,stime[i],ftime[i])
    end
    return windows
end
