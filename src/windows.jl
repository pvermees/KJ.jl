"""
setBwin!(run::Vector{Sample},
         bwin::AbstractVector;
         seconds::Bool=false)

Set the blank windows of the entire run.
bwin is a vector of tuples
"""
function setBwin!(run::Vector{Sample},
                  bwin::AbstractVector;
                  seconds::Bool=false)
    for i in eachindex(run)
        setBwin!(run[i],bwin;seconds=seconds)
    end
end

"""
setBwin!(samp::Sample,
         bwin::AbstractVector;
         seconds::Bool=false)
"""
function setBwin!(samp::Sample,
                  bwin::AbstractVector;
                  seconds::Bool=false)
    samp.bwin = seconds ? time2window(samp,bwin) : bwin
end

"""
setBwin!(run::Vector{Sample})

Automatically set the blank windows for an entire run.
"""
function setBwin!(run::Vector{Sample})
    for i in eachindex(run)
        setBwin!(run[i])
    end
end

"""
setBwin!(samp::Sample)

Automatically set the blank window for a single sample.
"""
function setBwin!(samp::Sample)
    bwin = autoWindow(samp,blank=true)
    setBwin!(samp,bwin)
end
export setBwin!

"""
setSwin!(run::Vector{Sample},
         swin::AbstractVector;
         seconds::Bool=false)

Set the signal windows of an entire run.
swin is a vector of tuples
"""
function setSwin!(run::Vector{Sample},
                  swin::AbstractVector;
                  seconds::Bool=false)
    for i in eachindex(run)
        setSwin!(run[i],swin;seconds=seconds)
    end
end

"""
setSwin!(samp::Sample,
         swin::AbstractVector;
         seconds::Bool=false)
"""
function setSwin!(samp::Sample,
                  swin::AbstractVector;
                  seconds::Bool=false)
    samp.swin = seconds ? time2window(samp,swin) : swin
end

"""
setSwin!(run::Vector{Sample})

Automatically set the signal windows for an entire run.
"""
function setSwin!(run::Vector{Sample})
    for samp in run
        setSwin!(samp)
    end
end

"""
setSwin!(samp::Sample)

Automatically set the signal window for a single sample.
"""
function setSwin!(samp::Sample)
    swin = autoWindow(samp,blank=false)
    setSwin!(samp,swin)
end
export setSwin!

"""
shift_windows!(run::Vector{Sample},
               shift::Number=0.0)

Shift the blank and signal windows to the left or the right
by a specified number of integrations.
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
bwinData(samp::Sample;
         add_xy::Bool=false)

Return the blank data from samp.
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
swinData(samp::Sample;
         add_xy::Bool=false)

Return the signal data from samp.
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

function string2windows(samp::Sample,text::AbstractString,single::Bool)
    @infiltrate
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
