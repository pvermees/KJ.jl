"""
load

Read mass spectrometer data

# Returns

- a vector of samples

# Methods

- `load(dname::AbstractString;
        instrument::AbstractString="Agilent",
        head2name::Bool=true)`
- `load(dfile::AbstractString,
        tfile::AbstractString;
        instrument::AbstractString="Agilent")`

# Arguments

- `dname`: directory containing mass spectrometer data files
- `instrument`: one of "Agilent" or "ThermoFisher"
- `head2name`: `true` if sample names should be read from the file headers.
               `false` if they should be extracted from the file names
- `dfile`: single data file
- `tfile`: laser timestamp file

# Examples
```julia
myrun = load("data/Lu-Hf";instrument="Agilent")
p = plot(myrun[1],["Hf176 -> 258","Hf178 -> 260"])
display(p)
```
"""
function load(dname::AbstractString;
              instrument::AbstractString="Agilent",
              head2name::Bool=true)
    fnames = readdir(dname)
    samples = Vector{Sample}(undef,0)
    datetimes = Vector{DateTime}(undef,0)
    ext = getExt(instrument)
    for fname in fnames
        if occursin(ext,fname)
            try
                pname = joinpath(dname,fname)
                samp = readFile(pname;
                                instrument=instrument,
                                head2name=head2name)
                push!(samples,samp)
                push!(datetimes,samp.datetime)
            catch e
                println("Failed to read "*fname)
            end
        end
    end
    order = sortperm(datetimes)
    sortedsamples = samples[order]
    sorteddatetimes = datetimes[order]
    dt = sorteddatetimes .- sorteddatetimes[1]
    runtime = Dates.value.(dt)
    duration = runtime[end] + sortedsamples[end].dat[end,1]
    for i in eachindex(sortedsamples)
        samp = sortedsamples[i]
        samp.dat.t = (samp.dat[:,1] .+ runtime[i])./duration
    end
    return sortedsamples
end
function load(dfile::AbstractString,
              tfile::AbstractString;
              instrument::AbstractString="Agilent")
    samples = Vector{Sample}(undef,0)
    datetimes = Vector{DateTime}(undef,0)
    dat = timestamps = DataFrame()
    try
        dat = readDat(dfile,instrument,false)[1]
        dat.t = dat[:,1] ./ dat[end,1]
    catch e
        println("Failed to read "*dfile)
    end
    try
        timestamps = CSV.read(tfile, DataFrame)
    catch e
        println("Failed to read "*tfile)
    end
    return parseData(dat,timestamps)
end
export load

function readFile(fname::AbstractString;
                  instrument::AbstractString="Agilent",
                  head2name::Bool=true)
    dat, sname, datetime = readDat(fname,instrument,head2name)
    return df2sample(dat,sname,datetime)
end

function df2sample(df::AbstractDataFrame,
                   sname::AbstractString,
                   datetime::DateTime)
    t = df[:,1]
    i0 = geti0(df[:,2:end])
    t0 = df[i0,1]
    bwin = autoWindow(t,t0;blank=true)
    swin = autoWindow(t,t0;blank=false)
    return Sample(sname,datetime,df,t0,bwin,swin,"sample")
end
function df2sample(df::AbstractDataFrame,
                   sname::AbstractString,
                   datetime::DateTime,
                   start::AbstractFloat,
                   stop::AbstractFloat,
                   on::AbstractFloat,
                   off::AbstractFloat)
    t = df[:,1]
    selection = (t.>=start .&& t.<=stop)
    absolute_buffer = 2.0
    relative_buffer = 0.1
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
    return Sample(sname,datetime,df[selection,:],on,bwin,swin,"sample")
end

function readDat(fname::AbstractString,
                 instrument::AbstractString="Agilent",
                 head2name::Bool=true)
    if instrument=="Agilent"
        sname, datetime, header, skipto, footerskip =
            readAgilent(fname,head2name)
    elseif instrument=="ThermoFisher"
        sname, datetime, header, skipto, footerskip =
            readThermoFisher(fname,head2name)
    else
        KJerror("unknownInstrument")
    end
    dat = CSV.read(
        fname,
        DataFrame;
        header = header,
        skipto = skipto,
        footerskip = footerskip,
        ignoreemptyrows = true,
        delim = ',',
    )
    select!(dat, [k for (k,v) in pairs(eachcol(dat)) if !all(ismissing, v)])
    return dat, sname, datetime
end

function readAgilent(fname::AbstractString,
                     head2name::Bool=true)

    lines = split(readuntil(fname, "Time [Sec]"), "\n")
    snamestring = head2name ? lines[1] : fname
    sname = split(split(snamestring,r"[\\/]")[end],".")[1]
    datetimeline = lines[3]
    from = findfirst(":",datetimeline)[1]+2
    to = findfirst("using",datetimeline)[1]-2
    datetime = automatic_datetime(datetimeline[from:to])
    header = 4
    skipto = 5
    footerskip = 3
    
    return sname, datetime, header, skipto, footerskip
    
end

function readThermoFisher(fname::AbstractString,
                          head2name::Bool=true)

    lines = split(readuntil(fname, "Time"), "\n")
    snamestring = head2name ? split(lines[1],":")[1] : fname
    sname = split(split(snamestring,r"[\\/]")[end],".")[1]
    datetimeline = lines[1]
    from = findfirst(":",datetimeline)[1]+1
    to = findfirst(";",datetimeline)[1]-1
    datetime = automatic_datetime(datetimeline[from:to])
    header = 14
    skipto = 16
    footerskip = 0
    
    return sname, datetime, header, skipto, footerskip
    
end

function parseData(data::AbstractDataFrame,
                   timestamps::AbstractDataFrame)
    run = Vector{Sample}(undef,0)
    nr,nc = size(data)
    ICPtime = data[:,1] # "Time [Sec]"
    signal = data[:,2:end]
    onoff = (cumsum(rle(timestamps[:,11])[2]).+1)[1:end-1] # Laser State
    # 1. select the low blank signals
    firstblank = collect(signal[1,:])
    blankrank = sortperm(firstblank)
    lq = firstblank[blankrank[floor(Int,nc/4)]]
    uq = firstblank[blankrank[ceil(Int,3*nc/4)]]
    iqr = uq - lq
    outliers = (firstblank .< lq - 1.5*iqr) .|| (firstblank .> uq + 1.5*iqr)
    lowblanksignal = signal[:,.!outliers]
    # 2. get the cumulative signal of the low blank signals
    total = sum.(eachrow(lowblanksignal))
    scaled = total./Statistics.mean(total)
    cs = cumsum(scaled)
    # 3. find the lag time between the ICP-MS and laser files
    ICPduration = ICPtime[end]
    LAduration = time_difference(timestamps[onoff[1],1],
                                 timestamps[onoff[end],1])
    lower = 0.0
    if LAduration > ICPduration
        @warn The laser session is longer than the ICP-MS session!
        upper = ICPduration
    else
        upper = ICPduration - LAduration
    end
    coverage = function(lag)
        i1 = argmin(abs.(ICPtime .- lag))
        i2 = argmin(abs.(ICPtime .< lag + LAduration))
        return cs[i2]-cs[i1]
    end
    t = ICPtime[ICPtime.>lower .&& ICPtime.<upper]
    lag_to_first_shot = t[argmax(coverage.(t))]
    wait_until_first_shot = time_difference(timestamps[1,1],
                                            timestamps[onoff[1],1])
    lag = lag_to_first_shot - wait_until_first_shot
    # 4. create vectors with the start and end of each sequence, and laser timings
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
    date_times = automatic_datetime.(timestamps[sequences,1])
    snames = timestamps[sequences,5] # "Comment"
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
    # 5. parse the data into samples
    for i in 1:ns
        samp = df2sample(data,snames[i],date_times[i],
                         ICPstart[i],ICPstop[i],ICPon[i],ICPoff[i])
        push!(run,samp)
    end
    return run
end

function dwelltime(run::Vector{Sample})
    channels = getChannels(run)
    out = Dict(channel => 1.0 for channel in channels)
    blk = pool(run;blank=true)
    for channel in channels
        values = sort(unique(blk[:,channel]))
        if length(values)>1
            dcps = minimum(values[2:end] .- values[1:end-1])
            out[channel] = minimum([1/dcps,out[channel]])
        end
    end
    return out
end
export dwelltime

"""
export2IsoplotR

Export isotopic ratio data to an IsoplotRgui json file

# Methods

- `export2IsoplotR(run::Vector{Sample},
                   dt::AbstractDict,
                   method::AbstractString,
                   channels::AbstractDict,
                   blank::AbstractDataFrame,
                   pars::NamedTuple;
                   PAcutoff=nothing,prefix=nothing,
                   fname::AbstractString="KJ.json")`
- `export2IsoplotR(ratios::AbstractDataFrame,
                   method::AbstractString;
                   fname::AbstractString="KJ.json")`

# Arguments

- `run`: the output of `load`
- `dt`: the dwell times (in seconds) of the mass channels
- `method`: a geochronometer (e.g., `Lu-Hf`, `Rb-Sr`, `U-Pb`)
- `channels`: dictionary of the type Dict("P" => "parent", "D" => "daughter", "d" => "sister")
- `blank`: the output of fitBlanks()
- `pars`: the output of fractionation() or process!()
- `PAcutoff`: the pulse-analog signal cutoff
- `fname`: path of the output file

# Examples
```julia
myrun = load("data/Lu-Hf",instrument="Agilent")
dt = dwelltime(myrun)
method = "Lu-Hf"
channels = Dict("d"=>"Hf178 -> 260",
                "D"=>"Hf176 -> 258",
                "P"=>"Lu175 -> 175")
standards = Dict("Hogsbo" => "hogsbo")
glass = Dict("NIST612" => "NIST612p")
cutoff = 1e7
blk, fit = process!(myrun,dt,method,channels,standards,glass;
                    PAcutoff=cutoff,nblank=2,ndrift=1,ndown=1)
selection = prefix2subset(ratios,"BP")
export2IsoplotR(selection,"Lu-Hf",fname="BP.json")
```
"""
function export2IsoplotR(run::Vector{Sample},
                         dt::AbstractDict,
                         method::AbstractString,
                         channels::AbstractDict,
                         blank::AbstractDataFrame,
                         pars::NamedTuple;
                         PAcutoff=nothing,
                         prefix=nothing,
                         fname::AbstractString="KJ.json")
    ratios = averat(run,dt,channels,blank,pars)
    if isnothing(prefix)
        export2IsoplotR(ratios,method;fname=fname)
    else
        export2IsoplotR(prefix2subset(ratios,prefix),method;fname=fname)
    end
end
function export2IsoplotR(ratios::AbstractDataFrame,
                         method::AbstractString;
                         fname::AbstractString="KJ.json")
    json = jsonTemplate()

    P, D, d = getPDd(method)

    datastring = "\"ierr\":1,\"data\":{"*
    "\""* P *"/"* D *"\":["*     join(ratios[:,2],",")*"],"*
    "\"err["* P *"/"* D *"]\":["*join(ratios[:,3],",")*"],"*
    "\""* d *"/"* D *"\":["*     join(ratios[:,4],",")*"],"*
    "\"err["* d *"/"* D *"]\":["*join(ratios[:,5],",")*"],"*
    "\"(rho)\":["*join(ratios[:,6],",")*"],"*
    "\"(C)\":[],\"(omit)\":[],"*
    "\"(comment)\":[\""*join(ratios[:,1],"\",\"")*"\"]"

    json = replace(json,"\""*method*"\":{}" =>
                   "\""*method*"\":{"*datastring*"}}")

    
    if method in ["Lu-Hf","Rb-Sr","K-Ca"]
                        
        old = "\"geochronometer\":\"U-Pb\",\"plotdevice\":\"concordia\""
        new = "\"geochronometer\":\""*method*"\",\"plotdevice\":\"isochron\""
        json = replace(json, old => new)
        
        old = "\""*method*"\":{\"format\":1,\"i2i\":true,\"projerr\":false,\"inverse\":false}"
        new = "\""*method*"\":{\"format\":2,\"i2i\":true,\"projerr\":false,\"inverse\":true}"
        json = replace(json, old => new)
        
    end
    
    file = open(fname,"w")
    write(file,json)
    close(file)
    
end
export export2IsoplotR
