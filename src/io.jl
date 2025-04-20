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
export readFile

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

function readDat(fname::AbstractString,
                 instrument::AbstractString="Agilent",
                 head2name::Bool=true)
    if instrument=="Agilent"
        sname, datetime, header, skipto, footerskip =
            readAgilent(fname,head2name)
    elseif instrument=="FIN2"
        sname, datetime, header, skipto, footerskip =
            readFIN(fname,head2name)
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

function readFIN(fname::AbstractString,
                 head2name::Bool=true)
    lines = split(readuntil(fname, "Time"), "\r\n")
    snamestring = head2name ? lines[3] : fname
    sname = split(split(snamestring,r"[\\/]")[end],".FIN")[1]
    datetime = Dates.DateTime(lines[2],"EEEE, U dd, yyyy HH:MM:SS")
    header = 8
    skipto = 9
    footerskip = 0
    return sname, datetime, header, skipto, footerskip
    
end

function parser_getLowBlankSignal(signal::AbstractDataFrame)
    nc = size(signal,2)
    firstblank = collect(signal[1,:])
    blankrank = sortperm(firstblank)
    lq = firstblank[blankrank[floor(Int,nc/4)]]
    uq = firstblank[blankrank[ceil(Int,3*nc/4)]]
    iqr = uq - lq
    outliers = (firstblank .< lq - 1.5*iqr) .|| (firstblank .> uq + 1.5*iqr)
    return signal[:,.!outliers]
end
function parser_getLaserLag(lowblanksignal::AbstractDataFrame,
                            ICPtime::AbstractVector,
                            timestamps::AbstractDataFrame)
    ICPduration = ICPtime[end]
    onoff = (cumsum(rle(timestamps[:,11])[2]).+1)[1:end-1] # Laser State
    total = sum.(eachrow(lowblanksignal))
    scaled = total./Statistics.mean(total)
    cs = cumsum(scaled)
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
    return lag_to_first_shot - wait_until_first_shot
end
function parser_ICPcutter(lag::AbstractFloat,
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
function parseData(data::AbstractDataFrame,
                   timestamps::AbstractDataFrame)
    ICPtime = data[:,1] # "Time [Sec]"
    lowblanksignal = parser_getLowBlankSignal(data[:,2:end])
    lag = parser_getLaserLag(lowblanksignal,ICPtime,timestamps)
    df = parser_ICPcutter(lag,ICPtime,timestamps)
    run = Vector{Sample}(undef,size(df,1))
    for i in eachindex(run)
        run[i] = df2sample(data,df.snames[i],df.date_times[i],
                           df.ICPstart[i],df.ICPstop[i],df.ICPon[i],df.ICPoff[i])
    end
    return run
end

function dwelltime(run::Vector{Sample})
    channels = getChannels(run)
    out = Dict(channel => 1.0 for channel in channels)
    all = pool(run)
    for channel in channels
        values = sort(unique(all[:,channel]))
        numval = minimum((5,length(values)))
        if numval>1
            dcps = minimum(values[2:numval] .- values[1:numval-1])
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
                         method::AbstractString,
                         channels::AbstractDict,
                         blank::AbstractDataFrame,
                         pars::NamedTuple;
                         PAcutoff::Union{AbstractFloat,Nothing}=nothing,
                         prefix=nothing,
                         fname::AbstractString="KJ.json")
    ratios = averat(run,channels,blank,pars)
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

    snames = ratios[:,1]
    PD = replace(ratios[:,2], NaN => "\"NA\"")
    sPD = replace(ratios[:,3], NaN => "\"NA\"")
    dD = replace(ratios[:,4], NaN => "\"NA\"")
    sdD = replace(ratios[:,5], NaN => "\"NA\"")
    rho = replace(ratios[:,6], NaN => "\"NA\"")

    datastring = "\"ierr\":1,\"data\":{"*
    "\""* P *"/"* D *"\":["*     join(PD,",")*"],"*
    "\"err["* P *"/"* D *"]\":["*join(sPD,",")*"],"*
    "\""* d *"/"* D *"\":["*     join(dD,",")*"],"*
    "\"err["* d *"/"* D *"]\":["*join(sdD,",")*"],"*
    "\"(rho)\":["*join(rho,",")*"],"*
    "\"(C)\":[],\"(omit)\":[],"*
    "\"(comment)\":[\""*join(snames,"\",\"")*"\"]"

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
