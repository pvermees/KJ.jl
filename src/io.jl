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

function df2sample(df::AbstractDataFrame,sname::AbstractString,datetime::DateTime)
    t = df[:,1]
    i0 = geti0(df[:,2:end])
    t0 = df[i0,1]
    bwin = autoWindow(t,t0;blank=true)
    swin = autoWindow(t,t0;blank=false)
    return Sample(sname,datetime,df,t0,bwin,swin,"sample")
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
    # 1. get the cumulative signal
    nr = size(data,1)
    runtime = data[:,1] # "Time [Sec]"
    signal = data[:,2:end]
    total = sum.(eachrow(signal))
    scaled = total./Statistics.mean(total)
    cs = cumsum(scaled)
    # 2. find the lag time between the laser and ICP-MS files
    ICPduration = runtime[end]
    lasertime = automatic_datetime.(timestamps[:,1]) # "Timestamp"
    start = lasertime[1]
    stop = lasertime[end]
    LAduration = Millisecond(stop - start).value/1000
    lower = 0.0
    if LAduration>ICPduration
        @warn The laser session is longer than the ICP-MS session!
        upper = ICPduration
    else
        upper = ICPduration - LAduration
    end
    misfit = function(lag)
        i1 = argmin(abs.(runtime .- lag))
        i2 = argmin(abs.(runtime .< lag + LAduration))
        log(cs[end]) - log(cs[i2]-cs[i1])
    end
    crude = argmin(misfit.(lower:1.0:upper))
    fit = Optim.optimize(misfit,runtime[crude-1],runtime[crude+1])
    lag = Optim.minimizer(fit)
    # 3. parse the signals into samples
    sequences = findall(!ismissing,timestamps[:,2]) # "Sequence Number"
    laser_on = findall(==("On"),timestamps[:,11]) # "Laser State"
    i1 = argmin(abs.(runtime .- lag))
    i2 = argmin(abs.(runtime .< lag + LAduration))
    for i in eachindex(sequences)
        datetime = from = lasertime[laser_on[i]]
        to = lasertime[laser_on[i]+2]
        t1 =  Millisecond(from-start).value/1000
        t2 =  Millisecond(to-start).value/1000
        first = maximum([1,floor(Int,(i2-i1)*t1/LAduration)])
        last = minimum([ceil(Int,(i2-i1)*t2/LAduration),nr])
        sname = timestamps[sequences[i],5] # "Comment"
        samp = df2sample(data[first:last,:],sname,datetime)
        push!(run,samp)
    end
    return run
end

"""
export2IsoplotR

Export isotopic ratio data to an IsoplotRgui json file

# Methods

- `export2IsoplotR(run::Vector{Sample},
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
- `method`: a geochronometer (e.g., `Lu-Hf`, `Rb-Sr`, `U-Pb`)
- `channels`: dictionary of the type Dict("P" => "parent", "D" => "daughter", "d" => "sister")
- `blank`: the output of fitBlanks()
- `pars`: the output of fractionation() or process!()
- `PAcutoff`: the pulse-analog signal cutoff
- `fname`: path of the output file

# Examples
```julia
myrun = load("data/Lu-Hf",instrument="Agilent")
method = "Lu-Hf"
channels = Dict("d"=>"Hf178 -> 260",
                "D"=>"Hf176 -> 258",
                "P"=>"Lu175 -> 175")
standards = Dict("Hogsbo" => "hogsbo")
glass = Dict("NIST612" => "NIST612p")
cutoff = 1e7
blk, fit = process!(myrun,method,channels,standards,glass;
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
                         PAcutoff=nothing,
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

    
    if method in ["Lu-Hf","Rb-Sr"]
                        
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
