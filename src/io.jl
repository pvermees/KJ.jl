"""
load(dname::AbstractString;
     format::AbstractString="Agilent",
     head2name::Bool=true)

Load MS data files
"""
function load(dname::AbstractString;
              format::AbstractString="Agilent",
              head2name::Bool=true)
    fnames = readdir(dname)
    samples = Vector{Sample}(undef,0)
    datetimes = Vector{DateTime}(undef,0)
    ext = getExt(format)
    for fname in fnames
        if occursin(ext,fname)
            try
                pname = joinpath(dname,fname)
                samp = readFile(pname;
                                format=format,
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
        samp.dat.outlier = falses(size(samp.dat,1))
        samp.dat.t = (samp.dat[:,1] .+ runtime[i])./duration
    end
    return sortedsamples
end
"""
load(dfile::AbstractString,
     tfile::AbstractString;
     format::AbstractString="Agilent")

Load and parse an ICP-MS file using time stamps
"""
function load(dfile::AbstractString,
              tfile::AbstractString;
              format::AbstractString="Agilent")
    dat = timestamps = DataFrame()
    try
        dat = readDat(dfile,format,false)[1]
        dat.outlier = falses(size(dat,1))
        dat.t = dat[:,1] ./ dat[end,1]
    catch e
        println("Failed to read "*dfile)
    end
    try
        timestamps = CSV.read(tfile, DataFrame)
    catch e
        println("Failed to read "*tfile)
    end
    return parser_main(dat,timestamps)
end
export load

function readFile(fname::AbstractString;
                  format::AbstractString="Agilent",
                  head2name::Bool=true)
    dat, sname, datetime = readDat(fname,format,head2name)
    return io_df2sample(dat,sname,datetime)
end

function readDat(fname::AbstractString,
                 format::AbstractString="Agilent",
                 head2name::Bool=true)
    if format=="Agilent"
        sname, datetime, header, skipto, footerskip =
            readAgilent(fname,head2name)
    elseif format=="FIN2"
        sname, datetime, header, skipto, footerskip =
            readFIN(fname,head2name)
    elseif format=="ThermoFisher"
        sname, datetime, header, skipto, footerskip =
            readThermoFisher(fname,head2name)
    else
        KJerror("unknownFormat")
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

function io_df2sample(df::AbstractDataFrame,
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


function export2IsoplotR(run::Vector{Sample},
                         method::Gmethod,
                         fit::Gfit;
                         prefix=nothing,
                         fname::AbstractString="KJ.json")
    ratios = averat(run,method,fit)
    if isnothing(prefix)
        export2IsoplotR(ratios,method;fname=fname)
    else
        export2IsoplotR(prefix2subset(ratios,prefix),
                        method;
                        fname=fname)
    end
end

function export2IsoplotR(ratios::AbstractDataFrame,
                         method::Gmethod;
                         fname::AbstractString="KJ.json")
    json = jsonTemplate()

    snames = ratios[:,1]
    PD = replace(ratios[:,2], NaN => "\"NA\"")
    sPD = replace(ratios[:,3], NaN => "\"NA\"")
    dD = replace(ratios[:,4], NaN => "\"NA\"")
    sdD = replace(ratios[:,5], NaN => "\"NA\"")
    rho = replace(ratios[:,6], NaN => "\"NA\"")

    datastring = "\"ierr\":1,\"data\":{"*
    "\""* method.P.ion *"/"* method.D.ion *"\":["*     join(PD,",")*"],"*
    "\"err["* method.P.ion *"/"* method.D.ion *"]\":["*join(sPD,",")*"],"*
    "\""* method.d.ion *"/"* method.D.ion *"\":["*     join(dD,",")*"],"*
    "\"err["* method.d.ion *"/"* method.D.ion *"]\":["*join(sdD,",")*"],"*
    "\"(rho)\":["*join(rho,",")*"],"*
    "\"(C)\":[],\"(omit)\":[],"*
    "\"(comment)\":[\""*join(snames,"\",\"")*"\"]"

    chronometer = method.name

    json = replace(json,"\""*chronometer*"\":{}" =>
                   "\""*chronometer*"\":{"*datastring*"}}")

    
    if chronometer in ["Lu-Hf","Rb-Sr","K-Ca","Re-Os"]
                        
        old = "\"geochronometer\":\"U-Pb\",\"plotdevice\":\"concordia\""
        new = "\"geochronometer\":\""*chronometer*"\",\"plotdevice\":\"isochron\""
        json = replace(json, old => new)
        
        old = "\""*chronometer*"\":{\"format\":1,\"i2i\":true,\"projerr\":false,\"inverse\":false}"
        new = "\""*chronometer*"\":{\"format\":2,\"i2i\":true,\"projerr\":false,\"inverse\":true}"
        json = replace(json, old => new)
        
    end
    
    file = open(fname,"w")
    write(file,json)
    close(file)
    
end
export export2IsoplotR
