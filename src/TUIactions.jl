function TUIinit!()
    _KJ["ctrl"] = TUIinit()
    return nothing
end
function TUIinit()
    return Dict(
        "priority" => Dict("load" => true, 
                           "method" => true,
                           "standards" => true, 
                           "glass" => false,
                           "process" => true),
        "history" => DataFrame(task=String[],action=String[]),
        "chain" => ["top"],
        "template" => false,
        "multifile" => true,
        "head2name" => true,
        "format" => "",
        "ICPpath" => "",
        "LApath" => "",
        "run" => nothing,
        "method" => nothing,
        "fit" => nothing,
        "i" => 1,
        "den" => nothing,
        "transformation" => "sqrt",
        "mapcolumn" => 2,
        "clims" => nothing,
        "log" => false,
        "cache" => nothing
    )
end

function TUIwelcome()
    version = string(pkgversion(@__MODULE__))
    title = " KJ "*version*" \n"
    width = Base.length(title)-1
    println('-'^width*"\n"*title*'-'^width)
end

function TUIcheck(ctrl::AbstractDict,action::AbstractString)
    return ctrl["priority"][action] ? "[*]" : ""
end

function TUIread(ctrl::AbstractDict)
    if ctrl["template"]
        if ctrl["multifile"]
            return "loadICPdir"
        else
            return "loadICPfile"
        end
    else
        return "format"
    end
end

function TUIformat!(ctrl::AbstractDict,
                    response::AbstractString)
    if response=="a"
        ctrl["format"] = "Agilent"
    elseif response=="t"
        ctrl["format"] = "ThermoFisher"
    elseif response=="f"
        ctrl["format"] = "FIN2"
    else
        @warn "Unsupported format"
        return "x"
    end
    return "dir|file"
end

function TUIloadICPdir!(ctrl::AbstractDict,
                        response::AbstractString)
    ctrl["run"] = load(response;
                       format=ctrl["format"],
                       head2name=ctrl["head2name"])
    ctrl["priority"]["load"] = false
    ctrl["multifile"] = true
    ctrl["ICPpath"] = response
    if ctrl["template"]
        TUIsetGroups!(ctrl,"standards")
        TUIsetGroups!(ctrl,"glass")
        return "x"
    else
        return "xxx"
    end
end

function TUIloadICPfile!(ctrl::AbstractDict,
                         response::AbstractString)
    ctrl["ICPpath"] = response
    return "loadLAfile"
end

function TUIloadLAfile!(ctrl::AbstractDict,
                        response::AbstractString)
    ctrl["LApath"] = response
    TUIloadICPLAdata!(ctrl::AbstractDict)
end

function TUIloadICPLAdata!(ctrl::AbstractDict)
    ctrl["run"] = load(ctrl["ICPpath"],ctrl["LApath"];
                       format=ctrl["format"])
    ctrl["priority"]["load"] = false
    ctrl["head2name"] = true
    ctrl["multifile"] = false
    if ctrl["template"]
        ctrl["priority"]["standards"] = false
        ctrl["priority"]["glass"] = false
        return "xx"
    else
        return "xxxx"
    end
end

function TUIchoosedir!(ctrl::AbstractDict,
                       response::AbstractString)
    ctrl["ICPpath"] = response
    files = readdir(response)
    for (i, file) in enumerate(files)
        @printf "%2d. %s\n" i file
    end
    return "pickICPLAfiles"
end

function TUIpickICPLAfiles!(ctrl::AbstractDict,
                            response::AbstractString)
    selection = split(response,',')
    ICPnum = parse(Int,selection[1])
    LAnum = parse(Int,selection[2])
    files = readdir(ctrl["ICPpath"];join=true)
    ctrl["ICPpath"] = files[ICPnum]
    ctrl["LApath"] = files[LAnum]
    return TUIloadICPLAdata!(ctrl)
end

function TUImethod!(ctrl::AbstractDict,
                    response::AbstractString)
    if response=="c"
        ctrl["method"] = Cmethod(ctrl["run"],Dict(),(nothing,nothing))
        return "internal"
    else
        i = parse(Int,response)
        methodname = _KJ["methods"].names[i]
        ctrl["method"] = Gmethod(methodname,Dict())
        return "columns"
    end
end

function TUItabulate(ctrl::AbstractDict)
    summarise(ctrl["run"];verbose=true)
    return nothing
end

function TUIinternal!(ctrl::AbstractDict,
                      response::AbstractString)
    i = parse(Int,response)
    ctrl["cache"] = ctrl["channels"][i]
    return "mineral"
end

function TUIstoichiometry!(ctrl::AbstractDict,
                           response::AbstractString)
    channel = ctrl["cache"]
    concentration = parse(Float64,response)
    ctrl["method"].internal = (channel,concentration)
    ctrl["priority"]["method"] = false
    ctrl["priority"]["standards"] = false
    return "xxxx"
end

function TUIchooseMineral!(ctrl::AbstractDict,
                           response::AbstractString)
    if response == "m"
        return "stoichiometry"
    else
        i = parse(Int,response)
        mineral = _KJ["stoichiometry"].names[i]
        channel = ctrl["cache"]
        ctrl["method"].internal = getInternal(mineral,channel)
        ctrl["priority"]["method"] = false
        ctrl["priority"]["standards"] = false
        return "xxx"
    end
end

function TUIsetChannels!(ctrl::AbstractDict,
                         response::AbstractString)
    labels = getChannels(ctrl["run"])
    selected = parse.(Int,split(response,","))
    setChannels!(ctrl["method"];
                 P=labels[selected[1]],
                 D=labels[selected[2]],
                 d=labels[selected[3]])
    equivocal = channels2proxies!(ctrl["method"])
    if equivocal
        return "setProxies"
    else
        ctrl["priority"]["method"] = false
        return "xx"
    end
end

function TUIsetProxies!(ctrl::AbstractDict,
                        response::AbstractString)
    selection = parse.(Int,split(response,","))
    ions = getIons(ctrl["method"])
    nuclides = ions2nuclidelist(ions)
    setProxies(ctrl["method"];
               P=nuclides[selection[1]],
               D=nuclides[selection[2]],
               d=nuclides[selection[3]])
    return "xxx"
end

function TUIchooseStandard!(ctrl::AbstractDict,
                            response::AbstractString)
    i = parse(Int,response)
    ctrl["cache"] = _KJ["refmat"][ctrl["method"]].names[i]
    return "addStandardGroup"
end

function TUIaddStandardsByPrefix!(ctrl::AbstractDict,
                                  response::AbstractString)
    refmat = ctrl["cache"]
    ctrl["method"].standards[refmat] = response
    setGroup!(ctrl["run"],response,refmat)
    ctrl["priority"]["standards"] = false
    return "xxx"
end

function TUIaddStandardsByNumber!(ctrl::AbstractDict,
                                  response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,ctrl["cache"])
    ctrl["priority"]["standards"] = false
    return "xxx"
end

function TUIremoveAllStandards!(ctrl::AbstractDict)
    standardnames = keys(ctrl["method"])
    for samp in ctrl["run"]
        if samp.group in standardnames
            setGroup!(samp,"sample")
        end
    end
    ctrl["method"].standards = Dict()
    ctrl["priority"]["standards"] = true
    return "x"
end

function TUIremoveStandardsByNumber!(ctrl::AbstractDict,
                                     response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,"sample")
    groups = unique(getGroups(ctrl["run"]))
    for k in keys(ctrl["standards"])
        if !in(k,groups)
            delete!(ctrl["method"].standards,k)
        end
    end
    ctrl["priority"]["standards"] = length(ctrl["standards"]) < 1
    return "xxx"
end

function TUIrefmatTab(ctrl::AbstractDict)
    for (key, value) in _KJ["refmat"][ctrl["method"]].dict
        print(key)
        print(": ")
        print_refmat_tx(value)
        print("y0=")
        print(value.y0[1])
        print("\n")
    end
    return nothing
end
function print_refmat_tx(entry::NamedTuple)
    if ismissing(entry.tx[1])
        nothing
    elseif entry.type == "isochron"
        print("t=")
        print(entry.tx[1])
        print("Ma, ")
    elseif entry.type == "point"
        print("x0=")
        print(entry.tx[1])
    else
        nothing
    end
end

function TUIchooseGlass!(ctrl::AbstractDict,
                         response::AbstractString)
    i = parse(Int,response)
    glass = _KJ["glass"].names[i]
    ctrl["cache"] = glass
    ctrl["method"].glass = Dict()
    return "addGlassGroup"
end

function TUIaddGlassByPrefix!(ctrl::AbstractDict,
                              response::AbstractString)
    setGroup!(ctrl["run"],response,ctrl["cache"])
    ctrl["glass"][ctrl["cache"]] = response
    ctrl["priority"]["glass"] = false
    return "xxx"
end

function TUIaddGlassByNumber!(ctrl::AbstractDict,
                              response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,ctrl["cache"])
    ctrl["priority"]["glass"] = false
    return "xxx"
end

function TUIremoveAllGlass!(ctrl::AbstractDict)
    glassnames = keys(ctrl["method"].glass)
    for samp in ctrl["run"]
        if samp.group in glassnames
            setGroup!(samp,"sample")
        end
    end
    ctrl["glass"] = Dict()
    ctrl["priority"]["glass"] = true
    return "x"
end

function TUIremoveGlassByNumber!(ctrl::AbstractDict,
                                 response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,"sample")
    groups = unique(getGroups(ctrl["run"]))
    for (k,v) in ctrl["glass"]
        if !in(k,groups)
            delete!(ctrl["glass"],k)
        end
    end
    ctrl["priority"]["glass"] = length(ctrl["glass"])<1
    return "xxx"
end

function TUIglassTab(ctrl::AbstractDict)
    for name in _KJ["glass"].names
        println(name)
    end
    return nothing
end

function TUIviewer(ctrl::AbstractDict)
    TUIplotter(ctrl)
    return "view"
end

function TUItitle(samp::Sample)
    return string(ctrl["i"]) * ". " * samp.sname * " [" * samp.group * "]"
end

function TUIplotter(ctrl::AbstractDict)
    samp = ctrl["run"][ctrl["i"]]
    p = plot(samp,
             ctrl["method"],
             ctrl["fit"];
             den=ctrl["den"],
             transformation=ctrl["transformation"],
             title=TUItitle(samp))
    if !isnothing(ctrl["method"].PAcutoff)
        TUIaddPAline!(p,ctrl["method"].PAcutoff)
    end
    display(p)
    return nothing
end

function TUIaddPAline!(p,cutoff::AbstractFloat)
    ylim = Plots.ylims(p)
    if  sqrt(cutoff) < 1.1*ylim[2]
        Plots.plot!(p,collect(Plots.xlims(p)),
                    fill(sqrt(cutoff),2),
                    seriestype=:line,label="",
                    linewidth=2,linestyle=:dash)
    end
    return nothing
end

function TUInext!(ctrl::AbstractDict)
    ctrl["i"] += 1
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    return TUIplotter(ctrl)
end

function TUIprevious!(ctrl::AbstractDict)
    ctrl["i"] -= 1
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    return TUIplotter(ctrl)
end

function TUIgoto!(ctrl::AbstractDict,
                  response::AbstractString)
    ctrl["i"] = parse(Int,response)
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    TUIplotter(ctrl)
    return "x"
end

function TUIratios!(ctrl::AbstractDict,
                    response::AbstractString)
    if response=="n"
        ctrl["den"] = nothing
    elseif response=="x"
        return "xx"
    else
        i = parse(Int,response)
        channels = getChannels(ctrl["run"])
        ctrl["den"] = channels[i]
    end
    TUIplotter(ctrl)
    return "x"
end

function TUIt0AutoOne!(ctrl::AbstractDict)
    samp = ctrl["run"][ctrl["i"]]
    sett0!(samp)
    setBwin!(samp)
    setSwin!(samp)
    TUIplotter(ctrl)
    return "x"
end

function TUIt0One!(ctrl::AbstractDict,
                   response::AbstractString)
    t0 = parse(Float64,response)
    samp = ctrl["run"][ctrl["i"]]
    sett0!(samp,t0)
    setBwin!(samp)
    setSwin!(samp)
    TUIplotter(ctrl)
    return "xx"
end

function TUIt0AutoAll!(ctrl::AbstractDict)
    sett0!(ctrl["run"])
    setBwin!(ctrl["run"])
    setSwin!(ctrl["run"])
    TUIplotter(ctrl)
    return "x"
end

function TUIt0All!(ctrl::AbstractDict,
                   response::AbstractString)
    t0 = parse(Float64,response)
    sett0!(ctrl["run"],t0)
    setBwin!(ctrl["run"])
    setSwin!(ctrl["run"])
    TUIplotter(ctrl)
    return "xx"
end

function TUIoneAutoBlankWindow!(ctrl::AbstractDict)
    setBwin!(ctrl["run"][ctrl["i"]])
    return TUIplotter(ctrl)
end

function TUIoneBlankWindow!(ctrl::AbstractDict,
                            response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    bwin = string2windows(samp,response,true)
    setBwin!(samp,bwin)
    TUIplotter(ctrl)
    return "xx"
end

function TUIallAutoBlankWindows!(ctrl::AbstractDict)
    setBwin!(ctrl["run"])
    return TUIplotter(ctrl)
end

function TUIallBlankWindows!(ctrl::AbstractDict,
                              response::AbstractString)
    for samp in ctrl["run"]
        bwin = string2windows(samp,response,true)
        setBwin!(samp,bwin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUIoneAutoSignalWindow!(ctrl::AbstractDict)
    setSwin!(ctrl["run"][ctrl["i"]])
    return TUIplotter(ctrl)
end

function TUIoneSignalWindow!(ctrl::AbstractDict,
                             response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    swin = string2windows(samp,response,true)
    setSwin!(samp,swin)
    TUIplotter(ctrl)
    return "xx"
end

function TUIallAutoSignalWindows!(ctrl::AbstractDict)
    setSwin!(ctrl["run"])
    return TUIplotter(ctrl)
end

function TUIallSignalWindows!(ctrl::AbstractDict,
                              response::AbstractString)
    for samp in ctrl["run"]
        swin = string2windows(samp,response,true)
        setSwin!(samp,swin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUImoveWin!(ctrl::AbstractDict,
                     response::AbstractString)
    shift_windows!(ctrl["run"],parse(Float64,response))
    TUIplotter(ctrl)
    return "x"
end

function TUItransformation!(ctrl::AbstractDict,
                            response::AbstractString)
    if response=="L"
        ctrl["transformation"] = "log"
    elseif response=="s"
        ctrl["transformation"] = "sqrt"
    else
        ctrl["transformation"] = nothing
    end
    TUIplotter(ctrl)
    return "x"
end

function TUIprocess!(ctrl::AbstractDict)
    println("Fitting the model...")
    ctrl["fit"] = process!(ctrl["run"],ctrl["method"])
    ctrl["priority"]["process"] = false
    println("Done")
    return nothing
end

function TUIexport2csv(ctrl::AbstractDict,
                       response::AbstractString)
    tab, out = TUIexportHelper(ctrl["run"],ctrl["method"],ctrl["fit"])
    fname = splitext(response)[1]*".csv"
    CSV.write(fname,tab)
    return out
end

function TUIexportHelper(run::Vector{Sample},method::Gmethod,fit::Gfit)
    tab = averat(ctrl["run"],ctrl["method"],ctrl["fit"])
    out = "xx"
    return tab, out
end
function TUIexportHelper(run::Vector{Sample},method::Cmethod,fit::Cfit)
    tab = concentrations(ctrl["run"],ctrl["method"],ctrl["fit"])
    out = "x"
    return tab, out
end

function TUIexport2json(ctrl::AbstractDict,
                        response::AbstractString)
    ratios = averat(ctrl["run"],ctrl["method"],ctrl["fit"])
    fname = splitext(response)[1]*".json"
    export2IsoplotR(ratios,ctrl["method"];fname=fname)
    return "xx"
end

function TUIimportLog!(ctrl::AbstractDict,
                       response::AbstractString;
                       verbose::Bool=false)
    TUIclear!(ctrl)
    ctrl["log"] = true
    history = CSV.read(response,DataFrame)
    for row in eachrow(history)
        try
            if verbose println(row) end
            dispatch!(ctrl;key=row[1],response=row[2],verbose=verbose)
        catch e
            println(e)
        end
    end
    ctrl["log"] = false
    return nothing
end

function TUIexportLog(ctrl::AbstractDict,
                      response::AbstractString)
    ctrl["history"] = ctrl["history"][1:end-1,:]
    CSV.write(response,ctrl["history"])
    return "xx"
end

function TUIopenTemplate!(ctrl::AbstractDict,
                          response::AbstractString)
    include(response)
    ctrl["format"] = format
    ctrl["head2name"] = head2name
    ctrl["multifile"] = multifile
    ctrl["method"] = method
    ctrl["priority"]["standards"] = (length(method.standards) == 0)
    ctrl["transformation"] = transformation
    ctrl["priority"]["method"] = false
    ctrl["template"] = true
    return "xx"
end

function TUIsaveTemplate(ctrl::AbstractDict,
                         response::AbstractString)
    PAcutoff = isnothing(ctrl["PAcutoff"]) ? "nothing" : string(ctrl["PAcutoff"])
    open(response, "w") do file
        write(file,"format = \"" * ctrl["format"] * "\"\n")
        write(file,"multifile = " * string(ctrl["multifile"]) * "\n")
        write(file,"head2name = " * string(ctrl["head2name"]) * "\n")
        write(file,"method = \"" * ctrl["method"] * "\"\n")
        write(file,"options = " * dict2string(ctrl["options"]) * "\n")
        write(file,"PAcutoff = " * PAcutoff * "\n")
        write(file,"transformation = \"" * ctrl["transformation"] * "\"\n")
        if length(ctrl["glass"])>0
            write(file,"glass = " * dict2string(ctrl["glass"]) * "\n")
        end
        if length(ctrl["standards"])>0
            write(file,"standards = " * dict2string(ctrl["standards"]) * "\n")
        end
        if ctrl["method"] == "concentrations"
            write(file,"channels = " * vec2string(ctrl["channels"]) * "\n")
        else
            write(file,"channels = " * dict2string(ctrl["channels"]) * "\n")
        end
        if isnothing(ctrl["internal"])
            write(file,"internal = nothing\n")
        else
            write(file,"internal = (\"" *
                  ctrl["internal"][1] * "\"," *
                  string(ctrl["internal"][2]) * ")")
        end
    end
    return "xx"
end

function TUIsetNblank!(ctrl::AbstractDict,
                       response::AbstractString)
    ctrl["method"].nblank = parse(Int,response)
    return "x"
end

function TUIsetNdrift!(ctrl::AbstractDict,
                       response::AbstractString)
    ctrl["method"].ndrift = parse(Int,response)
    return "x"    
end

function TUIsetNdown!(ctrl::AbstractDict,
                      response::AbstractString)
    ctrl["method"].ndown = parse(Int,response)
    return "x"    
end

function TUIPAlist(ctrl::AbstractDict)
    snames = getSnames(ctrl["run"])
    for i in eachindex(snames)
        dat = getSignals(ctrl["run"][i],ctrl["channels"])
        maxval = maximum(Matrix(dat))
        formatted = @sprintf("%.*e", 3, maxval)
        println(formatted*" ("*snames[i]*")")
    end
    return "x"
end

function TUIsetPAcutoff!(ctrl::AbstractDict,
                         response::AbstractString)
    cutoff = tryparse(Float64,response)
    ctrl["method"].PAcutoff = cutoff
    return "xx"
end

function TUIclearPAcutoff!(ctrl::AbstractDict)
    ctrl["method"].PAcutoff = nothing
    return "xx"
end

function TUIaddStandard!(ctrl::AbstractDict,
                         response::AbstractString)
    init_referenceMaterials!(response)
    return "x"
end

function TUIaddGlass!(ctrl::AbstractDict,
                      response::AbstractString)
    init_glass!(response)
    return "x"
end

function TUIhead2name!(ctrl::AbstractDict,
                       response::AbstractString)
    ctrl["head2name"] = response=="h"
    return "x"
end

function TUIrefresh!(ctrl::AbstractDict)
    if ctrl["multifile"]
        TUIloadICPdir!(ctrl,ctrl["ICPpath"])
    else
        TUIloadICPfile!(ctrl,ctrl["ICPpath"])
        TUIloadLAfile!(ctrl,ctrl["LApath"])
    end
    TUIsetGroups!(ctrl,"standards")
    TUIsetGroups!(ctrl,"glass")
    TUIprocess!(ctrl)
    return nothing
end

function TUIsetGroups!(ctrl::AbstractDict,std::AbstractString)
    for (refmat,prefix) in ctrl[std]
        if !isnothing(prefix)
            setGroup!(ctrl["run"],prefix,refmat)
        end
    end
end

function TUIclear!(ctrl::AbstractDict)
    default = TUIinit()
    for (k,v) in default
        ctrl[k] = v
    end
    for (i, extension) in enumerate(_KJ["extensions"])
        extension.extend!(_KJ)
    end
    return nothing
end

function TUInternochronViewer!(ctrl::AbstractDict)
    TUInternochron(ctrl)
    return "internoview"
end

function TUInternochron_next!(ctrl::AbstractDict)
    ctrl["i"] += 1
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    return TUInternochron(ctrl)
end

function TUInternochron_previous!(ctrl::AbstractDict)
    ctrl["i"] -= 1
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    return TUInternochron(ctrl)
end

function TUInternochron_goto!(ctrl::AbstractDict,
                              response::AbstractString)
    ctrl["i"] = parse(Int,response)
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    TUInternochron(ctrl)
    return "x"
end

function TUInternochron(ctrl::AbstractDict)
    p = internoplot(ctrl["run"][ctrl["i"]],ctrl["method"],ctrl["fit"];
                    title = string(ctrl["i"])*". "*samp.sname*" ["*samp.group*"]")
    display(p)
    return nothing
end

function internochron2csv(ctrl::AbstractDict,
                          fname::AbstractString)
    tab = internochron(ctrl["run"],ctrl["method"],ctrl["fit"])
    CSV.write(fname,tab)
    return "xx"
end

function CSVhelper(samp::Sample,method::Gmethod,fit::Gfit)
    return atomic(samp,method,fit;add_xy=true)
end
function CSVhelper(samp::Sample,method::Cmethod,fit::Cfit)
    return concentrations(samp,method,fit)
end

function TUItimeresolved2csv(ctrl::AbstractDict,
                             dname::AbstractString)
    for samp in ctrl["run"]
        fname = samp.sname * ".csv"
        path = joinpath(dname,fname)
        df = CSVhelper(samp,ctrl["method"],ctrl["fit"])
        CSV.write(path,df)
    end
    return "x"
end

function TUImapper(ctrl::AbstractDict)
    TUImap(ctrl)
    return "map"
end

function TUImap(ctrl::AbstractDict)
    ctrl["cache"] = CSVhelper(ctrl["run"][ctrl["i"]],
                              ctrl["method"],ctrl["fit"])
    colnames = names(ctrl["cache"])
    selected_column = colnames[ctrl["mapcolumn"]]
    p = plotMap(ctrl["cache"],
                selected_column;
                clims=ctrl["clims"])
    tit = string(ctrl["i"]) * ". " * ctrl["run"][ctrl["i"]].sname
    Plots.title!(p, tit)
    display(p)
    return nothing
end

function TUIchooseMapColumn!(ctrl::AbstractDict,
                             response::AbstractString)
    ctrl["mapcolumn"] = parse(Int,response)
    TUImap(ctrl)
    return "x"
end

function TUImapPrevious!(ctrl::AbstractDict)
    ctrl["i"] -= 1
    if ctrl["i"]<1 ctrl["i"] = length(ctrl["run"]) end
    return TUImap(ctrl)
end

function TUImapNext!(ctrl::AbstractDict)
    ctrl["i"] += 1
    if ctrl["i"]>length(ctrl["run"]) ctrl["i"] = 1 end
    return TUImap(ctrl)
end

function TUIchooseClims!(ctrl::AbstractDict,
                         response::AbstractString)
    if response == "r"
        ctrl["clims"] = nothing
    else
        parts = split(response,',')
        if length(parts)>1 && all(isdigit,parts[1]) && all(isdigit,parts[2])
            min = parse(Float64,parts[1])
            max = parse(Float64,parts[2])
            ctrl["clims"] = (min,max)
        end
    end
    TUImap(ctrl)
    return "x"
end
