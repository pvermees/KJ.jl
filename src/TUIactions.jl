function TUIinit!()
    _KJ["ctrl"] = TUIinit()
    return nothing
end

function TUIinit()
    return Dict(
        "priority" => Dict("load" => true, 
                           "method" => true,
                           "fractionation" => true,
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
        "den" => "",
        "transformation" => "sqrt",
        "mapcolumn" => 2,
        "clims" => (),
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
    TUIsetGroups!(ctrl)
    if ctrl["template"]
        return "x"
    else
        return "xxx"
    end
end

function TUIsetGroups!(ctrl::AbstractDict)
    if !isnothing(ctrl["method"])
        setGroup!(ctrl["run"],ctrl["method"])
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
        ctrl["priority"]["fractionation"] = false
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
        ctrl["method"] = Cmethod(ctrl["run"])
        return "internal"
    else
        i = parse(Int,response)
        methodname = _KJ["methods"].names[i]
        ctrl["method"] = Gmethod(name=methodname)
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
    ctrl["cache"] = getChannels(ctrl["run"])[i]
    return "mineral"
end

function TUIstoichiometry!(ctrl::AbstractDict,
                           response::AbstractString)
    channel = ctrl["cache"]
    concentration = parse(Float64,response)
    ctrl["method"].internal = (channel,concentration)
    ctrl["priority"]["method"] = false
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
        return "xxx"
    end
end

function TUIsetChannels!(ctrl::AbstractDict,
                         response::AbstractString)
    labels = getChannels(ctrl["run"])
    selected = parse.(Int,split(response,","))
    channels = labels[selected[1:3]]
    ctrl["method"].P.channel = channels[1]
    ctrl["method"].D.channel = channels[2]
    ctrl["method"].d.channel = channels[3]
    proxies = channel2proxy.(channels)
    if any(isnothing(proxies))
        return "setProxies"
    else
        ctrl["method"].P.proxy = proxies[1]
        ctrl["method"].D.proxy = proxies[2]
        ctrl["method"].d.proxy = proxies[3]
        ctrl["priority"]["method"] = false
        return "xx"
    end
end

function TUIgetChannels(ctrl::AbstractDict)
    if isnothing(ctrl["method"])
        channels = getChannels(ctrl["run"])
    else
        channels = getChannels(ctrl["method"])
    end
    return channels
end

function TUIsetProxies!(ctrl::AbstractDict,
                        response::AbstractString)
    selection = parse.(Int,split(response,","))
    ions = ctrl["method"].ions
    nuclides = ions2nuclidelist(ions)
    pr = ctrl["method"].proxies
    pr.P, pr.D, pr.d = nuclides[selection[1:3]]
    return "xxx"
end

function TUIchooseStandard!(ctrl::AbstractDict,
                            response::AbstractString)
    i = parse(Int,response)
    ctrl["cache"] = TUIgetRefmatName(ctrl["method"],i)
    return "addStandardGroup"
end

function TUIgetRefmatName(method::Gmethod,i::Int)
    return _KJ["refmat"][method.name].names[i]
end

function TUIgetRefmatName(method::Cmethod,i::Int)
    return _KJ["glass"].names[i]
end

function TUIsetStandard!(ctrl::AbstractDict,
                         group::AbstractString,
                         standard::AbstractString)
    push!(ctrl["method"].standards,group)
    ctrl["method"].groups[group] = standard
end

function TUIaddStandardsByPrefix!(ctrl::AbstractDict,
                                  response::AbstractString)
    TUIsetStandard!(ctrl,response,ctrl["cache"])
    setGroup!(ctrl["run"],ctrl["method"])
    ctrl["priority"]["fractionation"] = false
    return "xxx"
end

function TUIaddStandardsByNumber!(ctrl::AbstractDict,
                                  response::AbstractString)
    selection = parse.(Int,split(response,","))
    TUIsetStandard!(ctrl,ctrl["cache"],ctrl["cache"])
    setGroup!(ctrl["run"],selection,ctrl["cache"])
    ctrl["priority"]["fractionation"] = false
    return "xxx"
end

function TUIremoveAllStandards!(ctrl::AbstractDict)
    for group in ctrl["method"].standards
        delete!(ctrl["method"].groups,group)
    end
    setGroup!(ctrl["run"])
    setGroup!(ctrl["run"],ctrl["method"])
    ctrl["method"].standards = Set{String}()
    ctrl["priority"]["fractionation"] = true
    return "x"
end

function TUIremoveStandardsByNumber!(ctrl::AbstractDict,
                                     response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection)
    groups = unique(getGroups(ctrl["run"]))
    standards = ctrl["method"].standards
    for standard in standards
        if !in(standard,groups)
            delete!(ctrl["method"].groups,standard)
            pop!(standards,standard)
        end
    end
    ctrl["priority"]["fractionation"] = length(standards) < 1
    return "xx"
end

function TUIstandardsTab(ctrl::AbstractDict)
    return TUIstandardsTab(ctrl["method"])
end
function TUIstandardsTab(method::Cmethod)
    refmats = TUIgetStandards(method)
    for name in refmats.names
        print(name * "\n")
    end
    return nothing
end
function TUIstandardsTab(method::Gmethod)
    refmats = TUIgetStandards(method)
    for name in refmats.names
        refmat = get(refmats,name)
        print(name)
        print(" (")
        print(refmat.material)
        print("): ")
        print_refmat(refmat)
        print("\n")
    end
    return nothing
end
function print_refmat(entry::IsochronRefmat)
        print("t=")
        print(entry.t)
        print("Ma, ")
        print("y0=")
        print(entry.y0)
end
function print_refmat(entry::PointRefmat)
        print("x=")
        print(entry.x)
        print(", ")
        print("y=")
        print(entry.y)
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
    return "xxx"
end

function TUIaddGlassByNumber!(ctrl::AbstractDict,
                              response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,ctrl["cache"])
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

function TUItitle(ctrl::AbstractDict)
    samp = ctrl["run"][ctrl["i"]]
    return string(ctrl["i"]) * ". " * samp.sname * " [" * samp.group * "]"
end

function TUIplotter(ctrl::AbstractDict)
    p = plot(ctrl["run"][ctrl["i"]],
             ctrl["method"];
             fit=ctrl["fit"],
             den=ctrl["den"],
             transformation=ctrl["transformation"],
             title=TUItitle(ctrl))
    if ctrl["method"] isa Gmethod && !isnothing(ctrl["method"].PAcutoff)
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
                    linewidth=1,linestyle=:dash,
                    linecolor="black",
                    labels="P/A cutoff")
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
        ctrl["den"] = ""
    elseif response=="x"
        return "xx"
    else
        i = parse(Int,response)
        channels = TUIgetChannels(ctrl)
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
    TUIwindowHandler!(ctrl;all=false,single=true,blank=true)
end
function TUIallAutoBlankWindows!(ctrl::AbstractDict)
    TUIwindowHandler!(ctrl;all=true,single=true,blank=true)
end
function TUIoneAutoSignalWindow!(ctrl::AbstractDict)
    TUIwindowHandler!(ctrl;all=false,single=true,blank=false)
end
function TUIallAutoSignalWindows!(ctrl::AbstractDict)
    TUIwindowHandler!(ctrl;all=true,single=true,blank=false)
end
function TUIoneSingleBlankWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=false,single=true,blank=true)
end
function TUIoneMultiBlankWindow!(ctrl::AbstractDict,
                                 response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=false,single=false,blank=true)
end
function TUIallSingleBlankWindows!(ctrl::AbstractDict,
                                   response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=true,single=true,blank=true)
end
function TUIallMultiBlankWindows!(ctrl::AbstractDict,
                                  response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=true,single=false,blank=true)
end
function TUIoneSingleSignalWindow!(ctrl::AbstractDict,
                                   response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=false,single=true,blank=false)
end
function TUIoneMultiSignalWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=false,single=false,blank=false)
end
function TUIallSingleSignalWindows!(ctrl::AbstractDict,
                                    response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=true,single=true,blank=false)
end
function TUIallMultiSignalWindows!(ctrl::AbstractDict,
                                   response::AbstractString)
    TUIwindowHandler!(ctrl;response=response,all=true,single=false,blank=false)
end
function TUIwindowHandler!(ctrl::AbstractDict;
                           response::AbstractString="",
                           all::Bool=false,
                           single::Bool=false,
                           blank::Bool=false)
    target = ifelse(all,ctrl["run"],ctrl["run"][ctrl["i"]])
    fun! = ifelse(blank,setBwin!,setSwin!)
    if response==""
        fun!(target)
        next = "x"
    else
        win = string2windows(target,response,single)
        fun!(target,win)
        next = "xx"
    end
    TUIplotter(ctrl)
    return next
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
    tab, out = TUIexportHelper(ctrl["run"], ctrl["method"], ctrl["fit"])
    fname = splitext(response)[1]*".csv"
    CSV.write(fname,tab)
    return out
end

function TUIexportHelper(run::Vector{Sample},method::Gmethod,fit::Gfit)
    tab = averat(run,method,fit)
    out = "xx"
    return tab, out
end
function TUIexportHelper(run::Vector{Sample},method::Cmethod,fit::Cfit)
    tab = concentrations(run,method,fit)
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
    Base.invokelatest() do
        ctrl["refmats"] = refmats
        ctrl["format"] = format
        ctrl["head2name"] = head2name
        ctrl["multifile"] = multifile
        ctrl["method"] = method
        ctrl["transformation"] = transformation
        ctrl["priority"]["fractionation"] = (length(method.fractionation.standards) == 0)
        ctrl["priority"]["method"] = false
        ctrl["template"] = true
    end
    return "xx"
end

function TUIsaveTemplate(ctrl::AbstractDict,
                         response::AbstractString)
    open(response, "w") do file
        write(file,"format = \"" * ctrl["format"] * "\"\n")
        write(file,"multifile = " * string(ctrl["multifile"]) * "\n")
        write(file,"head2name = " * string(ctrl["head2name"]) * "\n")
        write(file,"transformation = \"" * ctrl["transformation"] * "\"\n")
        write(file,"refmats = " * TUIrefmats2text(ctrl["refmats"]) * "\n") 
        write(file,TUImethod2text(ctrl["method"]))
    end
    return "xx"
end

function TUImethod2text(method::Gmethod)
    F = method.fractionation
    i = F.ions
    p = F.proxies
    c = F.channels
    s = collect(F.standards)
    PA = method.PAcutoff
    out = "method = Gmethod(\"" * method.name * "\";\n"
    out *= "                 ions = (P=\"" * i.P * "\", D=\"" * i.D * "\", d=\"" * i.d * "\"),\n"
    out *= "                 proxies = (P=\"" * p.P * "\", D=\"" * p.D * "\", d=\"" * p.d * "\"),\n"
    out *= "                 channels = (P=\"" * c.P * "\", D=\"" * c.D * "\", d=\"" * c.d * "\"),\n"
    out *= "                 standards = [" * join("\"" .* s .* "\"", ", ") * "],\n"
    out *= "                 nblank = " * string(method.nblank) * ",\n"
    out *= "                 ndrift = " * string(method.ndrift) * ",\n"
    out *= "                 ndown = " * string(method.ndown) * ",\n"
    out *= "                 PAcutoff = " * ifelse(isnothing(PA),"nothing",PA) * ")\n"
    return out
end

function TUImethod2text(method::Cmethod)
    chunks = String[]
    for (channel,element) in pairs(method.elements)
        push!(chunks,"\"" * channel * "\" => " * "\"" * element * "\"")
    end
    out  = "elements = DataFrame(" * join(chunks,",\n                     ") * ")\n"
    out *= "standards = [" * join("\"" .* collect(method.standards) .* "\"", ", ") * "],\n"
    out *= "internal = (\"" * method.internal[1] * "\"," * string(method.internal[2]) * ")\n"
    out *= "method = Cmethod(elements,standards,internal," * string(method.nblank) * ")"
    return out
end

function TUIrefmats2text(refmats::AbstractDict)
    chunks = String[]
    for (k,v) in refmats
        push!(chunks, "\"" * k * "\" => \"" * v * "\"")
    end
    out = "Dict(" * join(chunks,",") * ")"
    return out
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
    for i in eachindex(ctrl["run"])
        samp = ctrl["run"][i]
        channels = getChannels(ctrl["method"])
        dat = getSignals(samp)[:,channels]
        maxval = maximum(Matrix(dat))
        formatted = @sprintf("%.*e", 3, maxval)
        println(formatted*" ("*samp.sname*")")
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
    init_referenceMaterials!(;isochrons=response)
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
    nold = length(ctrl["run"])
    if ctrl["ICPpath"] == ""
        return nothing
    end
    if ctrl["multifile"]
        TUIloadICPdir!(ctrl,ctrl["ICPpath"])
    else
        TUIloadICPfile!(ctrl,ctrl["ICPpath"])
        TUIloadLAfile!(ctrl,ctrl["LApath"])
    end
    nnew = length(ctrl["run"])
    if nnew > nold
        setGroup!(ctrl["run"][nold+1:nn],ctrl["refmats"])
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
    samp = ctrl["run"][ctrl["i"]]
    p = internoplot(samp,ctrl["method"],ctrl["fit"];
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

function CSVhelper(samp::Sample,
                   method::Gmethod,
                   fit::Gfit)
    a = atomic(samp,method,fit;add_xy=true)
    i = method.ions
    out = DataFrame(i.P => a.P, i.D => a.D, i.d => a.d)
    if !isnothing(a.x) out.x = x end
    if !isnothing(a.y) out.y = y end
    return out
end
function CSVhelper(samp::Sample,
                   method::Cmethod,
                   fit::Cfit)
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

function TUItodo!(ctrl::AbstractDict)
    println("This function has not yet been implemented, " * 
            "because KJ is still in development. " * 
            "Please check again later.")
    return nothing
end