function TUIinit!()
    _KJ["ctrl"] = TUIinit()
    return nothing
end
function TUIinit()
    return Dict(
        "priority" => Dict("load" => true, "method" => true,
                           "standards" => true, "glass" => true,
                           "process" => true),
        "history" => DataFrame(task=String[],action=String[]),
        "chain" => ["top"],
        "run" => nothing,
        "i" => 1,
        "den" => nothing,
        "multifile" => true,
        "head2name" => true,
        "method" => "",
        "instrument" => "",
        "ICPpath" => "",
        "LApath" => "",
        "channels" => nothing,
        "standards" => Dict{String,Union{Nothing,String}}(),
        "glass" => Dict{String,Union{Nothing,String}}(),
        "internal" => nothing,
        "options" => Dict("blank" => 2, "drift" => 1, "down" => 1),
        "PAcutoff" => nothing,
        "blank" => nothing,
        "par" => nothing,
        "cache" => nothing,
        "transformation" => "sqrt",
        "log" => false,
        "template" => false
    )
end

function TUI(key::AbstractString)
    return _KJ["ctrl"][key]
end
export TUI

function TUIwelcome()
    version = string(pkgversion(@__MODULE__))
    title = " Plasmatrace "*version*" \n"
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
        return "instrument"
    end
end

function TUIinstrument!(ctrl::AbstractDict,
                        response::AbstractString)
    if response=="a"
        ctrl["instrument"] = "Agilent"
    elseif response=="t"
        ctrl["instrument"] = "ThermoFisher"
    else
        @warn "Unsupported instrument"
        return "x"
    end
    return "dir|file"
end

function TUIloadICPdir!(ctrl::AbstractDict,
                        response::AbstractString)
    ctrl["run"] = load(response;
                       instrument=ctrl["instrument"],
                       head2name=ctrl["head2name"])
    if isnothing(ctrl["channels"])
        ctrl["channels"] = getChannels(ctrl["run"])
    end
    ctrl["priority"]["load"] = false
    ctrl["multifile"] = true
    ctrl["ICPpath"] = response
    println("Done")
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
    ctrl["run"] = load(ctrl["ICPpath"],ctrl["LApath"];
                       instrument=ctrl["instrument"])
    if isnothing(ctrl["channels"])
        ctrl["channels"] = getChannels(ctrl["run"]) # reset
    end
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

function TUImethod!(ctrl::AbstractDict,
                    response::AbstractString)
    methods = _KJ["methods"].method
    if response=="c"
        ctrl["method"] = "concentrations"
        return "internal"
    else
        i = parse(Int,response)
        if i > length(methods)
            return "x"
        else
            ctrl["method"] = methods[i]
            return "columns"
        end
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
    ctrl["internal"] = (channel,concentration)
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
        mineral = collect(keys(_KJ["stoichiometry"]))[i]
        channel = ctrl["cache"]
        ctrl["internal"] = getInternal(mineral,channel)
        ctrl["priority"]["method"] = false
        ctrl["priority"]["standards"] = false
        return "xxx"
    end
end

function TUIcolumns!(ctrl::AbstractDict,
                     response::AbstractString)
    labels = names(getSignals(ctrl["run"][1]))
    selected = parse.(Int,split(response,","))
    PDd = labels[selected]
    ctrl["channels"] = Dict("d" => PDd[3], "D" => PDd[2], "P" => PDd[1])
    ctrl["priority"]["method"] = false
    return "xx"
end

function TUIchooseStandard!(ctrl::AbstractDict,
                            response::AbstractString)
    i = parse(Int,response)
    standards = collect(keys(_KJ["refmat"][ctrl["method"]]))
    ctrl["cache"] = standards[i]
    ctrl["standards"][standards[i]] = nothing
    return "addStandardGroup"
end

function TUIaddStandardsByPrefix!(ctrl::AbstractDict,
                                  response::AbstractString)
    setGroup!(ctrl["run"],response,ctrl["cache"])
    ctrl["standards"][ctrl["cache"]] = response
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
    empty!(ctrl["standards"])
    setGroup!(ctrl["run"],"sample")
    ctrl["priority"]["standards"] = true
    return "xx"
end

function TUIremoveStandardsByNumber!(ctrl::AbstractDict,
                                     response::AbstractString)
    selection = parse.(Int,split(response,","))
    setGroup!(ctrl["run"],selection,"sample")
    groups = unique(getGroups(ctrl["run"]))
    for (k,v) in ctrl["standards"]
        if !in(k,groups)
            delete!(ctrl["standards"],k)
        end
    end
    ctrl["priority"]["standards"] = length(ctrl["standards"])<1
    return "xxx"
end

function TUIrefmatTab(ctrl::AbstractDict)
    for (key, value) in _KJ["refmat"][ctrl["method"]]
        print(key)
        print(": ")
        if !ismissing(value.t[1])
            print("t=")
            print(value.t[1])
            print("Ma, ")
        end
        print("y0=")
        print(value.y0[1])
        print("\n")
    end
    return "x"
end

function TUIchooseGlass!(ctrl::AbstractDict,
                         response::AbstractString)
    i = parse(Int,response)
    glass = collect(keys(_KJ["glass"]))
    ctrl["cache"] = glass[i]
    ctrl["glass"][glass[i]] = nothing
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
    empty!(ctrl["standards"])
    setGroup!(ctrl["run"],"sample")
    ctrl["priority"]["glass"] = true
    return "xx"
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
    for (key, value) in _KJ["glass"]
        println(key)
    end
    return "x"
end

function TUIviewer!(ctrl::AbstractDict)
    TUIplotter(ctrl)
    return "view"
end

function TUIplotter(ctrl::AbstractDict)
    samp = ctrl["run"][ctrl["i"]]
    if ctrl["method"] == "concentrations"
        p = TUIconcentrationPlotter(ctrl,samp)
    else
        p = TUIgeochronPlotter(ctrl,samp)
    end
    if !isnothing(ctrl["PAcutoff"])
        TUIaddPAline!(p,ctrl["PAcutoff"])
    end
    display(p)
    return nothing
end

function TUIconcentrationPlotter(ctrl::AbstractDict,samp::Sample)
    if (samp.group in keys(ctrl["glass"])) & !isnothing(ctrl["blank"])
        p = plot(samp,ctrl["blank"],ctrl["par"],ctrl["internal"][1];
                 den=ctrl["den"],transformation=ctrl["transformation"],i=ctrl["i"])
    else
        p = plot(samp;den=ctrl["den"],transformation=ctrl["transformation"],i=ctrl["i"])
    end
    return p
end

function TUIgeochronPlotter(ctrl::AbstractDict,samp::Sample)
    if isnothing(ctrl["blank"]) | (samp.group=="sample")
        p = plot(samp,ctrl["channels"];
                 den=ctrl["den"],transformation=ctrl["transformation"],i=ctrl["i"])
    else
        anchors = getAnchors(ctrl["method"],ctrl["standards"],ctrl["glass"])
        p = plot(samp,ctrl["method"],ctrl["channels"],ctrl["blank"],
                 ctrl["par"],ctrl["standards"],ctrl["glass"];
                 den=ctrl["den"],transformation=ctrl["transformation"],i=ctrl["i"])
    end
    return p
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
        if isa(ctrl["channels"],AbstractVector)
            channels = ctrl["channels"]
        elseif isa(ctrl["channels"],AbstractDict)
            channels = collect(values(ctrl["channels"]))
        else
            channels = getChannels(ctrl["run"])
        end
        ctrl["den"] = channels[i]
    end
    TUIplotter(ctrl)
    return "x"
end

function TUIoneAutoBlankWindow!(ctrl::AbstractDict)
    setBwin!(ctrl["run"][ctrl["i"]])
    return TUIplotter(ctrl)
end

function TUIoneSingleBlankWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    bwin = string2windows(samp,response,true)
    setBwin!(samp,bwin)
    TUIplotter(ctrl)
    return "xx"
end

function TUIoneMultiBlankWindow!(ctrl::AbstractDict,
                                 response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    bwin = string2windows(samp,response,false)
    setBwin!(samp,bwin)
    TUIplotter(ctrl)
    return "xx"
end

function TUIallAutoBlankWindow!(ctrl::AbstractDict)
    setBwin!(ctrl["run"])
    return TUIplotter(ctrl)
end

function TUIallSingleBlankWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    for i in eachindex(ctrl["run"])
        samp = ctrl["run"][i]
        bwin = string2windows(samp,response,true)
        setBwin!(samp,bwin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUIallMultiBlankWindow!(ctrl::AbstractDict,
                                 response::AbstractString)
    for i in eachindex(ctrl["run"])
        samp = ctrl["run"][i]
        bwin = string2windows(samp,response,false)
        setBwin!(samp,bwin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUIoneAutoSignalWindow!(ctrl::AbstractDict)
    setSwin!(ctrl["run"][ctrl["i"]])
    return TUIplotter(ctrl)
end

function TUIoneSingleSignalWindow!(ctrl::AbstractDict,
                                   response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    swin = string2windows(samp,response,true)
    setSwin!(samp,swin)
    TUIplotter(ctrl)
    return "xx"
end

function TUIoneMultiSignalWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    samp = ctrl["run"][ctrl["i"]]
    swin = string2windows(samp,response,false)
    setSwin!(samp,swin)
    TUIplotter(ctrl)
    return "x"
end

function TUIallAutoSignalWindow!(ctrl::AbstractDict)
    setSwin!(ctrl["run"])
    return TUIplotter(ctrl)
end

function TUIallSingleSignalWindow!(ctrl::AbstractDict,
                                   response::AbstractString)
    for i in eachindex(ctrl["run"])
        samp = ctrl["run"][i]
        swin = string2windows(samp,response,true)
        setSwin!(samp,swin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUIallMultiSignalWindow!(ctrl::AbstractDict,
                                  response::AbstractString)
    for i in eachindex(ctrl["run"])
        samp = ctrl["run"][i]
        swin = string2windows(samp,response,false)
        setSwin!(samp,swin)
    end
    TUIplotter(ctrl)
    return "xx"
end

function TUItransformation!(ctrl::AbstractDict,
                            response::AbstractString)
    if response=="L"
        ctrl["transformation"] = "log"
    elseif response=="s"
        ctrl["transformation"] = "sqrt"
    else
        ctrl["transformation"] = ""
    end
    TUIplotter(ctrl)
    return "x"
end

function TUIprocess!(ctrl::AbstractDict)
    println("Fitting blanks...")
    ctrl["blank"] = fitBlanks(ctrl["run"],nblank=ctrl["options"]["blank"])
    println("Fractionation correction...")
    if ctrl["method"] == "concentrations"
        ctrl["par"] = fractionation(ctrl["run"],
                                    ctrl["blank"],
                                    ctrl["internal"],
                                    ctrl["glass"])
    else
        ctrl["anchors"] = getAnchors(ctrl["method"],ctrl["standards"],ctrl["glass"])
        ctrl["par"] = fractionation(ctrl["run"],
                                    ctrl["method"],
                                    ctrl["blank"],
                                    ctrl["channels"],
                                    ctrl["standards"],
                                    ctrl["glass"];
                                    ndrift=ctrl["options"]["drift"],
                                    ndown=ctrl["options"]["down"],
                                    PAcutoff=ctrl["PAcutoff"])
    end
    ctrl["priority"]["process"] = false
    println("Done")
    return nothing
end

function TUIsubset!(ctrl::AbstractDict,
                    response::AbstractString)
    if response=="a"
        ctrl["cache"] = 1:length(ctrl["run"])
    elseif response=="s"
        ctrl["cache"] = findall(contains("sample"),getGroups(ctrl["run"]))
    elseif response=="x"
        return "x"
    else
        ctrl["cache"] = findall(contains(response),getSnames(ctrl["run"]))
    end
    if ctrl["method"] == "concentrations"
        return "csv"
    else
        return "format"
    end
end

function TUIexport2csv(ctrl::AbstractDict,
                       response::AbstractString)
    if ctrl["method"]=="concentrations"
        out = concentrations(ctrl["run"],ctrl["blank"],ctrl["par"],ctrl["internal"])
    else
        out = averat(ctrl["run"],ctrl["channels"],ctrl["blank"],ctrl["par"])
    end
    fname = splitext(response)[1]*".csv"
    CSV.write(fname,out[ctrl["cache"],:])
    if ctrl["method"] == "concentrations"
        return "xx"
    else
        return "xxx"
    end
end

function TUIexport2json(ctrl::AbstractDict,
                        response::AbstractString)
    ratios = averat(ctrl["run"],ctrl["channels"],ctrl["blank"],ctrl["par"])
    fname = splitext(response)[1]*".json"
    export2IsoplotR(ratios[ctrl["cache"],:],ctrl["method"];fname=fname)
    return "xxx"
end

function TUIimportLog!(ctrl::AbstractDict,
                       response::AbstractString)
    TUIclear!(ctrl)
    ctrl["log"] = true
    history = CSV.read(response,DataFrame)
    for row in eachrow(history)
        try
            dispatch!(ctrl;key=row[1],response=row[2])
        catch e
            println(e)
        end
    end
    ctrl["log"] = false
    return nothing
end

function TUIexportLog(ctrl::AbstractDict,
                      response::AbstractString)
    nh = size(ctrl["history"],1)
    deleteat!(ctrl["history"],[nh-1,nh])
    CSV.write(response,ctrl["history"])
    return "xx"
end

function TUIopenTemplate!(ctrl::AbstractDict,
                          response::AbstractString)
    include(response)
    ctrl["instrument"] = instrument
    ctrl["head2name"] = head2name
    ctrl["multifile"] = multifile
    ctrl["method"] = method
    ctrl["channels"] = channels
    ctrl["options"] = options
    ctrl["PAcutoff"] = PAcutoff
    if @isdefined(standards)
        ctrl["standards"] = standards
        ctrl["priority"]["standards"] = all(isnothing.(values(standards)))
    end
    if @isdefined(glass)
        ctrl["glass"] = glass
        ctrl["priority"]["glass"] = all(isnothing.(values(glass)))
    end
    ctrl["transformation"] = transformation
    ctrl["internal"] = internal
    ctrl["priority"]["method"] = false
    ctrl["template"] = true
    return "xx"
end

function TUIsaveTemplate(ctrl::AbstractDict,
                         response::AbstractString)
    PAcutoff = isnothing(ctrl["PAcutoff"]) ? "nothing" : string(ctrl["PAcutoff"])
    open(response, "w") do file
        write(file,"instrument = \"" * ctrl["instrument"] * "\"\n")
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
    ctrl["options"]["blank"] = parse(Int,response)
    return "x"
end

function TUIsetNdrift!(ctrl::AbstractDict,
                       response::AbstractString)
    ctrl["options"]["drift"] = parse(Int,response)
    return "x"    
end

function TUIsetNdown!(ctrl::AbstractDict,
                      response::AbstractString)
    ctrl["options"]["down"] = parse(Int,response)
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
    ctrl["PAcutoff"] = cutoff
    return "xx"
end

function TUIclearPAcutoff!(ctrl::AbstractDict)
    ctrl["PAcutoff"] = nothing
    return "xx"
end

function TUIaddStandard!(ctrl::AbstractDict,
                         response::AbstractString)
    setReferenceMaterials!(response)
    return "x"
end

function TUIaddGlass!(ctrl::AbstractDict,
                      response::AbstractString)
    setGlass!(response)
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
    snames = getSnames(ctrl["run"])
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
