function TUIwelcomeMessage(ctrl::AbstractDict)
    msg = "r: Read data files"*TUIcheck(ctrl,"load")*"\n"*
    "m: Specify the method"*TUIcheck(ctrl,"method")*"\n"*
    "t: Tabulate the samples\n"*
    "v: View and adjust each sample\n"*
    TUImethodSwitch(ctrl["method"],"i: Interferences\n","")*
    TUImethodSwitch(ctrl["method"],
                    "f: Fractionation",
                    "g: Glass")*TUIcheck(ctrl,"fractionation")*"\n"*
    TUImethodSwitch(ctrl["method"],"b: Mass bias\n","")*
    "p: Process the data"*TUIcheck(ctrl,"process")*"\n"*
    "e: Export the results\n"*
    "l: Logs and templates\n"*
    "o: Options\n"*
    "u: Update\n"*
    "c: Clear\n"*
    "a: Extra\n"*
    "x: Exit\n"*
    "?: Help"
    return msg
end

function TUImethodSwitch(method::Nothing,ifG::String,ifC::String)
    return ifG
end
function TUImethodSwitch(method::Gmethod,ifG::String,ifC::String)
    return ifG
end
function TUImethodSwitch(method::Cmethod,ifG::String,ifC::String)
    return ifC
end

function TUIshowMethods(ctrl::AbstractDict)
    methods = _KJ["methods"].names
    msg = ""
    for i in eachindex(methods)
        msg *= string(i)*": "*methods[i]*"\n"
    end
    msg *= "c: concentrations\n"
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIdirfileMessage(ctrl::AbstractDict)
    msg = 
    "d: Read a directory in which analysis is stored in a different file\n" *
    "b: Set the number of blocks per analysis (current value = " * string(ctrl["nblocks"]) * ")\n" *
    "p: Parse the data from a single file using a laser log (provide paths)\n" *
    "P: Parse the data from a single file using a laser log (choose from list)\n" *
    "x: Exit\n" *
    "?: Help"
    return msg
end

function TUIinternalMessage(ctrl::AbstractDict)
    msg = "Choose an internal standard from the following list of channels:\n"
    msg *= TUIlistChannels(ctrl)
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIlistChannels(ctrl::AbstractDict)
    msg = ""
    channels = getChannels(ctrl["run"])
    for i in eachindex(channels)
        msg *= string(i)*". "*channels[i]*"\n"
    end
    return msg
end

function TUImineralMessage(ctrl::AbstractDict)
    msg = "Automatically set the concentration of the internal standard " *
        "by selecting one of the following minerals, or specify a value manually:\n"
    minerals = _KJ["stoichiometry"].names
    for i in eachindex(minerals)
        msg *= string(i)*". "*minerals[i]*"\n"
    end
    msg *= "m. manual"
    return msg
end

function TUIstoichiometryMessage(ctrl::AbstractDict)
    element = channel2element(ctrl["cache"])[1]
    msg = "Specify the concentration (in wt%) of " * element * " in the sample:"
    return msg
end

function TUIcolumnMessage(ctrl::AbstractDict)
    msg = "Choose from the following list of channels:\n"
    msg *= TUIlistChannels(ctrl)
    msg *= "and select the channels corresponding to "*
    "the following isotopes or their proxies:\n"
    P, D, d = getPDd(ctrl["method"].name)
    msg *= P *", "* D *", "* d *"\n"
    msg *= "Specify your selection as a "*
    "comma-separated list of numbers:"
    return msg
end

function TUIions2isotopes(ions::AbstractVector)
    out = String[]
    for ion in ions
        element = channel2element(ion)
        isotopes = _KJ["nuclides"][element]
        for isotope in isotopes
            push!(out,element*string(isotope))
        end
    end
    return unique(out)
end

function TUIlistIsotopes(cache::ProxySelectionCache)
    msg = ""
    nc = length(cache.channels)
    ni = length(cache.isotopes)
    n = max(nc,ni)
    for i in 1:n
        if i <= nc
            msg *= string(i) * ". " * cache.channels[i] * ", "
        else            
            msg *= repeat(" ",5 + maximum(length.(cache.channels)))
        end
        if i <= ni
            msg *= string(i) * ". " * cache.isotopes[i]
        end
        msg *= "\n"
    end
    return msg
end

function TUIsetProxiesMessage(ctrl::AbstractDict)
    msg = "KJ is having trouble mapping channels to isotopes. " *
          "Inspect the following two lists:\n"
    msg *= TUIlistIsotopes(ctrl["cache"])
    msg *= "Specify your selection as a " *
           "bracketed list of paired numbers, e.g. " *
           "(1,2),(3,1),(2,3)"
    return msg
end

function TUIsetInterferenceProxyMessage(ctrl::AbstractDict)
    msg = "Choose the isotope that matches channel " *
          ctrl["cache"].interference.channel * " :\n"
    isotopes = TUIions2isotopes([ctrl["cache"].key])
    for i in eachindex(isotopes)
        msg *= string(i) * ". " * isotopes[i] * "\n"
    end
    return msg
end

function TUIlistIsotopesMessage(ctrl::AbstractDict)
    msg = "Choose the isotope that requires an interference correction:\n"
    pairings = TUIgetPairings(ctrl)
    for i in eachindex(pairings)
        msg *= string(i) * ". " * pairings[i].proxy * "\n"
    end
    msg *= "x: Exit\n" * "?: Help"
    return msg
end

function TUIchooseInterferenceIonMessage(ctrl::AbstractDict)
    proxy = ctrl["cache"].target.proxy
    msg = "Choose the isotope that interferes with " *
          proxy * " from the following list:\n"
    interferences = TUIgetInterferences(proxy)
    for i in eachindex(interferences)
        msg *= string(i)*". "*interferences[i]*"\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIchooseInterferenceProxyChannelMessage(ctrl::AbstractDict)
    msg = "Choose an interference-free proxy channel for the " *
          ctrl["cache"].key * "-interference on " *
          ctrl["cache"].target.channel * " from the following list:\n"
    msg *= TUIlistChannels(ctrl)
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIchooseREEInterferenceProxyChannelMessage(ctrl::AbstractDict)
    msg = 
    "Suppose that X is a REE whose oxide (XO, say) interferes with " *
    ctrl["cache"].target.proxy *
    ", then the interference correction is given by X x YO / Y, where Y is " * 
    "a non-interfering REE with known isotopic abundance relative to X and YO is its oxide. " *
    "Select the channels corresponding to X, Y, and YO as a comma-separated list of numbers:\n"
    msg *= TUIlistChannels(ctrl)
    msg *= "x: Exit\n" * "?: Help"
    return msg
end

function TUIchooseStandardMessage(ctrl::AbstractDict)
    msg = "Choose one of the following reference materials:\n"
    refmats = TUIgetStandards(ctrl["method"])
    for i in eachindex(refmats.names)
        refmat = get(refmats,i)
        msg *= string(i)*": "*refmats.names[i]*TUIprintRefmatInfo(refmat)*"\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIchooseGlassMessage(ctrl::AbstractDict)
    msg = "Choose one of the following reference glasses:\n"
    refmats = _KJ["glass"].names
    for i in eachindex(refmats)
        msg *= string(i)*": " * refmats[i] * "\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
end 

function TUIprintRefmatInfo(refmat::AbstractRefmat)
    return " ("*refmat.material*")"
end

function TUIprintRefmatInfo(refmat::DataFrameRow)
    return ""
end

function TUIgetStandardsHelper(method::Gmethod;
                               condition::Function= refmat -> true)
    refmats = _KJ["refmat"][method.name]
    out = OrderedDict()
    for i in eachindex(refmats.names)
        refmat = get(refmats,i)
        if condition(refmat)
            add2od!(out,refmats.names[i],refmat)
        end
    end
    return out
end

function TUIgetStandards(method::Cmethod)
    return _KJ["glass"]
end

function TUIgetStandards(method::Gmethod)
    condition = refmat -> !(refmat isa BiasRefmat)
    return TUIgetStandardsHelper(method; condition = condition)
end

function TUIgetBiasStandards(method::Gmethod)
    condition = refmat -> refmat isa PointRefmat || refmat isa BiasRefmat
    return TUIgetStandardsHelper(method; condition = condition)
end

function TUIaddByPrefixMessage(ctrl::AbstractDict)
    msg = "Specify the prefix of the " * 
    get_refmat_from_cache(ctrl["cache"]) * 
    " measurements (? for help, x to exit):" 
    return msg 
end 

function TUIaddByNumberMessage(ctrl::AbstractDict)
    msg = "Select the " * 
        get_refmat_from_cache(ctrl["cache"]) * 
        " measurements as a comma-separated list of numbers " *
        "(? for help, x to exit):"
    return msg
end

function TUIratioMessage(ctrl::AbstractDict)
    channels = TUIgetChannels(ctrl)
    msg = "Choose one of the following denominators:\n"
    for i in eachindex(channels)
        msg *= string(i)*": "*channels[i]*"\n"
    end
    msg *= "or\n"
    msg *= "n: No denominator. Plot the raw signals\n"
    msg *= "?: Help"
    return msg
end

function TUIchooseBiasElementMessage(ctrl::AbstractDict) 
    m = ctrl["method"]
    elements = TUIgetBiasElements(m)
    msg = "Choose the element for which you want to fit a mass bias correction :\n"
    for i in eachindex(elements)
        msg *= string(i) * ". " * elements[i] * "\n"
    end
    msg *= "x: Exit\n" * "?: Help"
    return msg
end

function TUIgetBiasElements(m::Gmethod)
    elements = [channel2element(m.D.proxy)]
    for pairing in (m.P,m.D,m.d)
        for (key,interference) in pairing.interferences
            if interference isa Interference
                element = channel2element(key)
                push!(elements,element)
            end
        end
    end
    return elements
end

function TUIcalibrationMessage(ctrl::AbstractDict)
    element = ctrl["cache"].element
    isotopes = element .* string.(_KJ["nuclides"][element])
    channels = getChannels(ctrl["run"])
    ni = length(isotopes)
    nc = length(channels)
    n = max(ni,nc)
    msg = 
    "Pair the isotope with the channel for the " *
    "numerator and denominator of your bias correction:\n"
    for i in 1:n
        if i <= ni
            msg *= string(i) * ". " * isotopes[i] *", "
        else
            msg *= repeat(" ",5 + maximum(length.(isotopes)))
        end
        if i <= nc
            msg *= string(i) * ". " * channels[i]
        end
        msg *= "\n"
    end
    msg *= 
    "For example (1,1),(2,3) would mean that " * 
    isotopes[1] * " is paired with channel " * channels[1] * " and " *
    isotopes[2] * " is paired with channel " * channels[3] * "."
    return msg
end

function TUIchooseBiasStandardMessage(ctrl::AbstractDict)
    msg = "Choose the reference material for the bias correction from the following list:\n"
    refmats = TUIgetBiasStandards(ctrl["method"])
    for i in eachindex(refmats.names)
        refmat = get(refmats,i)
        msg *= string(i)*": " * refmats.names[i] * TUIprintRefmatInfo(refmat) * "\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIexportFormatMessage(ctrl::AbstractDict)
    return TUIexportFormatMessage(ctrl["method"])
end
function TUIexportFormatMessage(method::Gmethod)
    msg = "c: Export to .csv\n" 
    msg *= "j: Export to .json\n"
    msg *= "x: Exit\n"
    msg *= "?: Help"
    return msg
end
function TUIexportFormatMessage(method::Cmethod)
    msg = "c: Export to .csv\n" 
    msg *= "x: Exit\n"
    msg *= "?: Help"
    return msg
end

function TUIsetNblankMessage(ctrl::AbstractDict)
    msg = "Enter a non-negative integer (current value = " *
    string(ctrl["method"].nblank) * ", ? for help, x to exit):"
    return msg
end

function TUIsetNdriftMessage(ctrl::AbstractDict)
    msg = "Enter a non-negative integer (current value = " *
    string(ctrl["method"].ndrift) * ", ? for help, x to exit)"
    return msg
end

function TUIsetNdownMessage(ctrl::AbstractDict)
    msg = "Enter a non-negative integer (current value = " *
    string(ctrl["method"].ndown) * ", ? for help, x to exit)"
    return msg
end

function TUIchooseColumnMessage(ctrl::AbstractDict)
    colnames = names(ctrl["cache"])[1:end-2]
    msg = "Choose one of the following columns:\n"
    for i in eachindex(colnames)
        msg *= string(i) * ". " * colnames[i] * "\n"
    end
    return msg
end

function TUIchooseClimsMessage(ctrl::AbstractDict)
    channel = ctrl["mapcolumn"]
    colnames = names(ctrl["cache"])
    values = ctrl["cache"][:,channel]
    minval = minimum(values[values.>0])
    minrounded = round(minval;sigdigits=3)
    maxval = maximum(values)
    maxrounded = round(maxval;sigdigits=3)
    m = string(minrounded)
    M = string(maxrounded)
    msg = "Choose the limits of the colour scale as a comma-separated " *
        "list of two numbers. FYI: the values of " *
        colnames[channel] * " range from " * m * " to " * M *
        ". Enter 'r' to reset to the default values."
    return msg
end