function TUIwelcomeMessage(ctrl::AbstractDict)
    msg = "r: Read data files"*TUIcheck(ctrl,"load")*"\n"*
    "m: Specify the method"*TUIcheck(ctrl,"method")*"\n"*
    "t: Tabulate the samples\n"*
    "v: View and adjust each sample\n"*
    "f: Fractionation"*TUIcheck(ctrl,"standards")*"\n"*
    "b: Mass bias\n"*
    "i: Interferences\n"*
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

function TUIlistIsotopes(ctrl::AbstractDict)
    msg = ""
    ions = getIons(ctrl["method"])
    nuclidelist = ions2nuclidelist(ions)
    for i in eachindex(nuclidelist)
        msg *= string(i) * ". " * nuclidelist[i] * "\n"
    end
    return msg
end

function ions2nuclidelist(ions::NamedTuple)
    out = []
    all_elements = string.(keys(_KJ["nuclides"]))
    for ion in ions
        matching_elements = filter(x -> occursin(x, ion), all_elements)
        for matching_element in matching_elements
            all_isotopes = string.(_KJ["nuclides"][matching_element])
            matching_isotopes = filter(x -> occursin(x, ion), all_isotopes)
            if length(matching_isotopes) > 0
                for isotope in all_isotopes
                    nuclide = matching_element * isotope
                    push!(out,nuclide)
                end
            end
        end
    end
    return unique(out)
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

function TUIsetProxiesMessage(ctrl::AbstractDict)
    msg = "KJ is having trouble mapping channels to isotopes. " *
          "Choose from the following list of isotopes:\n"
    msg *= TUIlistIsotopes(ctrl)
    msg *= "and select those corresponding to "*
    "the channels that you selected earlier:\n"
    channels = getChannels(ctrl["method"];as_tuple=true)
    msg *= channels.P *", "* channels.D *", "* channels.d *"\n"
    msg *= "Specify your selection as a "*
    "comma-separated list of numbers:"
    return msg
end

function TUIchooseStandardMessage(ctrl::AbstractDict)
    msg = "Choose one of the following reference materials:\n"
    refmats = TUIgetRefmats(ctrl["method"])
    for i in eachindex(refmats.names)
        refmat = get(refmats,i)
        msg *= string(i)*": "*refmats.names[i]*TUIprintRefmatInfo(refmat)*"\n"
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

function TUIgetRefmats(method::Gmethod)
    return _KJ["refmat"][method.name]
end

function TUIgetRefmats(method::Cmethod)
    return _KJ["glass"]
end

function TUIchooseGlassMessage(ctrl::AbstractDict)
    error("TODO")
end

function TUIaddByPrefixMessage(ctrl::AbstractDict)
    msg = "Specify the prefix of the " * ctrl["cache"] *
        " measurements (? for help, x to exit):"
    return msg
end

function TUIaddByNumberMessage(ctrl::AbstractDict)
    msg = "Select the " * ctrl["cache"] * " measurements " *
        "with a comma-separated list of numbers " *
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

function TUIexportFormatMessage(ctrl::AbstractDict)
    msg = "c: Export to .csv\n" 
    if ctrl["method"]!="concentrations"
        msg *= "j: Export to .json\n"
    end
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