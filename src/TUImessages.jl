function TUIwelcomeMessage(ctrl::AbstractDict)
    msg = "r: Read data files"*TUIcheck(ctrl,"load")*"\n"*
    "m: Specify the method"*TUIcheck(ctrl,"method")*"\n"*
    "t: Tabulate the samples\n"*
    "s: Mark mineral standards"*TUIcheck(ctrl,"standards")*"\n"*
    "g: Mark reference glasses"*TUIcheck(ctrl,"glass")*"\n"*
    "v: View and adjust each sample\n"*
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
    for i in eachindex(ctrl["channels"])
        msg *= string(i)*". "*ctrl["channels"][i]*"\n"
    end
    msg *= "x: Exit\n"*"?: Help"
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
    labels = names(getSignals(ctrl["run"][1]))
    for i in eachindex(labels)
        msg *= string(i)*". "*labels[i]*"\n"
    end
    msg *= "and select the channels corresponding to "*
    "the following isotopes or their proxies:\n"
    P, D, d = getPDd(ctrl["method"])
    msg *= P *", "* D *", "* d *"\n"
    msg *= "Specify your selection as a "*
    "comma-separated list of numbers:"
    return msg
end

function TUIchooseStandardMessage(ctrl::AbstractDict)
    msg = "Choose one of the following standards:\n"
    standards = _KJ["refmat"][ctrl["method"]].names
    for i in eachindex(standards)
        msg *= string(i)*": "*standards[i]*"\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
end

function TUIchooseGlassMessage(ctrl::AbstractDict)
    msg = "Choose one of the following reference glasses:\n"
    glasses = _KJ["glass"].names
    for i in eachindex(glasses)
        msg *= string(i)*": "*glasses[i]*"\n"
    end
    msg *= "x: Exit\n"*"?: Help"
    return msg
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
    if isa(ctrl["channels"],AbstractVector)
        channels = ctrl["channels"]
    elseif isa(ctrl["channels"],AbstractDict)
        channels = collect(values(ctrl["channels"]))
    else
        channels = getChannels(ctrl["run"])
    end
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
    string(ctrl["options"]["blank"]) * ", ? for help, x to exit):"
    return msg
end

function TUIsetNdriftMessage(ctrl::AbstractDict)
    msg = "Enter a non-negative integer (current value = " *
    string(ctrl["options"]["drift"]) * ", ? for help, x to exit)"
    return msg
end

function TUIsetNdownMessage(ctrl::AbstractDict)
    msg = "Enter a non-negative integer (current value = " *
    string(ctrl["options"]["down"]) * ", ? for help, x to exit)"
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
