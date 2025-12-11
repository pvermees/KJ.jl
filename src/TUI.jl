"""
KJ(extensions...;logbook::AbstractString="",reset=false)

KJ TUI

# Arguments

- `extensions...`: a comma separated list of KJ extensions
- `logbook`: the file path to a KJ log
- `reset`: hard reset of all inherited settings

# Examples
```julia
TUI(logbook="logs/test.log")
ctrl = getKJctrl()
samp = ctrl["run"][1]
p = plot(samp)
display(p)
```
"""
function TUI(extensions...;logbook::AbstractString="",reset=false)
    if reset
        init_KJ!()
    end
    _KJ["extensions"] = extensions
    TUIwelcome()
    if isnothing(_KJ["ctrl"])
        TUIinit!()
    else
        _KJ["ctrl"]["chain"] = ["top"]
    end
    for (i, extension) in enumerate(_KJ["extensions"])
        extension.extend!(_KJ)
    end
    if logbook != ""
        _KJ["tree"]["importLog"].action(_KJ["ctrl"],logbook)
    end
    while true
        if length(_KJ["ctrl"]["chain"])<1 return end
        try
            dispatch!(_KJ["ctrl"];verbose=false)
        catch e
            println(e)
            println("Press 'x' to quit or any other key to continue")
            if readline() == "x" break end
        end
    end
end
export TUI

function dispatch!(ctrl::AbstractDict;
                   key = ctrl["chain"][end],
                   response = nothing,
                   verbose::Bool=false)
    if verbose
        println(ctrl["chain"])
        print("key: " * key)
        println(", response: " * response)
    end
    (message,help,action) = _KJ["tree"][key]
    if isa(message,Function)
        println("\n" * message(ctrl))
    else
        println("\n" * message)
    end
    if isnothing(response)
        response = readline()
    end
    if response == "?"
        println(help)
        next = nothing
    elseif response in ["x","xx","xxx","xxxx"]
        next = response
    elseif isa(action,Function)
        next = action(ctrl,response)
    else
        next = action[response]
    end
    if isa(next,Function)
        final = next(ctrl)
    else
        final = next
    end
    if final == "x"
        if length(ctrl["chain"])<1 return end
        pop!(ctrl["chain"])
    elseif final == "xx"
        if length(ctrl["chain"])<2 return end
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
    elseif final == "xxx"
        if length(ctrl["chain"])<3 return end
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
    elseif final == "xxxx"
        if length(ctrl["chain"])<4 return end
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
        pop!(ctrl["chain"])
    elseif isnothing(final)
        if length(ctrl["chain"])<1 return end
    else
        push!(ctrl["chain"],final)
    end
    toskip = [TUItabulate,TUIclear!,"saveTemplate","exportLog","importLog"]
    if !(key in toskip) & !(next in toskip) & !(final in toskip)
        push!(ctrl["history"],[key,response])
    end
    if verbose
        print("key: " * key)
        print(", response: " * response)
        print(", next: "); print(next)
        print(", final: "); println(final)
        println(ctrl["history"])
    end
end

"""
KJtree!(tree::AbstractDict)

Replace the internal structure of TUI tree, including all the KJ settings.
"""
function KJtree!(tree::AbstractDict)
    _KJ["tree"] = tree
end
export KJtree!
"""
KJtree()

Returns the internal structure of TUI tree, including all the KJ settings.
"""
function KJtree()
    return _KJ["tree"]
end
export KJtree

function init_KJtree()
    Dict(
        "top" => (
            message = TUIwelcomeMessage,
            help = "This is the top-level menu. Asterisks (*) " * 
            "mark compulsory steps. Refresh reloads the data directory " * 
            "and only works if the asterisks are gone.",
            action = Dict(
                "r" => TUIread,
                "m" => "method",
                "t" => TUItabulate,
                "v" => TUIviewer,
                "f" => "standards",
                "b" => TUItodo!, # "glass",
                "i" => TUItodo!, # "interferences",
                "p" => TUIprocess!,
                "e" => "export",
                "l" => "log",
                "o" => "options",
                "a" => "extra",
                "u" => TUIrefresh!,
                "c" => TUIclear!
            )
        ),
        "format" => (
            message =
            "a: csv (Agilent)\n" * 
            "t: csv (ThermoFisher)\n" *
            "f: FIN2\n" * 
            "x: Exit\n" * 
            "?: Help",
            help = "Choose a file format. Email us if you don't " * 
            "find your format in this list.",
            action = TUIformat!
        ),
        "dir|file" => (
            message =
            "d: Read a directory of individual data files\n" * 
            "p: Parse the data from a single file using a laser log (provide paths)\n" *
            "P: Parse the data from a single file using a laser log (choose from list)\n" *
            "(? for help, x to exit):",
            help =
            "There are two ways to store mass spectrometer data. Each analysis " * 
            "(spot or raster) can be saved into its own file, or all analyses can be " * 
            "stored together in a single file. In the latter case, a laser log file is " * 
            "needed to parse the data into individual samples. Option p requires " *
            "the full paths of the data and log files. Option P requires the name of the " *
            "directory in which the data and laser log are stored.",
            action = Dict(
                "d" => "loadICPdir",
                "p" => "loadICPfile",
                "P" => "choosedir"
            )
        ),
        "loadICPdir" => (
            message =
            "Enter the full path of the data directory " *
            "(? for help, x to exit):",
            help =
            "KJ will read all the files in this folder. " * 
            "You don't need to select the files, just the folder.",
            action = TUIloadICPdir!
        ),
        "loadICPfile" => (
            message =
            "Enter the full path of the ICP-MS data file " * 
            "(? for help, x to exit):",
            help =
            "KJ will read the pooled data in this file and " * 
            "parse it using a laser log file, which is to be provided next.",
            action = TUIloadICPfile!
        ),
        "loadLAfile" => (
            message =
            "Enter the full path of the laser log file " * 
            "(? for help, x to exit):",
            help =
            "KJ will parse the ICP-MS data provided in the " * 
            "previous step using this log file.",
            action = TUIloadLAfile!
        ),
        "choosedir" => (
            message =
            "Specify the directory containing the data and laser log file " * 
            "(? for help, x to exit):",
            help =
            "Provide the full path of the directory that contains the " *
            "LA-ICP-MS data and the laser log file. In the next step, " *
            "KJ will list the content of this directory and ask you to " *
            "pick two files from it.",
            action = TUIchoosedir!
        ),
        "pickICPLAfiles" => (
            message =
            "Select the ICP-MS data and laser log file as " *
            "a comma-separated list of two numbers " *
            "(? for help, x to exit):",
            help =
            "For example, enter 3,5 to select the third and fifth " *
            "data file in the directory as the ICP-MS file and " *
            "laser log, respectively.",
            action = TUIpickICPLAfiles!
        ),
        "method" => (
            message = TUIshowMethods,
            help = "Choose a method. Email us if you cannot " * 
            "find the one you are looking for.",
            action = TUImethod!
        ),
        "internal" => (
            message = TUIinternalMessage,
            help = nothing,
            action = TUIinternal!
        ),
        "columns" => (
            message = TUIcolumnMessage,
            help = nothing,
            action = TUIsetChannels!
        ),
        "setProxies" => (
            message = TUIsetProxiesMessage,
            help = "KJ is having trouble mapping the mass spectrometer " *
            "channels to isotopic proxies. Please help it by identifying " *
            "the isotopes measured by the different channels by listing " *
            "them as a comma-separated list of numbers.",
            action = TUIsetProxies!
        ),
        "mineral" => (
            message = TUImineralMessage,
            help = nothing,
            action = TUIchooseMineral!
        ),
        "stoichiometry" => (
            message = TUIstoichiometryMessage,
            help = nothing,
            action = TUIstoichiometry!
        ),
        "standards" => (
            message =
            "a: Add a reference material\n" * 
            "r: Remove reference materials\n" * 
            "l: List the available reference materials\n" * 
            "t: Tabulate all the samples\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Choose one or more primary reference materials. " * 
            "Note that secondary reference materials should be " * 
            "treated as regular samples.",
            action = Dict(
                "a" => "chooseStandard",
                "r" => "removeStandard",
                "l" => TUIrefmatTab,
                "t" => TUItabulate
            )
        ),
        "chooseStandard" => (
            message = TUIchooseStandardMessage,
            help =
            "If you do not find your reference material in this list, " * 
            "then you can either specify your own reference " * 
            "material under 'options' in the top menu, or " * 
            "you can email us to add the material to the software.",
            action = TUIchooseStandard!
        ),
        "addStandardGroup" => (
            message =
            "p: Select samples by prefix\n" * 
            "n: Select samples by number\n" * 
            "t: Tabulate all the samples\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Select samples to label as primary standards, either " * 
            "by specifying the first few characters of their name, " * 
            "or by providing a comma-separated list of numbers. " * 
            "Enter 't' to remind yourself of the names and " * 
            "order of the samples.",
            action = Dict(
                "p" => "addStandardsByPrefix",
                "n" => "addStandardsByNumber",
                "t" => TUItabulate
            )
        ),
        "addStandardsByPrefix" => (
            message = TUIaddByPrefixMessage,
            help =
            "For example, suppose that Plesovice zircon reference " * 
            "materials are named STDCZ01, STDCZ02, ..., then you can " * 
            "select all the standards by entering STDCZ here. " * 
            "Enter 'x' to go up one level and tabulate the " * 
            "sample if you forgot the exact prefix of your standards.",
            action = TUIaddStandardsByPrefix!
        ),
        "addStandardsByNumber" => (
            message = TUIaddByNumberMessage,
            help =
            "For example, suppose that the analyses are labelled as " * 
            "G001, G002, ..., then it is not possible to identify " * 
            "the standards by prefix, but you can still select them " * 
            "by sequence number (e.g., 1,2,8,9,15,16,...).",
            action = TUIaddStandardsByNumber!
        ),
        "removeStandard" => (
            message =
            "a: Remove all standards\n" * 
            "s: Remove selected standards\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Remove or reset the samples that are marked as primary standards.",
            action = Dict(
                "a" => TUIremoveAllStandards!,
                "s" => "removeStandardsByNumber"
            )
        ),
        "removeStandardsByNumber" => (
            message =
            "Select the reference materials as a comma-separated list of numbers " * 
            "(? for help, x to exit):",
            help =
            "Select the samples that are to be removed as primary " * 
            "standards by providing a list of comma separated numbers " * 
            "(e.g., 1,2,8,9,15,16,...).",
            action = TUIremoveStandardsByNumber!
        ),
        "glass" => (
            message =
            "a: Add a reference glass\n" * 
            "r: Remove reference glasses\n" * 
            "l: List the available reference glasses\n" * 
            "t: Tabulate all the samples\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Choose one or more reference glasses. These are used " * 
            "to determine the mass dependent isotope fractionation factor.",
            action = Dict(
                "a" => "chooseGlass",
                "r" => "removeGlass",
                "l" => TUIglassTab,
                "t" => TUItabulate
            )
        ),
        "chooseGlass" => (
            message = TUIchooseGlassMessage,
            help =
            "If you do not find your reference glass in this list, " * 
            "then you can either specify your own under 'options' in " * 
            "the top menu, or you can email us to add the material to the software.",
            action = TUIchooseGlass!
        ),
        "addGlassGroup" => (
            message =
            "p: Select analyses by prefix\n" * 
            "n: Select analyses by number\n" * 
            "t: Tabulate all the analyses\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Select samples to label as reference glasses, either " * 
            "by specifying the first few characters of their name, " * 
            "or by providing a comma-separated list of numbers. " * 
            "Enter 't' to remind yourself of the names and " * 
            "order of the analyses.",
            action = Dict(
                "p" => "addGlassByPrefix",
                "n" => "addGlassByNumber",
                "t" => TUItabulate
            )
        ),
        "addGlassByPrefix" => (
            message = TUIaddByPrefixMessage,
            help =
            "For example, suppose that NIST-612 reference glasses " * 
            "are named GLASS01, GLASS02, ..., then you can " * 
            "select all the standards by entering GLASS here. " * 
            "Enter 'x' to go up one level and tabulate the " * 
            "sample if you forgot the exact prefix of your glasses." * 
            "Note that NIST glass is only used to determine chemical " * 
            "concentrations and mass fractionation factors, but is " * 
            "NOT used as an isotopic ratio standard. Conversely, " * 
            "minerals are only used as isotopic ratio standards " * 
            "but not as concentration standards. Use both to get the " * 
            "best results.",
            action = TUIaddGlassByPrefix!
        ),
        "addGlassByNumber" => (
            message = TUIaddByPrefixMessage,
            help =
            "For example, suppose that the analyses are labelled as " * 
            "G001, G002, ..., then it is not possible to identify " * 
            "the glasses by prefix, but you can still select them " * 
            "by sequence number (e.g., 1,2,8,9,15,16,...).",
            action = TUIaddGlassByNumber!
        ),
        "removeGlass" => (
            message =
            "a: Remove all glass measurements\n" * 
            "s: Remove selected analyses\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Unlabel the analyses that are marked as reference glass.",
            action = Dict(
                "a" => TUIremoveAllGlass!,
                "s" => "removeGlassByNumber"
            )
        ),
        "removeGlassByNumber" => (
            message =
            "Select the glass analyses as a comma-separated list of numbers " * 
            "(? for help, x to exit):",
            help =
            "Select the samples that are to be removed as reference " * 
            "glass by providing a list of comma separated numbers " * 
            "(e.g., 1,2,8,9,15,16,...).",
            action = TUIremoveGlassByNumber!
        ),
        "view" => (
            message = 
            "n: Next\n" * 
            "p: Previous\n" * 
            "g: Go to\n" * 
            "t: Tabulate all the samples in the session\n" * 
            "r: Plot signals or ratios?\n" *
            "0: Set 'time zero'\n" *
            "b: Select blank window(s)\n" * 
            "s: Select signal window(s)\n" *
            "m: Move selection windows\n" * 
            "d: Choose a data transformation\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "It is useful to view all the samples in your analytical " * 
            "sequence at least once, to modify blank and signal " * 
            "windows or checking the fit to the standards.",
            action = Dict(
                "n" => TUInext!,
                "p" => TUIprevious!,
                "g" => "goto",
                "t" => TUItabulate,
                "r" => "setDen",
                "0" => "t0",
                "b" => "Bwin",
                "s" => "Swin",
                "m" => "moveWin",
                "d" => "transformation"
            )
        ),
        "goto" => (
            message = "Enter the number of the sample to plot " * 
            "(? for help, x to exit):",
            help = "Jump to a specific analysis.",
            action = TUIgoto!
        ),
        "setDen" => (
            message = TUIratioMessage,
            help =
            "Plot the ratios of the channels relative to a common " * 
            "denominator. Zero values for the denominator are omitted.",
            action = TUIratios!
        ),
        "t0" => (
            message =
            "Choose an option to set the start of the signal ('time zero'):\n" * 
            "a: Automatic (current sample)\n" * 
            "s: Manually set t0 (current sample)\n" * 
            "A: Automatic (all samples)\n" * 
            "S: Manually set t0 (all samples)\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Specify t0 ('time zero'), marking the onset of laser ablation " * 
            "(in seconds) or trust KJ to choose t0 automatically.",
            action = Dict(
                "a" => TUIt0AutoOne!,
                "s" => "t0one",
                "A" => TUIt0AutoAll!,
                "S" => "t0all"
            )
        ),
        "t0one" => (
            message =
            "Enter the onset of laser ablation (t0) in seconds. " *
            "Type '?' for help and 'x' to exit.",
            help =
            "Specify the start of the signal (which equals the " *
            "the end of the blank) as a number of seconds.",
            action = TUIt0One!
        ),
        "t0all" => (
            message =
            "Enter the onset of laser ablation (t0) in seconds. " *
            "Type '?' for help and 'x' to exit.",
            help =
            "Specify the start of the signal (which equals the " *
            "the end of the blank) as a number of seconds.",
            action = TUIt0All!
        ),
        "Bwin" => (
            message =
            "Choose an option to set the blank window(s):\n" * 
            "a: Automatic (current sample)\n" * 
            "s: Manually set a one-part window (current sample)\n" * 
            "m: Manually set a multi-part window (current sample)\n" * 
            "A: Automatic (all samples)\n" * 
            "S: Manually set a one-part window (all samples)\n" * 
            "M: Manually set a multi-part window (all samples)\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Specify the blank (background signal) as pairs of time " * 
            "stamps (in seconds) or trust KJ to choose the " * 
            "blank windows automatically. The blanks of all the " * 
            "analyses (samples + standards) will be combined and " * 
            "interpolated under the signal.",
            action = Dict(
                "a" => TUIoneAutoBlankWindow!,
                "s" => "oneSingleBlankWindow",
                "m" => "oneMultiBlankWindow",
                "A" => TUIallAutoBlankWindows!,
                "S" => "allSingleBlankWindows",
                "M" => "allMultiBlankWindows"
            )
        ),
        "oneSingleBlankWindow" => (
            message =
            "Enter the start and end point of the selection window " * 
            "(in seconds). Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end point of the blank window " * 
            "as a comma-separated pair of numbers. For example: " * 
            "0,20 marks a blank window from 0 to 20 seconds.",
            action = TUIoneSingleBlankWindow!
        ),
        "oneMultiBlankWindow" => (
            message =
            "Enter the start and end points of the multi-part " * 
            "selection window (in seconds)\n" * 
            "Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end points of the blank window " * 
            "as a comma-separated list of bracketed pairs " * 
            " of numbers. For example: (0,20),(25,30) marks a two-part " * 
            "selection window from 0 to 20s, and from 25 to 30s.",
            action = TUIoneMultiBlankWindow!
        ),
        "allSingleBlankWindows" => (
            message =
            "Enter the start and end point of the selection window " * 
            "(in seconds). Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end point of the blank windows " * 
            "as a comma-separated pair of numbers. For example: " * 
            "0,20 marks a blank window from 0 to 20 seconds.",
            action = TUIallSingleBlankWindows!
        ),
        "allMultiBlankWindows" => (
            message =
            "Enter the start and end points of the multi-part " * 
            "selection window (in seconds). " * 
            "Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end points of the blank windows " * 
            "as a comma-separated list of bracketed pairs " * 
            " of numbers. For example: (0,20),(25,30) marks a two-part " * 
            "selection window from 0 to 20s, and from 25 to 30s.",
            action = TUIallMultiBlankWindows!
        ),
        "Swin" => (
            message =
            "Choose an option to set the signal window(s):\n" * 
            "a: Automatic (current sample)\n" * 
            "s: Manually set a one-part window (current sample)\n" * 
            "m: Manually set a multi-part window (current sample)\n" * 
            "A: Automatic (all samples)\n" * 
            "S: Manually set a one-part window (all samples)\n" * 
            "M: Manually set a multi-part window (all samples)\n" * 
            "x: Exit\n" * 
            "?: Help",
            help = 
            "Specify the signal as pairs of time stamps (in seconds) " * 
            "or trust KJ to choose the signal windows " * 
            "automatically. The signals of the reference materials " * 
            "are used to define the drift and fractionation " * 
            "corrections, which are then applied to the samples.",
            action = Dict(
                "a" => TUIoneAutoSignalWindow!,
                "s" => "oneSingleSignalWindow",
                "m" => "oneMultiSignalWindow",
                "A" => TUIallAutoSignalWindows!,
                "S" => "allSingleSignalWindows",
                "M" => "allMultiSignalWindows"
            )
        ),
        "oneSingleSignalWindow" => (
            message =
            "Enter the start and end point of the selection window " * 
            "(in seconds). Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end points of the signal window " * 
            "as a comma-separated pair of numbers. For example: " * 
            "30,60 marks a signal window from 30 to 60 seconds.",
            action = TUIoneSingleSignalWindow!
        ),
        "oneMultiSignalWindow" => (
            message =
            "Enter the start and end points of the multi-part " * 
            "selection window (in seconds). " * 
            "Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end points of the signal window " * 
            "as a comma-separated list of bracketed pairs " * 
            "of numbers. For example: (40,45),(50,60) marks a two-part " * 
            "signal window from 40 to 45s, and from 50 to 60s.",
            action = TUIoneMultiSignalWindow!
        ),
        "allSingleSignalWindows" => (
            message =
            "Enter the start and end point of the selection windows " * 
            "(in seconds). Type '?' for help and 'x' to exit.",
            help =
            "Specify the start and end points of the signal windows " * 
            "as a comma-separated pair of numbers. For example: " * 
            "30,60 marks a signal window from 30 to 60 seconds.",
            action = TUIallSingleSignalWindows!
        ),
        "allMultiSignalWindows" => (
            message =
            "Enter the start and end points of the multi-part " * 
            "selection windows (in seconds)\n" * 
            "Type '?' for help.",
            help =
            "Specify the start and end points of the signal windows " * 
            "as a comma-separated list of bracketed pairs " * 
            "of numbers. For example: (40,45),(50,60) marks a two-part " * 
            "signal window from 40 to 45s, and from 50 to 60s.",
            action = TUIallMultiSignalWindows!
        ),
        "moveWin" => (
            message =
            "Enter the shift in seconds:\n" * 
            "(? for help, x to exit):",
            help =
            "Moves all the time zero values and selection windows " *
            "by a fixed number of seconds (can either be positive or negative). " *
            "This can be useful to correct a misaligned laser log.",
            action = TUImoveWin!
        ),
        "transformation" => (
            message =
            "Choose a data transformation for the y-axis:\n" * 
            "l: Linear\n" * 
            "L: Logarithmic\n" * 
            "s: Square root\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Using a square root or log-transform makes it easier to compare " * 
            "signal strengths or signal ratios that vary over several " * 
            "orders of magnitude.",
            action = TUItransformation!
        ),
        "export" => (
            message = TUIexportFormatMessage,
            help =
            "Export to a comma-separated variable format with " * 
            "samples names as row names, for futher processing " * 
            "in Excel, say; or export to a .json format that " * 
            "can be opened in an online or offline instance " * 
            "of IsoplotRgui.",
            action = Dict(
                "c" => "csv",
                "j" => "json"
            )
        ),
        "csv" => (
            message = "Enter the path and name of the .csv " * 
            "file (? for help, x to exit):",
            help = "Provide the file name with or without " * 
            "the .csv extension.",
            action = TUIexport2csv
        ),
        "json" => (
            message = "Enter the path and name of the " * 
            ".json file (? for help, x to exit):",
            help = "Provide the file name with or without " * 
            "the .json extension.",
            action = TUIexport2json
        ),
        "log" => (
            message =
            "Choose an option:\n" * 
            "i: Import a session log\n" * 
            "e: Export the session log\n" * 
            "o: Open a template\n" * 
            "s: Save a template\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Session logs are an easy and useful way to " * 
            "save intermediate results, increase throughput, " * 
            "develop FAIR data processing chains, and report bugs.\n" * 
            "Methods are Julia scripts that save default settings " * 
            "whilst omitting sample-specific details.",
            action = Dict(
                "i" => "importLog",
                "e" => "exportLog",
                "o" => "openTemplate",
                "s" => "saveTemplate"
            )
        ),
        "importLog" => (
            message = "Enter the path and name of the log " * 
            "file (? for help, x to exit):",
            help = "Open a previous log and continue where you left off.",
            action = TUIimportLog!
        ),
        "exportLog" => (
            message = "Enter the path and name of the log " * 
            "file (? for help, x to exit):",
            help = "Save the current KJ session so that you " * 
            "can replicate your results later",
            action = TUIexportLog
        ),
        "openTemplate" => (
            message = "Enter the path and name of the template " * 
            "file (? for help, x to exit):",
            help = "Open default settings in a template file to avoid " * 
            "repetitive entry of the method details.",
            action = TUIopenTemplate!
        ),
        "saveTemplate" => (
            message = "Enter the path and name of the template " * 
            "file (? for help, x to exit):",
            help = "Save the current KJ method for use " * 
            "in a future session.",
            action = TUIsaveTemplate
        ),
        "options" => (
            message =
            "b: Set the polynomial order of the blank correction\n" * 
            "d: Set the polynomial order of the drift correction\n" * 
            "h: Set the polynomial order of the down hole fractionation correction\n" * 
            "p: Subset the data by P/A cutoff\n" * 
            "l: List the available age standards\n" * 
            "a: Define new age standards\n" * 
            "r: List the available reference glasses\n" * 
            "g: Define new reference glasses\n" * 
            "n: Configure the file name reader\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "Advanced settings to override the default behaviour " * 
            "of KJ. See the individual options for " * 
            "further information.",
            action = Dict(
                "b" => "setNblank",
                "d" => "setNdrift",
                "h" => "setNdown",
                "p" => "PA",
                "l" => TUIrefmatTab,
                "a" => "addStandard",
                "r" => TUIglassTab,
                "g" => "addGlass",
                "n" => "head2name"
            )
        ),
        "setNblank" => (
            message = TUIsetNblankMessage,
            help =
            "The blank is fitted by the following equation: " * 
            "b = a[1] + a[2]*t[1] + ... + a[n]*t^(n-1). " *
            "Here you can specify the value of n.",
            action = TUIsetNblank!
        ),
        "setNdrift" => (
            message = TUIsetNdriftMessage,
            help =
            "The session drift is fitted by the following equation: " * 
            "d = exp( a[1] + a[2]*t + ... + a[n]*t^(n-1) ). " * 
            "Here you can specify the value of n.",
            action = TUIsetNdrift!
        ),
        "setNdown" => (
            message = TUIsetNdownMessage,
            help =
            "The down-hole drift is fitted by the following equation: " * 
            "d = exp( a[1]*t + a[2]*t^2 + ... + a[n]*t^n ). " * 
            "Here you can specify the value of n.",
            action = TUIsetNdown!
        ),
        "PA" => (
            message =
            "l: List the maximum signal strength for each file\n" * 
            "a: Add a pulse/analog cutoff\n" * 
            "r: Remove the pulse/analog cutoff\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "On single collector ICP-MS instruments, low intensity " * 
            "ion beams are measured in 'pulse' (P) mode and high " * 
            "intensity ion beams are measured in 'analog' (A) mode. " * 
            "The intercalibration of these two modes is not always " * 
            "perfect. Here you can set or remove the cutoff value " * 
            "between P and A mode, so that KJ can use the " * 
            "standards to estimate a correction factor between the " * 
            "two modes.",
            action = Dict(
                "l" => TUIPAlist,
                "a" => "setPAcutoff",
                "r" => TUIclearPAcutoff!
            )
        ),
        "setPAcutoff" => (
            message =
            "Enter the cutoff between pulse and analog mode in cps " * 
            "(? for help, x to exit):",
            help =
            "After entering this value, samples and standards will be " * 
            "split into two groups, corresponding to pulse mode (maximum " * 
            "signal below the cutoff) and analog mode (maximum signal " * 
            "above the cutoff).",
            action = TUIsetPAcutoff!
        ),
        "addStandard" => (
            message =
            "Enter the path to the standards file (? for help, x to exit)",
            help =
            "Add new isotopic reference materials to KJ " * 
            "by specifying the path to a .csv file that is " * 
            "formatted as follows:\n\n" * 
            "method,name,t,st,y0,sy0\n" * 
            "Lu-Hf,Hogsbo,1029,1.7,3.55,0.05\n" * 
            "Lu-Hf,BP,1745,5.0,3.55,0.05\n\n" * 
            "where 't' is the age of the sample, 'y0' is the " * 
            "y-intercept of the inverse isochron, and 'st' and 'sy0' " * 
            "are their standard errors, which are not used for the " * 
            "calculations in this version of the sofware.",
            action = TUIaddStandard!
        ),
        "addGlass" => (
            message =
            "Enter the path to the glass file (? for help, x to exit)",
            help =
            "Add new concentration standards to KJ " * 
            "by specifying the path to a .csv file that is " * 
            "formatted as follows:\n\n" * 
            "SRM,Be,B,F,Mg,...,Bi,Th,U,Hf177Hf176,Rb87Rb86,Pb207Pb206\n" * 
            "NIST612,476,350,304,432,...,384,457.2,461.5,3.544566,1.410312,0.9097\n" * 
            "NIST610,37.5,34.3,80,68,...,30.2,37.79,37.38,3.544842,1.409344,0.9074",
            action = TUIaddGlass!
        ),
        "head2name" => (
            message =
            "h: Extract the sample names from the file headers\n" * 
            "f: Extract the sample names from the file names\n" * 
            "x: Exit\n" * 
            "?: Help",
            help =
            "In most cases, the name of each ablation spot is " * 
            "registered in the headers of the input files. However, " * 
            "occasionally, the headers contain generic names (e.g., " * 
            "spot01, spot02, ...) and the spot name must be inferred " * 
            "from the file name (e.g. /path/to/data/Plesovice-01.csv)",
            action = TUIhead2name!
        ),
        "extra" => (
            message =
            "i: Internal isochrons\n" *
            "t: Time-resolved estimates\n" *
            "m: Maps\n" *
            "x: Exit\n" * 
            "?: Help",
            help =
            "If you supplied any extensions when starting KJ, " *
            "they will tend to put their options in this menu",
            action = Dict(
                "i" => "internochron",
                "t" => "timeresolved",
                "m" => TUImapper
            )
        ),
        "internochron" => (
            message =
            "p: Plot\n" *
            "e: Export\n" *
            "x: Exit\n" *
            "?: Help",
            help =
            "Plot the internal isochron for a single laser spot " *
            "or export all the spots to a .csv file.",
            action = Dict(
                "p" => TUInternochronViewer!,
                "e" => "internocsv"
            )
        ),
        "internoview" => (
            message =
            "n: Next\n" *
            "p: Previous\n" *
            "g: Go to\n" *
            "t: Tabulate all the samples in the session\n" *
            "x: Exit\n" *
            "?: Help",
            help =
            "It is useful to view the outcome of the internal isochron " *
            "regression to ensure that the results are sensible and that the " *
            "algorithm did not get stuck in a local minimum.",
            action = Dict(
                "n" => TUInternochron_next!,
                "p" => TUInternochron_previous!,
                "g" => "internogoto",
                "t" => TUItabulate
            )
        ),
        "internogoto" => (
            message =
            "Enter the number of the sample to plot " *
            "(? for help, x to exit):",
            help = "Jump to a specific analysis.",
            action = TUInternochron_goto!
        ),
        "internocsv" => (
            message =
            "Enter the path and name of the .csv " *
            "file (? for help, x to exit):",
            help =
            "Provide the file name with or without " *
            "the .csv extension.",
            action = internochron2csv
        ),
        "timeresolved" => (
            message = "Enter the full path of the directory into " *
            "to save the csv files with the time-resolved results " *
            "(? for help, x to exit):",
            help =
            "KJ will export the all time-resolved results in a single " *
            "folder as a sequence of csv files.",
            action = TUItimeresolved2csv
        ),
        "map" => (
            message =
            "c: choose a column to plot\n" *
            "p: previous\n" *
            "n: next\n" *
            "g: goto\n" *
            "s: set the limits of the colour scale\n" *
            "x: Exit\n" *
            "?: Help",
            help =
            "Plots time resolved data in x,y-space. Only works " *
            "for data that come with a laser log.",
            action = Dict(
                "c" => "chooseMapColumn",
                "p" => TUImapPrevious!,
                "n" => TUImapNext!,
                "s" => "clims"
            )
        ),
        "chooseMapColumn" => (
            message = TUIchooseColumnMessage,
            help = "Choose a channel to form a colour scale for the map.",
            action = TUIchooseMapColumn!
        ),
        "clims" => (
            message = TUIchooseClimsMessage,
            help = "Adjust the colour scale by marking its minimum and maxmium value.",
            action = TUIchooseClims!
        )
    )
end
