function PT(debug=false)
    welcome()
    control = Dict(
        "priority" => ["load","standards","process","instrument","read"],
        "history" => DataFrame(task=String[],action=String[]),
        "chain" => ["top"]
    )
    while true
        key = control["chain"][end]
        (message,action) = tree(key,control)
        println(message)
        response = readline()
        next = action[response]
        if next=="x"
            pop!(control["chain"])
        elseif isa(next,Tuple)
            next[1](next[2],control)
        elseif isa(next,AbstractString)
            push!(control["chain"],next)
        end
        push!(control["history"],[key,response])
        if debug
            println(control["history"])
            println(control["chain"])
        end
        if length(control["chain"])==0 return end
    end
end
export PT

function tree(key::AbstractString,ctrl::AbstractDict)
    branches = Dict(
        "top" => (
            message =
            "r: Read the data files"*check(ctrl["priority"],"load")*"\n"*
            "s: Mark the standards"*check(ctrl["priority"],"standards")*"\n"*
            "b: Bulk settings\n"*
            "v: View and adjust each sample\n"*
            "p: Process the data"*check(ctrl["priority"],"process")*"\n"*
            "e: Export the results\n"*
            "l: Import/export a session log\n"*
            "x: Exit",
            action = Dict(
                "r" => (dispatch!,"load"),
                "s" => "standards",
                "b" => "bulk",
                "v" => "view",
                "p" => "process",
                "e" => "export",
                "l" => "log",
                "x" => "x"
            )
        )
    )
    return branches[key]
end

function welcome()
    version = string(pkgversion(@__MODULE__))
    title = " Plasmatrace "*version*" \n"
    width = Base.length(title)-1
    println('-'^width*"\n"*title*'-'^width)
end

function check(pl::AbstractVector,action::AbstractString)
    return action in pl ? "[*]" : ""
end

function dispatch!(task::AbstractString,ctrl::AbstractDict)
    println("TODO: "*task)
    pop!(ctrl["chain"])
end
