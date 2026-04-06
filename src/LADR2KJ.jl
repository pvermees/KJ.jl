"""
    LADR2KJ(ifile::AbstractString; ofile::AbstractString=joinpath(pwd(), "windows.log"))

Convert a LADR time-window CSV file to a KJ signal-window log.

The input table is expected to contain `Start`, `End`, and `Join` columns.
Rows with `Join == 1` are merged with the preceding row and written as a
single `oneMultiSignalWindow` entry; all other rows are written as
`oneSingleSignalWindow` entries. The resulting command sequence is saved to
`ofile` and can be replayed in the KJ TUI.

# Arguments
- `ifile`: Path to the LADR CSV file.
- `ofile`: Path to the output KJ window log file.
"""
function LADR2KJ(ifile::AbstractString;
                 ofile::AbstractString=joinpath(pwd(), "windows.log"))
    
    df = CSV.read(ifile, DataFrame)
    Join = findall(df[:,"Join"] .== 1)
    multi = sort(unique([Join; Join.-1]))
    out = ""
    nr = nrow(df)
    rnum = 1
    snum = 1
    while rnum <= nr
        out *= "view,g\n" *
               "goto," * string(snum) * "\n" *
               "view,s\n"
        if rnum in multi
            out *= "Swin,m\n" *
                   "oneMultiSignalWindow,\""
            windows = String[]
            while rnum in multi
                window = "(" * string(df[rnum,"Start"]) * 
                         "," * string(df[rnum,"End"]) * ")"
                push!(windows, window)
                rnum += 1
            end
            out *= join(windows, ",") * "\"\n"
        else
            out *= "oneSingleSignalWindow,\"" * 
                   string(df[rnum,"Start"]) * "," * 
                   string(df[rnum,"End"]) * "\"\n"
            rnum += 1
        end
        snum += 1
    end
    write(ofile, out)
end
export LADR2KJ