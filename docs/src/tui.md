# Text-Based User Interface (TUI)

The TUI provides an interactive, menu-driven workflow for loading data, defining methods, processing analyses, and exporting results.

## Start the TUI

```julia
using KJ
TUI()
```

## Example Session (Menu-Driven)

```
julia> using KJ
julia> TUI()
----------
 KJ 0.8.2
----------

r: Read data files[*]
m: Specify the method[*]
t: Tabulate the samples
v: View and adjust each sample
i: Interferences
f: Fractionation[*]
b: Mass bias
p: Process the data[*]
e: Export the results
l: Logs and templates
o: Options
u: Update
c: Clear
a: Extra
x: Exit
?: Help
r

a: csv (Agilent)
t: csv (ThermoFisher)
f: FIN2
x: Exit
?: Help
a

d: Read a directory in which analysis is stored in a different file
b: Set the number of blocks per analysis (current value = 1)
p: Parse the data from a single file using a laser log (provide paths)
P: Parse the data from a single file using a laser log (choose from list)
x: Exit
?: Help
d

Enter the full path of the data directory (? for help, x to exit):
data/Lu-Hf

r: Read data files
m: Specify the method[*]
t: Tabulate the samples
v: View and adjust each sample
i: Interferences
f: Fractionation[*]
b: Mass bias
p: Process the data[*]
e: Export the results
l: Logs and templates
o: Options
u: Update
c: Clear
a: Extra
x: Exit
?: Help
v
```

![Example plot](assets/plot.png)

## Logs and Templates

The TUI includes a built-in system for saving session logs and reusable templates. These are accessed from the main menu via `l: Logs and templates`. A `KJ` log is a text file that records all the key strokes of a TUI session. For example, the log for the previous section looks like this:
```
task,action
top,r
format,a
dir|file,d
loadICPdir,data/Lu-Hf
top,v
view,x
```

### Export a Session Log

Use logs to save a processing session and replay it later:

1. From the main menu, press `l`.
2. Choose `e: Export the session log`.
3. Enter a file path, e.g. `logs/Lu-Hf.log`.

You can then replay the session with:

```julia
TUI(logbook="logs/Lu-Hf.log")
```

Alternatively, you can also continue a previous session from within the TUI:

1. Press `l`.
2. Choose `i: Import a session log`.
3. Provide the path to the saved log file.

### Save a Template

Templates store default settings (format, method definition, transformations, and other defaults) in a Julia script so you can avoid re-entering them.

1. Press `l`.
2. Choose `s: Save a template`.
3. Provide a path such as `templates/Lu-Hf.tmp`.

Example template file:

```julia
format = "Agilent"
multifile = true
head2name = true
transformation = "log"
nblocks = 1
method = Gmethod(name="Lu-Hf",
                 P=Pairing(ion="Lu176", proxy="Lu175", channel="Lu175 -> 175"),
                 D=Pairing(ion="Hf176", proxy="Hf176", channel="Hf176 -> 258"),
                 d=Pairing(ion="Hf177", proxy="Hf178", channel="Hf178 -> 260"),
                 nblank=2,
                 ndrift=2,
                 ndown=1,
                 PAcutoff=1.0e7)
method.groups = Dict("hogsbo" => "Hogsbo", "NIST612p" => "NIST612")
method.standards = Set(["hogsbo"])
```

### Open a Template

To apply a template during a session:

1. Press `l`.
2. Choose `o: Open a template`.
3. Provide the path to the template file.

## Extensions

The TUI uses a decision tree that is saved as a Dict. Extensions are Julia packages that modify this Dict to alter the default behaviour of the TUI. Extensions can be used to add or remove branches of TUI tree, or to swap out functions. `KJgui` is a proof-of-concept example of a simple extension, which adds GUI elements to the TUI, including interactive file choosers. This extension can be obtained from [https://github.com/pvermees/KJgui.jl](https://github.com/pvermees/KJgui.jl) and used as follows:
```julia
julia> using KJ, KJgui
julia> TUI(KJgui)
```

For detailed function documentation, see the [API reference](api.md).
