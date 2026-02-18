# KJ

## Julia package for LA-ICP-MS data reduction

KJ ("Kasper Julia") is still in development and has not yet been added
to the [Julia](https://julialang.org/) package repository. However, if
you want to play around with the current functionality, then you can
install the package from GitHub. First, make sure that you have Julia
installed on your system by downloading it from
[here](https://julialang.org/downloads/#current_stable_release). Then,
at the Julia REPL:

```
import Pkg; Pkg.add(url="https://github.com/pvermees/KJ.jl.git")
```

To save yourself some time when you run `KJ` later on, it is useful
(but not strictly necessary) to *precompile* the package now:

```
Pkg.build("KJ")
```

There are four ways to interact with KJ:

1. [TUI](#tui-text-based-user-interface): an interactive text-based user interface
2. [REPL](#repl-command-line-interface): the command-line interface
3. [Hybrid](#tui-repl): combining the TUI and REPL

## TUI (text-based user interface)

Here is a short example of a menu-driven KJ session:

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

<img src="./img/plot.png" width="480px">

## REPL (command-line interface)

Advanced users can interact with KJ via Julia's command line interface
or REPL ("read-eval-print loop"). Here is an example of a carbonate
U-Pb data reduction using WC-1 for time-dependent elemental
fractionation correction between U and Pb and NIST-612 for
mass-dependent fractionation correction of the Pb-isotopes. The script
exports all the aliquots of the "Duff" sample to a JSON file that can
be opened in IsoplotR:

```
myrun = load("data/U-Pb";format="Agilent",head2name=false)
method = Gmethod(name="U-Pb",groups=Dict("STDCZ" => "Plesovice"))
fit = process!(myrun,method)
export2IsoplotR(myrun,method,fit;
                fname="output/U-Pb.json")
```

Type `?load`, `?process!`, `?export2IsoplotR` or `?KJ` at the REPL to
view the documentation.

## TUI + REPL

You can seamlessly switch from the TUI to the REPL and back. The
following example stores the TUI settings into a variable called
`ctrl` (type `ctrl` at the REPL to view its contents). You can
manipulate the contents of `ctrl` and sync it with the TUI using the
`setKJctrl()` function.

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
x

julia> ctrl = getKJctrl();
julia> KJ.plot(ctrl["run"][1])
```
