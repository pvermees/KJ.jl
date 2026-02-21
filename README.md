# KJ

## Julia package for LA-ICP-MS data reduction

📖 **Documentation:** <https://pvermees.github.io/KJ.jl/>

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
Pkg.precompile("KJ")
```

There are three ways to interact with KJ:

1. [TUI](#tui-text-based-user-interface): an interactive text-based user interface
2. [REPL](#repl-command-line-interface): the command-line interface
3. [Hybrid](#tui-repl): combining the TUI and REPL

## TUI (text-based user interface)

The `TUI()` function starts a menu-driven KJ session:

```julia
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

```

The workflow goes from the top to the bottom, with asterisks marking compulsory steps. To read ICP-MS data files, type `r` and proceed.

## REPL (command-line interface)

Advanced users can interact with KJ via Julia's command line interface
or REPL ("read-eval-print loop"). Here is an example of a carbonate
U-Pb data reduction using WC-1 for time-dependent elemental
fractionation correction between U and Pb and NIST-612 for
mass-dependent fractionation correction of the Pb-isotopes. The script
exports all the aliquots of the "Duff" sample to a JSON file that can
be opened in [IsoplotR](https://isoplotr.es.ucl.ac.uk):

```julia
myrun = load("data/carbonate", format="Agilent")
method = Gmethod(name="U-Pb", groups=Dict("WC1"=>"WC1"))
fit = process!(myrun, method)
export2IsoplotR(myrun, method, fit;
                prefix="Duff", fname="output/Duff.json")
```

Type `?load`, `?process!`, `?export2IsoplotR` or `?KJ` at the REPL to
view the documentation.

## TUI + REPL

You can seamlessly switch from the TUI to the REPL and back. The
following example stores the TUI settings into a variable called
`ctrl` (type `ctrl` at the REPL to view its contents). You can
manipulate the contents of `ctrl` and sync it with the TUI using the
`setKJctrl!()` function.

```julia
julia> using KJ
julia> TUI()
# ... use TUI ...
# press 'x' to exit

julia> ctrl = getKJctrl();
julia> KJ.plot(ctrl["run"][1])
```

## Building the Documentation Locally

```
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The generated site will be placed in `docs/build/`.
