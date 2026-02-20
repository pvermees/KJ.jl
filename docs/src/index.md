# KJ.jl

KJ ("Kasper Julia") is a Julia package for **LA-ICP-MS data reduction** (Laser Ablation Inductively Coupled Plasma Mass Spectrometry).

## Installation

First, ensure [Julia](https://julialang.org/downloads/#current_stable_release) is installed. Then at the Julia REPL:

```julia
import Pkg; Pkg.add(url="https://github.com/pvermees/KJ.jl.git")
```

Optionally, precompile the package for faster startup:

```julia
Pkg.build("KJ")
```

## Ways to Interact with KJ

### 1. TUI (Text-Based User Interface)

KJ provides an interactive text-based user interface for data processing workflows.

```julia
using KJ
TUI()
```

This launches an interactive menu where you can:
- Load LA-ICP-MS data files
- Specify analysis methods
- View and adjust samples
- Set interferences and fractionation corrections
- Process data and export results

The TUI provides options for reading data files, specifying methods, adjusting samples, handling interferences, applying fractionation and mass bias corrections, processing data, and exporting results.

### 2. REPL (Command-Line Interface)

Advanced users can use Julia's command-line interface for direct scripting and analysis control. Here is a simple example of a U-Pb data processing session using test data that is packaged with `KJ`:

```julia
using KJ

# Load LA-ICP-MS data
myrun = load("data/U-Pb"; format="Agilent", head2name=false)

# Define the analysis method
method = Gmethod(name="U-Pb", groups=Dict("STDCZ" => "Plesovice"))

# Process the data
fit = process!(myrun, method)

# Export results to IsoplotR format
export2IsoplotR(myrun, method, fit; fname="output/U-Pb.json")
```

Type `?load`, `?process!`, `?export2IsoplotR`, or `?KJ` at the REPL for detailed documentation.

### 3. Hybrid: TUI + REPL

You can seamlessly switch between the TUI and REPL. From the TUI, press `x` twice to exit and return to the Julia REPL:

```julia
julia> TUI()
# ... use TUI ...
# press 'x' to exit

julia> ctrl = getKJctrl()
julia> samp = ctrl["run"][1]
julia> plot(samp)
```

You can also synchronize changes made in the REPL back to the TUI:

```julia
julia> setKJctrl!(ctrl)
julia> TUI()  # Resume with updated settings
```
Seamlessly switch between the TUI and REPL. Store TUI settings into a variable and manipulate them programmatically using `getKJctrl()` and `setKJctrl()` functions.

See the [API documentation](api.md) for complete function reference.
