# KJ.jl

KJ ("Kasper Julia") is a Julia package for **LA-ICP-MS data reduction** (Laser Ablation Inductively Coupled Plasma Mass Spectrometry). The package is still in development and currently available via GitHub.

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
An interactive menu-driven interface for performing analyses:

```
julia> using KJ
julia> TUI()
```

The TUI provides options for reading data files, specifying methods, adjusting samples, handling interferences, applying fractionation and mass bias corrections, processing data, and exporting results.

### 2. REPL (Command-Line Interface)
Advanced users can use Julia's command-line interface for direct scripting and analysis control (see Quick Start below).

### 3. Hybrid: TUI + REPL
Seamlessly switch between the TUI and REPL. Store TUI settings into a variable and manipulate them programmatically using `getKJctrl()` and `setKJctrl()` functions.

## Quick Start

Load and process data using the REPL:

```julia
using KJ

# Load LA-ICP-MS data
myrun = load("data/U-Pb"; format="Agilent")

# Define the analysis method
method = Gmethod(name="U-Pb", groups=Dict("STDCZ" => "Plesovice"))

# Process the data
fit = process!(myrun, method)

# Export results to IsoplotR format
export2IsoplotR(myrun, method, fit; fname="output/U-Pb.json")
```

Type `?load`, `?process!`, `?export2IsoplotR`, or `?KJ` at the REPL for detailed documentation.

See the [API documentation](api.md) for complete function reference.
