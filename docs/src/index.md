# KJ.jl

KJ ("Kasper Julia") is a Julia package for LA-ICP-MS (Laser Ablation Inductively Coupled Plasma Mass Spectrometry) data reduction.

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

See the [API documentation](api.md) for complete function reference.
