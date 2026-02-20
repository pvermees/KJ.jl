# Command-line interface (REPL)

Advanced users can interact directly with KJ via the Julia REPL ("read-eval-print loop") for maximum flexibility and control over data processing pipelines.

## Loading Data

Load LA-ICP-MS data from a directory:

```julia
using KJ

# Load data in Agilent format
myrun = load("data/Lu-Hf"; format="Agilent")

# Load data in ThermoFisher format
myrun = load("data/iCap"; format="ThermoFisher")

# Load data with timestamps from separate file
myrun = load("data/timestamp/Moreira_data.csv",
             "data/timestamp/Moreira_timestamps.csv";
             format="Agilent")

# View a summary of the loaded data
summarise(myrun)
```

## Defining Analysis Methods

### Simple Method Definition

Create a method for standard isotopic ratio measurements:

```julia
# U-Pb dating method
method = Gmethod(name="U-Pb", groups=Dict("STDCZ" => "Plesovice"))

# Lu-Hf method with custom pairings
method = Gmethod(name="Lu-Hf",
                 groups=Dict("hogsbo" => "Hogsbo", "NIST612p" => "NIST612"),
                 P=Pairing(ion="Lu176", proxy="Lu175", channel="Lu175 -> 175"),
                 D=Pairing(ion="Hf176", channel="Hf176 -> 258"),
                 d=Pairing(ion="Hf177", proxy="Hf178", channel="Hf178 -> 260"),
                 standards=Set(["hogsbo"]))
```

### Advanced Method Configuration

Add calibrations and interference corrections:

```julia
# Add mass bias calibration
Calibration!(method; standards=Set(["NIST612p"]))

# Define interference for a particular ion
method.D.interferences["Lu176"] = Interference(proxy="Lu175", channel="Lu175 -> 257")

# Re-Os example with complex interferences
method = Gmethod(name="Re-Os",
                 groups=Dict("Nis3" => "NiS-3",
                             "Nist_massbias" => "NIST610",
                             "Nist_REEint" => "NIST610",
                             "Qmoly" => "QMolyHill"),
                 P=Pairing(ion="Re187", proxy="Re185", channel="Re185 -> 185"),
                 D=Pairing(ion="Os187", channel="Os187 -> 251"),
                 d=Pairing(ion="Os188", proxy="Os189", channel="Os189 -> 253"),
                 standards=Set(["Qmoly"]))

# Add REE interference correction
method.P.interferences["Tm169 -> 185"] = REEInterference(REE="Lu175 -> 191",
                                                        REEO="Ir191 -> 191",
                                                        standards=Set(["Nist_REEint"]))
```

### Concentration Method

For elemental concentration determinations:

```julia
# Define an internal standard-based method
method = Cmethod(myrun;
                 groups=Dict("NIST612p" => "NIST612"),
                 internal=("Al27 -> 27", 1.2e5))
```

## Setting Groups and Standards

Assign reference material group classifications to samples:

```julia
# Set groups based on filename patterns
groups = Dict("Nis3" => "NiS-3",
              "Nist_massbias" => "NIST610",
              "Qmoly" => "QMolyHill")
setGroup!(myrun, collect(keys(groups)))

# Or set directly via method
setGroup!(myrun, method)
```

## Processing Data

### Window Setting

Define background and signal windows:

```julia
# Set automatic windows based on data features
setBwin!(myrun)  # Background windows
setSwin!(myrun)  # Signal windows

# Manually override for specific sample
sample_index = 2
setSwin!(myrun[sample_index], [(70,90), (100,140)])
setBwin!(myrun[sample_index], [(0,22)]; seconds=true)

# Visualize windows
p = KJ.plot(myrun[sample_index];
            channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)
```

### Basic Processing

Process the entire dataset with a single command:

```julia
# Define method and process
method = Gmethod(name="U-Pb", groups=Dict("STDCZ" => "Plesovice", "91500" => "91500"))
fit = process!(myrun, method)

# Access drift and downhole fractionation parameters
println("Drift: ", fit.drift, ", Downhole: ", fit.down)
```

## Visualization

### Basic Plotting

View individual analyses:

```julia
# Plot a sample with specified channels
p = KJ.plot(myrun[1]; channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)

# Plot with log transformation
p = KJ.plot(myrun[1];
            channels=["Lu175 -> 175", "Hf176 -> 258", "Hf178 -> 260"],
            den="Hf178 -> 260",
            transformation="log")
display(p)
```

### Processing Results

Visualize fitted results:

```julia
# Plot a standard with fitted results
p = KJ.plot(myrun[2], method; fit=fit,
            transformation="log",
            den=method.D.channel)
display(p)

# Plot predicted isochron
p = KJ.plot(samp, method; fit=fit,
            den=method.D.channel,
            transformation="log",
            return_offset=true)
```

### Specialized Plots

Plot various data features:

```julia
# Plot internal isochron
p = internoplot(myrun[7], method, fit)
display(p)

# Plot elemental concentration map
p = plotMap(conc, "ppm[U] from U238"; clims=(0, 500))
display(p)
```

## Outlier Detection

Detect and flag anomalies:

```julia
# Detect outliers in specified channels
channels = ["K39 -> 39", "Ca40 -> 59", "Ca44 -> 63"]
detect_outliers!(myrun; include_samples=true, channels=channels)

# Visualize outliers
p = KJ.plot(myrun[1]; channels=channels, transformation="log")
display(p)
```

## Data Analysis and Statistics

### Extract Isotopic Ratios

Compute mean isotopic ratios for all samples:

```julia
# Average ratios across all measurements
ratios = averat(myrun, method, fit)

# Filter by sample prefix and export
selection = prefix2subset(ratios, "hogsbo")
```

### Atomic-scale Analysis

Compute spatial variations of isotopic ratios:

```julia
# Get position-dependent isotopic data
a = atomic(myrun[2], method, fit)

# Add x-y position information (for maps)
a = atomic(myrun[10], method, fit; add_xy=true)

# Convert to DataFrame for further analysis
df = DataFrame(a)
```

### Internal Isochron Calculation

Compute internal (within-grain) isochrons:

```julia
# Calculate internal isochron from line transect
isochron = internochron(myrun, method, fit)

# Save to CSV
CSV.write("output/isochron.csv", isochron)
```

## Exporting Results

### IsoplotR Format

Export results for use in IsoplotR:

```julia
# Export all samples
export2IsoplotR(myrun, method, fit; fname="output/U-Pb.json")

# Export with prefix filter
export2IsoplotR(myrun, method, fit;
                prefix="Duff",
                fname="output/Duff.json")

# Export pre-computed ratios
export2IsoplotR(ratios, method; fname="output/Lu-Hf.json")
```

### CSV Export

Save ratio tables to CSV:

```julia
using CSV

# Export averaged ratios
CSV.write("output/Lu-Hf.csv", selection)

# Export concentrations
conc = concentrations(myrun, method, fit)
CSV.write("output/concentrations.csv", conc)
```

For detailed function documentation, see the [API reference](api.md).
