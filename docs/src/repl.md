# Examples (REPL)

Advanced users can interact directly with KJ via the Julia REPL ("read-eval-print loop") for maximum flexibility and control over data processing pipelines.

## Loading Data

Some labs prefer to store each analysis in a separate data file. Others prefer to put all analyses into a single file, which is parsed into chunks according to a separately provided laser log file. A third option is to divide each analysis into several file corresponding to background, signal and washout, respectively. Sample names are either stored in the data files, or have to be derived from the file names. `KJ` accomodates all these methods of operation. They are illustrated here with example data files that are packaged with the software.

Load LA-ICP-MS data from a directory of `.csv` files from an Agilent instrument:
```julia
using KJ
LuHf_run = load("data/Lu-Hf"; format="Agilent")
ReOs_run = load("data/Re-Os"; format="Agilent")
KCa_run = load("data/K-Ca"; format="Agilent")
```

Load data in a ThermoFisher format:
```julia
iCap_run = load("data/iCap"; format="ThermoFisher")
```

Load UPb data with sample names extracted from file names:
```julia
UPb_run = load("data/U-Pb"; format="Agilent", head2name=false)
```

Load data from a Neptune MC-ICP-MS with background, signal and washout split in different files:
```julia
MC_run = load("data/FIN2"; format="FIN2", nblocks=3)
```

Load data from one massive `.csv` file with timestamps provided in a separate laser log file:
```julia
map_run = load("data/timestamp/NHM_data.csv",
               "data/timestamp/NHM_timestamps.csv";
               format="Agilent")
```

## View a summary of the loaded data

Show all the samples in a run:
```julia
summarise(UPb_run)
```

## Defining Analysis Methods

### Simple Geochronology

`KJ` implements two types of methods `Gmethod` for geochronology, and `Cmethod` for chemical concentration measurements. The simplest definition of a `Gmethod` requires just the name of the geochronometer of interest:

```julia
UPb_method = Gmethod(name="U-Pb")
```

Additional arguments can be used to specify the prefixes of reference materials as `groups`. Prefixes that are missing from this Dict will be treated as samples. For example, suppose that a run uses two primary reference materials: Plesovice and 91500; which are labelled as `STDCZ` and `91500`, respectively:

```julia
UPb_method = Gmethod(name="U-Pb", 
                     groups=Dict("STDCZ" => "Plesovice", "91500" => "91500"),
                     standards=Set(["STDCZ","91500"]))
```

The `.csv` files with raw data contain column headers ("channels") that correspond to different mass-to-charge ratios in the mass spectrometer. For the U-Pb dataset shown above, it is easy to match these column headers with the geochronologically relevant isotopes ("U238", "Pb206" and "Pb207"). This is not always the case. In the case of &beta;-chronometry (Lu-Hf, Re-Os, K-Ca, Rb-Sr), the column headers include chemical information that is not so easy for `KJ` to disentangle. These methods are based on the following age equation:

```math
\frac{D}{d} = \left[\frac{D}{d}\right]_0 + \frac{P}{d} \left(e^{\lambda t} - 1\right)
```

where `P`, `D` and `d` refer to the radioactive "parent", radiogenic "daughter" and non-radiogenic "sister" isotope. An additional complexity arises because `P` and `d` are often measured by proxy, using a different isotope. These complexities are accommodated by providing the `Gmethod` constructor with `Pairing`s:

```julia
LuHf_method = Gmethod(name="Lu-Hf",
                      groups=Dict("hogsbo" => "Hogsbo", "NIST612p" => "NIST612"),
                      P=Pairing(ion="Lu176", proxy="Lu175", channel="Lu175 -> 175"),
                      D=Pairing(ion="Hf176", channel="Hf176 -> 258"),
                      d=Pairing(ion="Hf177", proxy="Hf178", channel="Hf178 -> 260"),
                      standards=Set(["hogsbo"]))
```

### Advanced Geochronology

Add a mass bias correction for <sup>176</sup>Hf/<sup>178</sup>Hf based on the known ratio of NIST-612 standard glass, and a correction for the isobaric interference of <sup>176</sup>Lu on <sup>176</sup>Hf, measured by proxy using the reaction products of <sup>175</sup>Lu and the known <sup>176</sup>Lu/<sup>175</sup>Lu ratio of the Solar System:

```julia
Calibration!(LuHf_method; standards=Set(["NIST612p"]))
LuHf_method.D.interferences["Lu176"] = Interference(proxy="Lu175", channel="Lu175 -> 257")
```

A Re-Os example with complex interferences, including:

- an elemental fractionation correction using QMolyHill molybdenite as a primary reference material
- a mass bias correction for <sup>187</sup>Os<sup>188</sup>Os measured in the NiS-3 reference material
- an interference correction of <sup>187</sup>Re on <sup>187</sup>Os measured using NIST-610 glass
- a Rare Earth interference of TmO on <sup>185</sup>Re, using LuO/Lu ratio as a proxy

```julia
ReOs_method = Gmethod(name="Re-Os",
                      groups=Dict("Nis3" => "NiS-3",
                                  "Nist_massbias" => "NIST610",
                                  "Nist_REEint" => "NIST610",
                                  "Qmoly" => "QMolyHill"),
                      P=Pairing(ion="Re187", proxy="Re185", channel="Re185 -> 185"),
                      D=Pairing(ion="Os187", channel="Os187 -> 251"),
                      d=Pairing(ion="Os188", proxy="Os189", channel="Os189 -> 253"),
                      standards=Set(["Qmoly"]))
Calibration!(ReOs_method;standards=Set(["Nis3"]))
Re_bias = Calibration(num=(ion="Re187",channel="Os187 -> 251"),
                      den=(ion="Re185",channel="Re185 -> 249"),
                      standards=Set(["Nist_massbias"]))
ReOs_method.D.interferences["Re187"] = Interference(proxy="Re185",
                                                    channel="Re185 -> 249",
                                                    bias=Re_bias)
ReOs_method.P.interferences["Tm169 -> 185"] = REEInterference(REE="Lu175 -> 191",
                                                              REEO="Ir191 -> 191",
                                                              standards=Set(["Nist_REEint"]))
```

### Concentrations

To measure elemental concentration, `KJ` uses stoichiometric elements as internal references. It then compares all other elements with this inernal standard, and with a concentration standard such as NIST glass:

```julia
map_method = Cmethod(map_run;
                     groups=Dict("NIST612" => "NIST612"),
                     internal=getInternal("zircon","Si29"))
```

## Processing Data

### Window Setting

Define background and signal windows:

```julia
# Set automatic windows based on data features
setBwin!(LuHf_run)  # Background windows
setSwin!(LuHf_run)  # Signal windows

# Manually override for specific sample
sample_index = 2
setSwin!(LuHf_run[sample_index], [(70,90), (100,140)])
setBwin!(LuHf_run[sample_index], [(0,22)]; seconds=true)

# Visualise windows
p = KJ.plot(LuHf_run[sample_index];
            channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)
```

### Basic Processing

Process the entire dataset with a single command:

```julia
UPb_fit = process!(UPb_run, UPb_method)

# Access drift and downhole fractionation parameters
println("Drift: ", UPb_fit.drift, ", Downhole: ", UPb_fit.down)
```

## Visualisation

### Basic Plotting

View individual analyses:

```julia
# Plot the first sample of the LuHf run with specified channels
p = KJ.plot(LuHf_run[1]; channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)

# Plot with log transformation
p = KJ.plot(LuHf_run[1];
            channels=["Lu175 -> 175", "Hf176 -> 258", "Hf178 -> 260"],
            den="Hf178 -> 260",
            transformation="log")
display(p)
```

### Processing Results

Visualise fitted results:

```julia
# Plot a standard with fitted results
p = KJ.plot(UPb_run[2], UPb_method; fit=UPb_fit,
            transformation="log")
display(p)
```

### Specialized Plots

Plot an internal isochron:

```julia
lines_run = load("data/lines",format="Agilent")
lines_method = Gmethod(name="Lu-Hf",
                       groups=Dict("Hog" => "Hogsbo"),
                       P=Pairing(ion="Lu176",proxy="Lu175",channel="Lu175 -> 175"),
                       D=Pairing(ion="Hf176",proxy="Hf176",channel="Hf176 -> 258"),
                       d=Pairing(ion="Hf177",proxy="Hf178",channel="Hf178 -> 260"),
                       standards=Set(["Hog"]),
                       ndown=0,
                       ndrift=1)
lines_fit = process!(lines_run,lines_method)
p = internoplot(lines_run[7],lines_method,lines_fit)
display(p)
```

Plot an elemental concentration map

```julia
map_fit = process!(map_run,map_method)
conc = concentrations(map_run[10],map_method,map_fit)
p = plotMap(conc,"ppm[U] from U238";
            clims=(0,500),
            colorbar_scale=:identity)
display(p)
```

## Outlier Detection

Detect and flag anomalies:

```julia
# Detect outliers in specified channels
channels = ["K39 -> 39", "Ca40 -> 59", "Ca44 -> 63"]
detect_outliers!(KCa_run; include_samples=true, channels=channels)

# Visualise outliers
p = KJ.plot(KCa_run[1]; channels=channels, transformation="log")
display(p)
```

## Data Analysis and Statistics

### Estimate atomic abundances

Compute spatial variations of isotopic ratios:

```julia
# Reload the earlier concentration data for geochronological purposes:
UPb_map_run = load("data/timestamp/NHM_data.csv",
                   "data/timestamp/NHM_timestamps.csv";
                   format="Agilent")
UPb_map_method = Gmethod(name="U-Pb", groups=Dict("91500"=>"91500"))
UPb_map_fit = process!(UPb_map_run,UPb_map_method)

# Estimate the atomic abundances along an x,y-grid
a = atomic(UPb_map_run[10],UPb_map_method,UPb_map_fit;add_xy=true)

# Convert to a DataFrame for further analysis
using DataFrames
df = DataFrame(a)
```

### Extract Isotopic Ratios

Compute mean isotopic ratios for all samples:

```julia
# Average ratios across all measurements
LuHf_ratios = averat(LuHf_run, LuHf_method, LuHf_fit)

# Filter by sample prefix and export
selection = prefix2subset(LuHf_ratios, "hogsbo")
```

### Internal Isochron Calculation

Compute internal (within-grain) isochrons:

```julia
# Calculate internal isochrons from line transects
LuHf_fit = process!(LuHf_run, LuHf_method)
isochron = internochron(LuHf_run, LuHf_method, LuHf_fit)

# Save the ages and intercepts to CSV
using CSV
CSV.write("output/isochron.csv", isochron)
```

## Exporting Results

### IsoplotR Format

Export results for use in [IsoplotR](https://isoplotr.es.ucl.ac.uk):

```julia
# Export all samples
export2IsoplotR(UPb_run, UPb_method, UPb_fit; fname="output/U-Pb.json")

# Export with prefix filter
export2IsoplotR(UPb_run, UPb_method, UPb_fit;
                prefix="GJ1",
                fname="output/GJ1.json")

# Export pre-computed ratios
export2IsoplotR(LuHf_ratios, LuHf_method; fname="output/Lu-Hf.json")
```

### CSV Export

Save ratio tables to CSV:

```julia
using CSV

# Export averaged ratios
CSV.write("output/Lu-Hf.csv", selection)

# Export concentrations
conc = concentrations(map_run, map_method, map_fit)
CSV.write("output/concentrations.csv", conc)
```

For detailed function documentation, see the [API reference](api.md).
