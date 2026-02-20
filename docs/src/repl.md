# Examples (REPL)

Advanced users can interact directly with KJ via the Julia REPL ("read-eval-print loop") for maximum flexibility and control over data processing pipelines.

## Loading Data

Some labs prefer to store each analysis in a separate data file. Others prefer to put all analyses into a single file, which is parsed into chunks according to a separately provided laser log file. A third option is to divide each analysis into several files corresponding to background, signal and washout, respectively. Sample names are either stored in the data files, or have to be derived from the file names. `KJ` accommodates all these methods of operation. They are illustrated here with example data files that are packaged with the software.

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
- a mass bias correction for <sup>187</sup>Os/<sup>188</sup>Os measured in the NiS-3 reference material
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

To measure elemental concentration, `KJ` uses stoichiometric elements as internal references. It then compares all other elements with this internal standard, and with a concentration standard such as NIST glass:

```julia
map_method = Cmethod(map_run;
                     groups=Dict("NIST612" => "NIST612"),
                     internal=getInternal("zircon","Si29"))
```

## Processing Data

### Selection Windows

A typical laser ablation analysis consists of three phases:

1. an initial phase in which the laser fires into a shutter whilst the ICP-MS measures the *background*
2. opening of the shutter ("time zero") followed by the measurements of the *signal*
3. switching off the laser and *washout* of the signal from the tubing of the instrument

`KJ` aims to automatically partition these measurements into these three phases. However, the user can also manually override its decisions. To view the windows for the second sample in the Lu-Hf run:

```julia
sample_index = 2
p = KJ.plot(LuHf_run[sample_index];
            channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)
```

Now change the blank window for all the samples from 0 to 22 seconds:

```julia
setBwin!(LuHf_run, [(0,22)]; seconds=true)
```

Create a two-part signal window for the selected sample:

```julia
setSwin!(LuHf_run[sample_index], [(70,90), (100,140)])
```

Automatically set the background and signal windows for all samples:

```julia
setBwin!(LuHf_run)
setSwin!(LuHf_run)
```

### Basic Processing

Process the entire dataset with a single command. This automatically places samples into the correct groups, applies corrections for background, interferences and mass bias, before parameterising the elemental fractionation.

```julia
UPb_fit = process!(UPb_run, UPb_method)
```

## Visualisation

### Basic Plotting

Plot the first sample of the LuHf run with specified channels:
```julia
p = KJ.plot(LuHf_run[1]; channels=["Hf176 -> 258", "Hf178 -> 260"])
display(p)
```

Take ratios and plot on a log scale:
```julia
p = KJ.plot(LuHf_run[1];
            channels=["Lu175 -> 175", "Hf176 -> 258", "Hf178 -> 260"],
            den="Hf178 -> 260",
            transformation="log")
display(p)
```

### Processing Results

The output of `process!` can be added to the plot for the reference materials. This allows a visual inspection of the goodness-of-fit:
```julia
p = KJ.plot(UPb_run[2], UPb_method; fit=UPb_fit,
            transformation="log")
display(p)
```

### Specialised Plots

One of the unique features of `KJ` is its ability to fit *internal isochrons* based on time resolved atomic abundance estimates. This can be used to correct the non-radiogenic component of individual crystals without the need to assume an isochron intercept. Applying this idea to a laser ablation line across a large garnet:

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

`KJ` also provides very simple functionality to plot elemental concentration maps:
```julia
map_fit = process!(map_run,map_method)
conc = concentrations(map_run[10],map_method,map_fit)
p = plotMap(conc,"ppm[U] from U238";
            clims=(0,500),
            colorbar_scale=:identity)
display(p)
```

The output of the `concentrations(...)` function can also be saved to a data frame or output file for more sophisticated plotting in other Julia packages, or in external software:
```julia
using DataFrames
df = DataFrame(a)
```

## Outlier Detection

By default, `KJ` automatically identifies "spikes" in the time resolved data, using a multivariate outlier detection algorithm that fits a multivariate normal distribution to the raw data using the robust minimum covariance determinant (MCD) estimation method. Given this fitted distribution, outliers are flagged by applying "Tukey's fences" to the Mahalanobis distance. Applying this algorithm to some K-Ca data:

```julia
channels = ["K39 -> 39", "Ca40 -> 59", "Ca44 -> 63"]
detect_outliers!(KCa_run; include_samples=true, channels=channels)
p = KJ.plot(KCa_run[1]; channels=channels, transformation="log")
display(p)
```

## Data Analysis and Statistics

### Estimate atomic abundances

`KJ` converts the time resolved signals (which can be recorded in V, A or Hz) to time resolved estimates of atomic abundances in arbitrary units. This low level output can be extracted from `KJ` and processed in Julia or other software.
```julia
using Plots
a = atomic(LuHf_run[1],LuHf_method,LuHf_fit)
p = Plots.scatter(a.D,a.d,label="",xlabel="P",ylabel="d")
display(p)
```

### Extract Isotopic Ratios

Compute mean isotopic ratios for all samples:
```julia
using CSV
LuHf_fit = process!(LuHf_run, LuHf_method)
LuHf_ratios = averat(LuHf_run, LuHf_method, LuHf_fit)
selection = prefix2subset(LuHf_ratios, "hogsbo")
CSV.write("output/hogsbo.csv",selection)
```

### Internal Isochron Calculation

Save the ages and intercepts of all internal isochrons in `LuHf_run` to a `.csv` file:
```julia
isochron = internochron(LuHf_run, LuHf_method, LuHf_fit)
CSV.write("output/isochron.csv", isochron)
```

## Exporting to IsoplotR

Export the U-Pb results for use in [IsoplotR](https://isoplotr.es.ucl.ac.uk):

```julia
export2IsoplotR(UPb_run, UPb_method, UPb_fit; fname="output/U-Pb.json")
```

Select specific samples with a prefix filter:

```julia
export2IsoplotR(UPb_run, UPb_method, UPb_fit;
                prefix="GJ1",
                fname="output/GJ1.json")
```

Export pre-computed ratios:

```julia
export2IsoplotR(LuHf_ratios, LuHf_method; fname="output/Lu-Hf.json")
```

For detailed function documentation, see the [API reference](api.md).
