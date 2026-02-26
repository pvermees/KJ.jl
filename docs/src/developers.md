# Developers

This section provides further details about the internal structure of `KJ`.

## Types

### `Sample`
Represents a single LA-ICP-MS (laser ablation inductively coupled plasma mass spectrometry) measurement. Contains the sample name, acquisition date/time, time-series measurement data, ablation onset time, and signal/blank windows.

### `KJmethod`
Abstract base type for LA-ICP-MS data reduction methods. Concrete subtypes include:

- **`Gmethod`**: Geochronology method for parent-daughter isotope dating systems (U-Pb, Rb-Sr, Lu-Hf, K-Ca, Re-Os). Specifies group definitions, isotope pairings (P, D, d), calibration, standards, and polynomial orders for blank fitting, drift correction, downhole fractionation, and bias correction.
- **`Cmethod`**: Concentration method for quantitative element analysis. Specifies elements, group definitions, internal standard (channel and concentration), and blank polynomial order.

### `Pairing`
Isotope pairing specification that links ions, proxies, and channels. Each pairing can have associated interference corrections. Fields include ion notation (e.g., "Pb206"), a proxy isotope for calculations, the measured channel name, and a dictionary of interference corrections.

### `AbstractInterference`
Abstract base type for interference correction specifications. Concrete subtypes include:

- **`Interference`**: General isobaric interference correction with a proxy isotope, corrected channel, and associated bias correction settings
- **`REEInterference`**: Rare earth element oxide interference correction with REE channel, REE oxide channel, and reference material standards

### `Calibration`
Settings for mass bias correction. Specifies numerator and denominator isotopes (with ion and channel names) and the set of reference materials used for calibration.

### `AbstractRefmat`
Abstract base type for reference material definitions. Concrete subtypes include:

- **`IsochronRefmat`**: Reference material with isochron-based age information (material name, age, uncertainty, intercept, intercept uncertainty)
- **`PointRefmat`**: Reference material with point-based estimates (material name, x ratio, x uncertainty, y ratio, y uncertainty)
- **`BiasRefmat`**: Reference material for bias correction (material name, intercept, intercept uncertainty)

### `AbstractAnchor`
Abstract base type for anchor point definitions in isotope ratio space. Anchor points represent known isotopic compositions of reference materials used for fractionation correction. Concrete subtypes include:

- **`IsochronAnchor`**: Isochron-based anchor with x₀, y₀, and y₁ coordinates
- **`PointAnchor`**: Point-based anchor with x and y coordinates
- **`BiasAnchor`**: Simple anchor with a single y coordinate (used for bias anchoring)

### `AbstractBias`
Abstract base type for mass bias correction parameters. Concrete subtypes include:

- **`Bias`**: Standard isotope ratio bias with numerator mass, denominator mass, and polynomial correction coefficients
- **`REEBias`**: Rare earth element oxide bias with polynomial correction coefficients

### `KJfit`
Abstract base type for fitted LA-ICP-MS correction parameters. Concrete subtypes include:

- **`Gfit`**: Fitted parameters for geochronology methods with blank correction coefficients, drift correction, downhole fractionation, analog mode drift, covariance matrix, and mass bias corrections by element
- **`Cfit`**: Fitted parameters for concentration methods with blank correction coefficients and sensitivity factors for each element

### `OrderedDict`
Utility type that maintains a collection of names and associated data in an ordered dictionary structure, used for preserving element ordering in the TUI.

## Global settings

KJ stores global configuration data in the `_KJ` dictionary, which is initialized by the `init_KJ!()` function at package load time. This dictionary contains system-wide settings and reference data that are used throughout the package for data reduction and calibration.

### Contents of `_KJ`

Most of the settings are loaded from CSV files stored in the `settings/` folder:

- **`methods`**: Available geochronology methods (U-Pb, Rb-Sr, Lu-Hf, K-Ca, Re-Os) with their default parent (P), daughter (D), and normalizing (d) isotopes. Loaded from `settings/methods.csv`.

- **`lambda`**: Decay constants and uncertainties for each geochronology method. Loaded from `settings/lambda.csv`.

- **`iratio`**: Natural isotopic abundances for all elements, organized by element and isotope. Loaded from `settings/iratio.csv`.

- **`nuclides`**: Lists of all isotopes available for each element. Loaded from `settings/nuclides.csv`.

- **`refmat`**: Reference material definitions for all supported geochronology methods, including isochron-based, point-based, and bias correction standards. Loaded from:
  - `settings/standards/isochron.csv`
  - `settings/standards/point.csv`
  - `settings/standards/bias.csv`

- **`glass`**: Standard reference material compositions for concentration measurements (e.g., NIST SRM glass standards). Loaded from `settings/glass.csv`.

- **`stoichiometry`**: Stoichiometric compositions of common minerals (e.g., apatite, zircon, calcite) used for internal standardization in concentration calculations. Loaded from `settings/stoichiometry.csv`.

- **`tree`**: The KJ tree structure used by the Text User Interface (TUI) for navigation and data organization. Initialized programmatically via `init_KJtree()`.

- **`ctrl`**: Control structure for the TUI state management (initialized as `nothing`).

- **`extensions`**: Extension system for additional functionality (initialized as `nothing`).

Individual settings can be reloaded from custom CSV files using functions like `init_methods!()`, `init_referenceMaterials!()`, `init_glass!()`, and `init_stoichiometry!()`, allowing users to customize the package configuration without modifying the source code.

## Key functions

### Data Loading and Management

#### `load`
Load LA-ICP-MS data from files in various formats (Agilent CSV, ThermoFisher CSV, or FIN2). Can read from a directory of files, a single file with timestamps, or files split into blocks. Automatically sorts samples by acquisition time and creates `Sample` objects with initial blank and signal window estimates.

#### `setGroup!`
Assign group labels to samples in a run. Groups can be assigned based on a method definition, by selecting specific sample indices, or by matching prefixes in sample names. Used to identify which samples are standards versus unknowns.

#### `summarise`/`summarize`
Create a summary DataFrame showing sample names, acquisition dates, and group assignments for all samples in a run. Useful for quickly inspecting the contents of a dataset.

### Selection Windows

#### `setBwin!`
Set blank (background) windows for one or more samples. Can be set automatically or manually specified as time ranges (in seconds) or index ranges. The blank windows define which data points are used to fit the background signal.

#### `setSwin!`
Set signal windows for one or more samples. Can be set automatically or manually specified. The signal windows define which data points contain the ablation signal of interest.

#### `sett0`
Set "time zero" (the onset of laser ablation) for samples. Can be determined automatically by detecting the signal increase or set manually. Used to properly align blank and signal windows.

### Data Processing

#### `fitBlanks`
Fit polynomial models to blank window data for each channel. Returns a DataFrame of polynomial coefficients that describe the time-dependent background signal. Default is quadratic (order 2) but can be adjusted.

#### `Gmethod`
Constructor for geochronology method definitions. Specifies the parent-daughter isotope system (U-Pb, Rb-Sr, Lu-Hf, K-Ca, or Re-Os), reference materials, isotope pairings, and polynomial orders for various corrections. Core configuration object for geochronology workflows.

#### `Cmethod`
Constructor for concentration method definitions. Specifies which elements to quantify, reference materials, and the internal standard (channel and concentration). Used for quantitative elemental analysis.

#### `Calibration!`
Initialize mass bias calibration for a geochronology method. Sets up the calibration structure using isotope ratios (typically sister/daughter isotopes) that will be fitted against reference materials.

#### `Interference`
Define an isobaric interference correction. Specifies the proxy isotope used to estimate the interference, the channel being corrected, and optional bias correction settings.

#### `REEInterference`
Define a rare earth element (REE) oxide interference correction. Specifies the REE channel, REE oxide channel, and reference materials used for calibration.

#### `getAnchor`
Retrieve the anchor coordinates for a reference material in isotope ratio space. Anchors represent the known composition of standards and are used for fractionation correction. Returns `IsochronAnchor`, `PointAnchor`, or `BiasAnchor` objects depending on the standard type.

#### `process!`
Execute the complete data reduction workflow: assign groups, detect outliers, fit blanks, fit bias corrections (for geochronology), and fit drift/downhole fractionation. Returns a `KJfit` object containing all correction parameters.

### Data Analysis

#### `atomic`
Convert raw signal intensities to corrected atomic abundances. Applies all corrections (blank subtraction, interference correction, bias correction, drift, and downhole fractionation) to produce time-resolved estimates of parent (P), daughter (D), and normalizing (d) isotope abundances.

#### `averat`
Calculate weighted mean isotope ratios for samples. For geochronology methods, returns P/D and d/D ratios with uncertainties and correlations. Uses a maximum likelihood approach to optimally combine time-resolved measurements while accounting for covariances.

#### `internochron`
Calculate internal isochron coordinates (x₀, y₀) by fitting an isochron through the time-resolved atomic abundances of individual samples. Returns intercept values with uncertainties, which can be used to correct for initial daughter isotope ratios without assuming a common value.

### Export and Visualization

#### `export2IsoplotR`
Export processed data in JSON format compatible with IsoplotR for further analysis and visualization. Can export ratios or raw data, with optional sample filtering by prefix.

#### `getInternal`
Retrieve the internal standard specification for a mineral. Returns a tuple of (channel, concentration) based on stoichiometric compositions stored in the global settings. Used for concentration calculations.

#### `concentrations`
Calculate element concentrations from calibrated LA-ICP-MS data. For single samples, returns time-resolved concentrations. For runs, returns summary statistics (mean ± standard error) for each sample.

#### `plot`
Create time-series plots of LA-ICP-MS signals or ratios. Can plot raw signals, ratios, or transformed data (log, sqrt). Displays blank and signal windows with different markers for outliers. Can overlay fitted corrections for standard samples.

#### `plotFitted!`
Overlay fitted model predictions onto an existing plot. Shows the expected signal based on the fractionation model, anchor points, and corrections. Used to visualize model quality and check for misfits.

#### `internoplot`
Create an isochron plot for internal isochron data. Plots time-resolved P/D vs d/D data with a fitted isochron line and uncertainty envelope. For U-Pb, also overlays the concordia curve. Displays calculated ages with uncertainties.

#### `plotMap`
Create a spatial map visualization of concentration or isotope ratio data. Plots measured values as colored points at their x,y coordinates (if available). Useful for visualizing spatial heterogeneity in laser ablation maps.

## TUI design

### Architecture Overview

The text-based user interface (TUI) is implemented as a three-layer architecture consisting of menu definition, message generation, and action execution:

- **TUI.jl**: Defines the menu structure and navigation logic
- **TUImessages.jl**: Generates user-facing prompts and messages
- **TUIactions.jl**: Implements the actions executed when users select menu options

### Role of _KJ["ctrl"]

The `_KJ["ctrl"]` dictionary is the central state management system for the TUI. It is initialized by `TUIinit!()` and contains all persistent state information for a TUI session:

**Session Navigation:**
- `chain`: Stack of menu nodes representing the current navigation path through the menu tree
- `template`: Boolean indicating if using a saved template workflow

**Data Objects:**
- `run`: Vector of `Sample` objects loaded from files
- `method`: Current `KJmethod` (either `Gmethod` or `Cmethod`)
- `fit`: Current `KJfit` object containing fitted parameters

**User Preferences:**
- `format`: File format (Agilent, ThermoFisher, or FIN2)
- `head2name`: Whether to extract sample names from file headers
- `multifile`: Whether data is in multiple files or single file
- `nblocks`: Number of blocks per analysis (for MC-ICP-MS data)
- `den`: Denominator channel for plotting ratios
- `transformation`: Data transformation for plotting (log, sqrt, or linear)
- `mapcolumn`/`clims`: Parameters for map visualization

**Task Prioritization:**
- `priority`: Dictionary tracking which workflow steps are required (marked with * in menu)
  - `load`: Must load data files
  - `method`: Must specify analysis method
  - `fractionation`: Must select reference materials
  - `process`: Must process data with fitted parameters

**Session Persistence:**
- `history`: DataFrame recording all menu selections for log file export
- `log`: Boolean indicating if currently replaying a log file
- `ICPpath`/`LApath`: File paths for data and laser log files
- `i`: Current sample index for viewing
- `cache`: Temporary storage for intermediate values during menu navigation

### Menu System

The menu system is driven by `_KJ["tree"]`, a dictionary where each key is a menu node containing:
- `message`: String or function generating the prompt displayed to the user
- `help`: String explaining the menu option in detail
- `action`: Either a function to execute, another node name, or a dictionary mapping user inputs to next nodes

The `dispatch!()` function implements the core navigation logic:

1. Retrieves the current menu node from `_KJ["tree"]` using the last entry in `ctrl["chain"]`
2. Displays the message (evaluating functions like those from `TUImessages.jl` to show context-dependent prompts)
3. Reads user input
4. Executes the corresponding action from `TUIactions.jl`, which modifies `ctrl` state
5. Updates `ctrl["chain"]` to navigate to the next menu node

**Example workflow:**
- User types "r" at top menu → action function `TUIread()` checks template status and returns next node
- Next node asks for file format → action `TUIformat!()` sets `ctrl["format"]` and returns guidance on next step
- User loads data in `loadICPdir` → action `TUIloadICPdir!()` calls `load()`, populates `ctrl["run"]`, and marks load priority as complete

### Message Generation

`TUImessages.jl` contains functions that generate dynamic prompts by:
- Accessing current state from `ctrl` (e.g., `ctrl["method"]`, `ctrl["nblocks"]`)
- Calling helper functions like `TUIcheck()` to show priority markers (`[*]`) for incomplete tasks
- Using `TUImethodSwitch()` to display different options based on whether a geochronology or concentration method is selected
- Building error messages and status displays

This separation allows messages to adapt automatically based on the current workflow state.

