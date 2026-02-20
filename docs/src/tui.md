# Text-Based User Interface (TUI)

KJ provides an interactive text-based user interface for data processing workflows.

## Starting the TUI

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

## Switching Between TUI and REPL

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

For detailed function documentation, see the [API reference](api.md).
