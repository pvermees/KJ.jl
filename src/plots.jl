"""
    plot(samp::Sample, method::KJmethod; channels=getChannels(method), fit=nothing, num="", den="", transformation="", ...)
    plot(samp::Sample; channels=getChannels(samp), num="", den="", transformation="", offset=..., ...)

Plot time-resolved LA-ICP-MS data for a sample.

Displays signal intensities vs time with blank and signal windows marked.
Optionally overlays fitted corrections if a fit object is provided.

# Arguments
- `samp`: Sample to plot
- `method`: Method definition (optional, for determining channels and fit)
- `channels`: Channels to display
- `fit`: Fitted corrections to overlay (optional)
- `num`, `den`: Numerator/denominator for ratio plots (default: plot raw signals)
- `transformation`: Apply transformation ("log", "sqrt", or "")
- `xlim`, `ylim`: Plot limits
- `title`: Plot title
- `legend`: Legend position
- `cpalette`: Color palette
- `ms`, `ma`: Marker size and alpha
- Additional plotting options

# Returns
- Plots.jl plot object
"""
function plot(samp::Sample,
              method::KJmethod;
              channels::AbstractVector=getChannels(method),
              fit::Union{Nothing,KJfit}=nothing,
              num::AbstractString="",
              den::AbstractString="",
              transformation::AbstractString="",
              ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              linecolor::Symbol=:black,
              linestyle::Symbol=:solid,
              title::AbstractString=samp.sname*" ["*samp.group*"]",
              legend=:topleft,
              cpalette=:viridis,
              titlefontsize::Number=10,
              return_offset::Bool=false)

    offset = get_offset(samp;method=method,fit=fit,
                        channels=channels,
                        transformation=transformation,
                        num=num,den=den)

    p = plot(samp;
             channels=channels,
             num=num,den=den,
             transformation=transformation,
             offset=offset,
             ms=ms,ma=ma,xlim=xlim,ylim=ylim,
             title=title,legend=legend,cpalette=cpalette,
             titlefontsize=titlefontsize)

    if !isnothing(fit)
        if samp.group !== "sample"
            plotFitted!(p,samp,method,fit;
                        channels=channels,num=num,den=den,
                        transformation=transformation,
                        offset=offset,linecolor=linecolor,
                        linestyle=linestyle)
        end
        plotFittedBlank!(p,samp,method,fit;
                         channels=channels,num=num,den=den,
                         transformation=transformation,offset=offset,
                         linecolor=linecolor,linestyle=linestyle)
    end

    if return_offset
        return p, offset
    else
        return p
    end
end

function plot(samp::Sample;
              channels::AbstractVector=getChannels(samp),
              num::AbstractString="",
              den::AbstractString="",
              transformation::AbstractString="",
              offset::Number=get_offset(samp;transformation=transformation,num=num,den=den),
              ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              title::AbstractString=samp.sname*" ["*samp.group*"]",
              legend=:topleft,
              cpalette=:viridis,
              titlefontsize=10,
              padding::Number=0.1,
              return_offset::Bool=false)

    x, y, xlab, ylab = prep_plot(samp,channels;
                                 offset=offset,num=num,den=den,
                                 transformation=transformation)
    if ylim == :auto
        ylim = get_ylim(y,samp.swin)
    end
    p = Plots.plot(xlimits=xlim,ylimits=ylim,legend=legend)
    channels = names(y)
    cols = Plots.palette(cpalette,length(channels))
    marker = fill(:circle,length(x))
    marker[samp.dat.outlier] .= :xcross
    for i in eachindex(channels)
        Plots.scatter!(p,x,y[:,channels[i]];
                       ms=ms,ma=ma,marker=marker,
                       label=String(channels[i]),
                       markercolor=cols[i])
    end
    Plots.xlabel!(xlab)
    Plots.ylabel!(ylab)
    Plots.title!(title;titlefontsize=titlefontsize)
    buffer = (ylim[2]-ylim[1])*padding/2
    dy_win = (ylim[1] + buffer, ylim[2] - buffer)
    # plot t0:
    Plots.plot!(p,[samp.t0,samp.t0],collect(dy_win[[1,2]]);
                linecolor=:grey,linestyle=:dot,label="")
    # plot selection windows:
    for win in [samp.bwin,samp.swin]
        for w in win
            from = x[w[1]]
            to = x[w[2]]
            Plots.plot!(p,[from,from,to,to,from],
                        collect(dy_win[[1,2,2,1,1]]);
                        linecolor=:black,linestyle=:dot,label="")
        end
    end
    
    if return_offset
        return p, offset
    else
        return p
    end
end
export plot

function get_ylim(dat::AbstractDataFrame,
                  window::AbstractVector;
                  padding::Number=0.1)
    selection, x, y = windows2selection(window)
    miny, maxy = extrema(Matrix(dat[selection,:]))
    buffer = (maxy-miny)*padding
    return (miny-buffer,maxy+buffer)
end

"""
    prep_plot(samp::Sample, channels; offset=0.0, num="", den="", transformation="")

Prepare data for plotting by extracting signals and applying transformations.

# Arguments
- `samp`: Sample to extract data from
- `channels`: Channels to include
- `offset`: Offset to add before transformation
- `num`, `den`: For ratio plots
- `transformation`: Transformation to apply

# Returns
- Tuple of (x, y, xlab, ylab) where x is time, y is DataFrame of processed signals
"""
function prep_plot(samp::Sample,
                   channels::AbstractVector;
                   offset::Number=0.0,
                   num::AbstractString="",
                   den::AbstractString="",
                   transformation::AbstractString="")
    xlab = names(samp.dat)[1]
    x = samp.dat[:,xlab]
    meas = samp.dat[:,channels]
    y = transformeer(meas,transformation;
                     num=num,den=den,offset=offset)
    ratsig = (num=="" && den=="") ? "signal" : "ratio"
    if transformation == ""
        ylab = ratsig
    elseif offset > 0.0
        ylab = transformation * "(" * ratsig * "+offset)"
    else
        ylab = transformation * "(" * ratsig * ")"
    end
    return x, y, xlab, ylab
end
export prep_plot

"""
    plotFitted!(p, samp::Sample, method::KJmethod, fit::KJfit; channels=..., num="", den="", transformation="", offset=0.0, linecolor=:black, linestyle=:solid)

Overlay fitted signal predictions on an existing plot.

# Arguments
- `p`: Existing plot to overlay on
- `samp`: Sample being plotted
- `method`: Method definition
- `fit`: Fit object with correction parameters
- Additional arguments control what is plotted and how it appears
"""
function plotFitted!(p,
                     samp::Sample,
                     method::KJmethod,
                     fit::KJfit;
                     channels::AbstractVector=getChannels(method),
                     num::AbstractString="",
                     den::AbstractString="",
                     transformation::AbstractString="",
                     offset::Number=0.0,
                     linecolor::Symbol=:black,
                     linestyle=:solid)
    pred = predict(samp,method,fit;generic_names=false)
    fitted_channels = intersect(channels,names(pred))
    if length(fitted_channels) > 0
        plotFitted!(p,samp,pred[:,fitted_channels];
                    num=num,den=den,
                    transformation=transformation,
                    offset=offset,linecolor=linecolor,
                    linestyle=linestyle)
    end
end
function plotFitted!(p,
                     samp::Sample,
                     pred::AbstractDataFrame;
                     blank::Bool=false,
                     num::AbstractString="",
                     den::AbstractString="",
                     transformation::AbstractString="",
                     offset::Number=0.0,
                     linecolor::Symbol=:black,
                     linestyle=:solid)
    dat = ifelse(blank,bwinData(samp),swinData(samp))
    good = .!dat.outlier
    x = dat[good,1]
    y = transformeer(pred[good,:],transformation;
                     num=num,den=den,offset=offset)
    for yi in eachcol(y)
        Plots.plot!(p,x,yi;linecolor=linecolor,linestyle=linestyle,label="")
    end
end
export plotFitted!

"""
    plotFittedBlank!(p, samp::Sample, method::KJmethod, fit::KJfit; channels=..., num="", den="", transformation="", offset=0.0, linecolor=:black, linestyle=:solid)

Overlay fitted blank predictions on an existing plot.

# Arguments
- `p`: Existing plot to overlay on
- `samp`: Sample being plotted
- `method`: Method definition
- `fit`: Fit object with blank parameters
- Additional arguments control what is plotted and how it appears
"""
function plotFittedBlank!(p,
                          samp::Sample,
                          method::KJmethod,
                          fit::KJfit;
                          channels::AbstractVector=getChannels(method),
                          num::AbstractString="",
                          den::AbstractString="",
                          transformation::AbstractString="",
                          offset::Number=0.0,
                          linecolor::Symbol=:black,
                          linestyle::Symbol=:solid)
    pred = predict(samp,fit.blank)
    plotFitted!(p,samp,pred[:,channels];
                blank=true,num=num,den=den,
                transformation=transformation,offset=offset,
                linecolor=linecolor,linestyle=linestyle)
end
export plotFittedBlank!

"""
    plotMap(df::AbstractDataFrame, column::AbstractString; clims=(), markersize=2, markershape=:square, colorbar_scale=:log10, aspect_ratio=:equal, color=:viridis, ignore_negative=false)

Create a spatial map visualization of a data column.

Requires that the DataFrame contains x and y coordinate columns.

# Arguments
- `df`: DataFrame containing data with x, y coordinates
- `column`: Column name to visualize
- `clims`: Color scale limits (optional)
- `markersize`: Size of map points (default: 2)
- `markershape`: Shape of markers (default: :square)
- `colorbar_scale`: Scale for colorbar, typically :log10 or :identity (default: :log10)
- `aspect_ratio`: Plot aspect ratio (default: :equal)
- `color`: Color palette (default: :viridis)
- `ignore_negative`: Exclude negative values from plot (default: false)

# Returns
- Plots.jl plot object, or nothing if x,y coordinates are missing
"""
function plotMap(df::AbstractDataFrame,
                 column::AbstractString;
                 clims::Tuple=(),
                 markersize::Number=2,
                 markershape::Symbol=:square,
                 colorbar_scale::Symbol=:log10,
                 aspect_ratio::Symbol=:equal,
                 color::Symbol=:viridis,
                 ignore_negative::Bool=false)
    has_x = "x" in names(df) && !any(isnothing.(df[:,"x"]))
    has_y = "y" in names(df) && !any(isnothing.(df[:,"y"]))
    if has_x & has_y
        if ignore_negative
            selection = df[:,column] .> 0
        else
            selection = fill(true,size(df,1))
        end
        z = df[selection,column]
        if length(clims) == 0
            clims = (minimum(z),maximum(z))
        end
        p = Plots.scatter(df.x[selection],
                          df.y[selection];
                          marker_z=z,
                          color=color,
                          aspect_ratio=aspect_ratio,
                          legend=false, 
                          colorbar=true,
                          colorbar_scale=colorbar_scale,
                          clims=clims,
                          markersize=markersize,
                          markershape=markershape,
                          markerstrokewidth=0,
                          colorbar_title=column)
    else
        @warn "This dataset does not contain x and y coordinates."
        return nothing
    end
end
export plotMap
