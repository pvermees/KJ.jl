# ratios
function plot(samp::Sample,
              method::AbstractString,
              channels::AbstractDict,
              blank::AbstractDataFrame,
              pars::NamedTuple,
              standards::Union{AbstractDict,AbstractVector},
              glass::Union{AbstractDict,AbstractVector};
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,xlim=:auto,ylim=:auto,
              linecol="black",linestyle=:solid,
              i=nothing,legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              kw...)
    anchors = getAnchors(method,standards,glass)
    return plot(samp,channels,blank,pars,anchors;
                num=num,den=den,transformation=transformation,
                seriestype=seriestype,
                ms=ms,ma=ma,xlim=xlim,ylim=ylim,i=i,
                legend=legend,cpalette=cpalette,
                show_title=show_title,
                titlefontsize=titlefontsize,
                kw...)
end
function plot(samp::Sample,
              channels::AbstractDict,
              blank::AbstractDataFrame,
              pars::NamedTuple,
              anchors::AbstractDict;
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,
              xlim=:auto,ylim=:auto,
              linecol="black",
              linestyle=:solid,
              i=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              return_offset::Bool=false,
              kw...)

    channelvec = collect(values(channels))

    x, y, xlab, ylab, ylim, offset = prep_plot(samp,channelvec;
                                               num=num,den=den,
                                               ylim=ylim,
                                               transformation=transformation)
    
    p = plot(samp,x,y;
             ms=ms,ma=ma,seriestype=seriestype,
             label=permutedims(names(y)),
             xlimits=xlim,ylimits=ylim,
             i=i,legend=legend,cpalette=cpalette,
             show_title=show_title,
             titlefontsize=titlefontsize,
             xlab=xlab,ylab=ylab,kw...)
    
    if samp.group != "sample"

        plotFitted!(p,samp,blank,pars,channels,anchors;
                    num=num,den=den,transformation=transformation,
                    offset=offset,linecolor=linecol,linestyle=linestyle)
        
    end

    plotFittedBlank!(p,samp,blank,channelvec;
                     num=num,den=den,
                     transformation=transformation,offset=offset,
                     linecolor=linecol,linestyle=linestyle)

    if return_offset
        return p, offset
    else
        return p
    end
end
# concentrations
function plot(samp::Sample,
              blank::AbstractDataFrame,
              pars::AbstractVector,
              elements::AbstractDataFrame,
              internal::AbstractString;
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,xlim=:auto,ylim=:auto,
              linecol="black",linestyle=:solid,i=nothing,
              legend=:topleft,cpalette=:viridis,
              show_title=true,
              titlefontsize=10,kw...)

    channelvec = collect(getChannels(samp))

    x, y, xlab, ylab, ylim, offset = prep_plot(samp,channelvec;
                                               num=num,den=den,
                                               ylim=ylim,
                                               transformation=transformation)
    p = plot(samp;
             num=num,den=den,transformation=transformation,
             seriestype=seriestype,ms=ms,ma=ma,
             xlim=xlim,ylim=ylim,i=i,
             legend=legend,cpalette=cpalette,
             show_title=show_title,
             titlefontsize=titlefontsize,
             xlab=xlab,ylab=ylab,kw...)

    if samp.group != "sample"

        plotFitted!(p,samp,blank,pars,elements,internal;
                    num=num,den=den,transformation=transformation,
                    offset=offset,linecolor=linecol,linestyle=linestyle)
        
    end
    
    plotFittedBlank!(p,samp,blank;
                     num=num,den=den,transformation=transformation,
                     offset=offset,linecolor=linecol,linestyle=linestyle)

    return p
end
function plot(samp::Sample,
              blank::AbstractDataFrame,
              pars::AbstractVector,
              internal::AbstractString;
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,xlim=:auto,ylim=:auto,
              linecol="black",linestyle=:solid,i=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,kw...)
    elements = channels2elements(samp)
    return plot(samp,blank,pars,elements,internal;
                num=num,den=den,
                transformation=transformation,
                seriestype=seriestype,ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                linecol=linecol,linestyle=linestyle,i=i,
                legend=legend,cpalette=cpalette,
                show_title=show_title,
                titlefontsize=titlefontsize,kw...)
end
function plot(samp::Sample;
              channels::AbstractVector=getChannels(samp),
              num::Union{Nothing,AbstractString}=nothing,
              den::Union{Nothing,AbstractString}=nothing,
              transformation::Union{Nothing,AbstractString}=nothing,
              seriestype=:scatter,ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              i::Union{Nothing,Integer}=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              padding::Number=0.1,
              kw...)

    x, y, xlab, ylab, ylim, offset = prep_plot(samp,channels;
                                               num=num,den=den,ylim=ylim,
                                               transformation=transformation)

    p = plot(samp,x,y;
             seriestype=seriestype,ms=ms,ma=ma,
             xlim=xlim,ylim=ylim,i=i,legend=legend,
             cpalette=cpalette,show_title=show_title,
             titlefontsize=titlefontsize,
             padding=padding,xlab=xlab,ylab=ylab,kw...)

    return p
end
function plot(samp::Sample,
              x::AbstractVector,
              y::AbstractDataFrame;
              xlab::AbstractString="x",
              ylab::AbstractString="y",
              ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              i::Union{Nothing,Integer}=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              padding::Number=0.1,
              kw...)
    p = Plots.plot(xlimits=xlim,ylimits=ylim,legend=legend)
    channels = names(y)
    cols = Plots.palette(cpalette,length(channels))
    for i in eachindex(channels)
        Plots.scatter!(p,x,y[:,channels[i]];
                       ms=ms,ma=ma,label=channels[i],
                       markercolor=cols[i],kw...)
    end
    Plots.xlabel!(xlab)
    Plots.ylabel!(ylab)
    if show_title
        title = samp.sname*" ["*samp.group*"]"
        if !isnothing(i)
            title = string(i) * ". " * title
        end
        Plots.title!(title;titlefontsize=titlefontsize)
    end
    if ylim == :auto
        dy_win = collect(Plots.ylims(p))
    else
        buffer = (ylim[2]-ylim[1])*padding/2
        dy_win = (ylim[1] + buffer, ylim[2] - buffer)
    end
    # plot t0:
    Plots.plot!(p,[samp.t0,samp.t0],collect(dy_win[[1,2]]);
                linecolor="grey",linestyle=:dot,label="")
    # plot selection windows:
    for win in [samp.bwin,samp.swin]
        for w in win
            from = x[w[1]]
            to = x[w[2]]
            Plots.plot!(p,[from,from,to,to,from],
                        collect(dy_win[[1,2,2,1,1]]);
                        linecolor="black",linestyle=:dot,label="")
        end
    end
    return p
end
export plot

function prep_plot(samp::Sample,
                   channels::AbstractVector;
                   num::Union{Nothing,AbstractString}=nothing,
                   den::Union{Nothing,AbstractString}=nothing,
                   ylim=:auto,
                   transformation::Union{Nothing,AbstractString}=nothing,
                   padding::Number=0.1)
    xlab = names(samp.dat)[1]
    x = samp.dat[:,xlab]
    meas = samp.dat[:,channels]
    ratsig = isnothing(den) ? "signal" : "ratio"
    y = (ratsig == "signal") ? meas : formRatios(meas,num,den)    
    arg = nothing
    min_val = minimum(Matrix(y))
    if isnothing(transformation)
        ylab = ratsig
    elseif (transformation == "log" && min_val <= 0) ||
        (transformation == "sqrt" && min_val < 0)
        ylab = transformation * "(" * ratsig * "+offset)"
    else
        ylab = transformation*"("*ratsig*")"
    end
    ty, offset = transformeer(y,transformation)
    if ylim == :auto && ratsig == "ratio"
        ylim = get_ylim(ty,samp.swin)
    end
    return x, ty, xlab, ylab, ylim, offset
end
export prep_plot

function get_ylim(dat::AbstractDataFrame,
                  window::AbstractVector;
                  padding::Number=0.1)
    selection, x, y = windows2selection(window)
    miny, maxy = extrema(Matrix(dat[selection,:]))
    buffer = (maxy-miny)*padding
    return (miny-buffer,maxy+buffer)
end

# minerals
function plotFitted!(p,
                     samp::Sample,
                     blank::AbstractDataFrame,
                     pars::NamedTuple,
                     channels::AbstractDict,
                     anchors::AbstractDict;
                     num::Union{Nothing,AbstractString}=nothing,
                     den::Union{Nothing,AbstractString}=nothing,
                     transformation::Union{Nothing,AbstractString}=nothing,
                     offset::Number=0.0,
                     linecolor="black",
                     linestyle=:solid,
                     debug::Bool=false)
    pred = predict(samp,pars,blank,channels,anchors)
    rename!(pred,[channels[i] for i in names(pred)])
    plotFitted!(p,samp,pred;
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,
                linestyle=linestyle,debug=debug)
end
# concentrations
function plotFitted!(p,
                     samp::Sample,
                     blank::AbstractDataFrame,
                     pars::AbstractVector,
                     elements::AbstractDataFrame,
                     internal::AbstractString;
                     num::Union{Nothing,AbstractString}=nothing,
                     den::Union{Nothing,AbstractString}=nothing,
                     transformation::Union{Nothing,AbstractString}=nothing,
                     offset::Union{Nothing,Number}=nothing,
                     linecolor="black",
                     linestyle=:solid,
                     debug::Bool=false)
    pred = predict(samp,pars,blank,elements,internal)
    plotFitted!(p,samp,pred;
                num=num,den=den,
                transformation=transformation,offset=offset,
                linecolor=linecolor,linestyle=linestyle)
end
# helper
function plotFitted!(p,
                     samp::Sample,
                     pred::AbstractDataFrame;
                     blank::Bool=false,signal::Bool=true,
                     num::Union{Nothing,AbstractString}=nothing,
                     den::Union{Nothing,AbstractString}=nothing,
                     transformation::Union{Nothing,AbstractString}=nothing,
                     offset::Union{Nothing,Number}=nothing,
                     linecolor="black",linestyle=:solid,
                     debug::Bool=false)
    x = windowData(samp,blank=blank,signal=signal)[:,1]
    y = formRatios(pred,num,den)
    ty, offset = transformeer(y,transformation;
                              offset=offset,debug=debug)
    for tyi in eachcol(ty)
        Plots.plot!(p,x,tyi;linecolor=linecolor,linestyle=linestyle,label="")
    end
end
export plotFitted!

# minerals
function plotFittedBlank!(p,
                          samp::Sample,
                          blank::AbstractDataFrame,
                          channels::AbstractVector;
                          num::Union{Nothing,AbstractString}=nothing,
                          den::Union{Nothing,AbstractString}=nothing,
                          transformation::Union{Nothing,AbstractString}=nothing,
                          offset::Union{Nothing,Number}=0.0,
                          linecolor="black",
                          linestyle::Symbol=:solid)
    pred = predict(samp,blank[:,channels])
    plotFitted!(p,samp,pred;
                blank=true,signal=false,
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,linestyle=linestyle)
end
# concentrations
function plotFittedBlank!(p,
                          samp::Sample,
                          blank::AbstractDataFrame;
                          num::Union{Nothing,AbstractString}=nothing,
                          den::Union{Nothing,AbstractString}=nothing,
                          transformation::Union{Nothing,AbstractString}=nothing,
                          offset::Union{Nothing,Number}=0.0,
                          linecolor::AbstractString="black",
                          linestyle::Symbol=:solid)
    pred = predict(samp,blank)
    plotFitted!(p,samp,pred;
                blank=true,signal=false,
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,linestyle=linestyle)
end
export plotFittedBlank!

function plotMap(df::AbstractDataFrame,
                 column::AbstractString;
                 clims::Union{Nothing,Tuple}=nothing,
                 markersize::Number=2,
                 markershape::Symbol=:square,
                 colorbar_scale::Symbol=:log10,
                 aspect_ratio::Symbol=:equal,
                 color::Symbol=:viridis,
                 ignore_negative::Bool=true,
                 kw...)
    has_x = "x" in names(df) && !any(isnothing.(df[:,"x"]))
    has_y = "y" in names(df) && !any(isnothing.(df[:,"y"]))
    if has_x & has_y
        if ignore_negative
            selection = df[:,column] .> 0
        else
            selection = fill(true,size(df,1))
        end
        z = df[selection,column]
        if isnothing(clims)
            clims = (minimum(z),maximum(z))
        end
        p = Plots.plot(df.x[selection],
                       df.y[selection];
                       seriestype=:scatter,
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
                       colorbar_title=column,
                       kw...)
    else
        @warn "This dataset does not contain x and y coordinates."
        return nothing
    end
end
export plotMap
