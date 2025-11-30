function plot(samp::Sample,
              method::KJmethod,
              fit::KJfit;
              num=nothing,den=nothing,
              transformation=nothing,
              ms=2,ma=0.5,
              xlim=:auto,ylim=:auto,
              linecol="black",
              linestyle=:solid,
              i=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              return_offset::Bool=false)

    channelvec = collect(values(getChannels(method)))
    
    p, offset = plot(samp;
                     channels=channelvec,
                     num=num,den=den,
                     transformation=transformation,
                     ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                     i=i,legend=legend,cpalette=cpalette,
                     show_title=show_title,
                     titlefontsize=titlefontsize,
                     return_offset=true)
    
    if samp.group != "sample"

        plotFitted!(p,samp,method,fit;
                    num=num,den=den,transformation=transformation,
                    offset=offset,linecolor=linecol,linestyle=linestyle)
        
    end

    plotFittedBlank!(p,samp,method,fit;
                     num=num,den=den,
                     transformation=transformation,offset=offset,
                     linecolor=linecol,linestyle=linestyle)

    if return_offset
        return p, offset
    else
        return p
    end
end

function plot(samp::Sample;
              channels::AbstractVector=getChannels(samp),
              num::Union{Nothing,AbstractString}=nothing,
              den::Union{Nothing,AbstractString}=nothing,
              transformation::Union{Nothing,AbstractString}=nothing,
              ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              i::Union{Nothing,Integer}=nothing,
              legend=:topleft,
              cpalette=:viridis,
              show_title=true,
              titlefontsize=10,
              padding::Number=0.1,
              return_offset::Bool=false)

    x, y, xlab, ylab, offset = prep_plot(samp,channels;
                                         num=num,den=den,ylim=ylim,
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
    if show_title
        title = samp.sname*" ["*samp.group*"]"
        if !isnothing(i)
            title = string(i) * ". " * title
        end
        Plots.title!(title;titlefontsize=titlefontsize)
    end
    buffer = (ylim[2]-ylim[1])*padding/2
    dy_win = (ylim[1] + buffer, ylim[2] - buffer)
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
    ratsig = isnothing(num) & isnothing(den) ? "signal" : "ratio"
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
    return x, ty, xlab, ylab, offset
end

function plotFitted!(p,
                     samp::Sample,
                     method::KJmethod,
                     fit::KJfit;
                     num::Union{Nothing,AbstractString}=nothing,
                     den::Union{Nothing,AbstractString}=nothing,
                     transformation::Union{Nothing,AbstractString}=nothing,
                     offset::Union{Nothing,Number}=nothing,
                     linecolor="black",
                     linestyle=:solid)
    pred = predict(samp,method,fit)
    rename!(pred,getChannels(method))
    plotFitted!(p,samp,pred;
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,
                linestyle=linestyle)
end

# helper
function plotFitted!(p,
                     samp::Sample,
                     pred::AbstractDataFrame;
                     blank::Bool=false,
                     num::Union{Nothing,AbstractString}=nothing,
                     den::Union{Nothing,AbstractString}=nothing,
                     transformation::Union{Nothing,AbstractString}=nothing,
                     offset::Union{Nothing,Number}=nothing,
                     linecolor="black",linestyle=:solid)
    dat = ifelse(blank,bwinData(samp),swinData(samp))
    good = .!dat.outlier
    x = dat[good,1]
    y = formRatios(pred[good,:],num,den)
    ty, offset = transformeer(y,transformation;offset=offset)
    for tyi in eachcol(ty)
        Plots.plot!(p,x,tyi;linecolor=linecolor,linestyle=linestyle,label="")
    end
end
export plotFitted!

function plotFittedBlank!(p,
                          samp::Sample,
                          method::KJmethod,
                          fit::KJfit;
                          num::Union{Nothing,AbstractString}=nothing,
                          den::Union{Nothing,AbstractString}=nothing,
                          transformation::Union{Nothing,AbstractString}=nothing,
                          offset::Union{Nothing,Number}=0.0,
                          linecolor="black",
                          linestyle::Symbol=:solid)
    channels = collect(values(getChannels(method)))
    pred = predict(samp,fit.blank[:,channels])
    plotFitted!(p,samp,pred;
                blank=true,num=num,den=den,transformation=transformation,
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
                 ignore_negative::Bool=true)
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
