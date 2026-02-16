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

    p, offset = plot(samp;
                     channels=channels,
                     num=num,den=den,
                     transformation=transformation,
                     ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                     title=title,legend=legend,cpalette=cpalette,
                     titlefontsize=titlefontsize,
                     return_offset=true)

    if !isnothing(fit)
        if samp.group !== "sample"
            plotFitted!(p,samp,method,fit;
                        channels=channels,
                        num=num,den=den,
                        transformation=transformation,
                        offset=offset,linecolor=linecolor,
                        linestyle=linestyle)
        end
        plotFittedBlank!(p,samp,method,fit;
                         num=num,den=den,
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
              ms::Number=2,ma::Number=0.5,
              xlim=:auto,ylim=:auto,
              title::AbstractString=samp.sname*" ["*samp.group*"]",
              legend=:topleft,
              cpalette=:viridis,
              titlefontsize=10,
              padding::Number=0.1,
              return_offset::Bool=false)

    x, y, xlab, ylab, offset = prep_plot(samp,channels;
                                         num=num,den=den,
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

function prep_plot(samp::Sample,
                   channels::AbstractVector;
                   num::AbstractString="",
                   den::AbstractString="",
                   transformation::AbstractString="")
    xlab = names(samp.dat)[1]
    x = samp.dat[:,xlab]
    meas = samp.dat[:,channels]
    offset = get_offset(meas;
                        transformation=transformation,
                        num=num,den=den)
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
    return x, y, xlab, ylab, offset
end
export prep_plot

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

function plotFittedBlank!(p,
                          samp::Sample,
                          method::KJmethod,
                          fit::KJfit;
                          num::AbstractString="",
                          den::AbstractString="",
                          transformation::AbstractString="",
                          offset::Number=0.0,
                          linecolor::Symbol=:black,
                          linestyle::Symbol=:solid)
    channels = getChannels(method)
    pred = predict(samp,fit.blank[:,channels])
    plotFitted!(p,samp,pred;
                blank=true,num=num,den=den,
                transformation=transformation,offset=offset,
                linecolor=linecolor,linestyle=linestyle)
end
export plotFittedBlank!

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
