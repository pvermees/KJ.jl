# ratios
function plot(samp::Sample,
              method::AbstractString,
              channels::AbstractDict,
              blank::AbstractDataFrame,
              pars::NamedTuple,
              standards::Union{AbstractDict,AbstractVector},
              glass::Union{AbstractDict,AbstractVector};
              dt::Union{AbstractDict,Nothing}=nothing,
              dead::AbstractFloat=0.0,
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,xlim=:auto,ylim=:auto,
              linecol="black",linestyle=:solid,
              i=nothing,legend=:topleft,
              show_title=true,
              titlefontsize=10,
              kw...)
    Sanchors = getAnchors(method,standards,false)
    Ganchors = getAnchors(method,glass,true)
    anchors = merge(Sanchors,Ganchors)
    return plot(samp,channels,blank,pars,anchors;
                dt=dt,dead=dead,
                num=num,den=den,transformation=transformation,
                seriestype=seriestype,
                ms=ms,ma=ma,xlim=xlim,ylim=ylim,i=i,
                legend=legend,show_title=show_title,
                titlefontsize=titlefontsize,
                kw...)
end
function plot(samp::Sample,
              channels::AbstractDict,
              blank::AbstractDataFrame,
              pars::NamedTuple,
              anchors::AbstractDict;
              dt::Union{AbstractDict,Nothing}=nothing,
              dead::AbstractFloat=0.0,
              num=nothing,den=nothing,
              transformation=nothing,
              seriestype=:scatter,
              ms=2,ma=0.5,
              xlim=:auto,ylim=:auto,
              linecol="black",
              linestyle=:solid,
              i=nothing,
              legend=:topleft,
              show_title=true,
              titlefontsize=10,
              kw...)
    if samp.group == "sample"

        p = plot(samp;
                 channels=collect(values(channels)),
                 num=num,den=den,transformation=transformation,
                 seriestype=seriestype,ms=ms,ma=ma,
                 xlim=xlim,ylim=ylim,i=i,
                 legend=legend,show_title=show_title,
                 titlefontsize=titlefontsize,kw...)
        
    else

        offset = getOffset(samp,channels,blank,pars,anchors,transformation;
                           dt=dt,dead=dead,num=num,den=den)

        p = plot(samp;
                 channels=collect(values(channels)),
                 num=num,den=den,transformation=transformation,offset=offset,
                 seriestype=seriestype,ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                 i=i,legend=legend,show_title=show_title,
                 titlefontsize=titlefontsize,kw...)

        plotFitted!(p,samp,blank,pars,channels,anchors;
                    dt=dt,dead=dead,num=num,den=den,transformation=transformation,
                    offset=offset,linecolor=linecol,linestyle=linestyle)
        
    end
    return p
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
              legend=:topleft,show_title=true,
              titlefontsize=10,kw...)
    if samp.group == "sample"

        p = plot(samp;
                 num=num,den=den,transformation=transformation,
                 seriestype=seriestype,ms=ms,ma=ma,
                 xlim=xlim,ylim=ylim,i=i,
                 legend=legend,show_title=show_title,
                 titlefontsize=titlefontsize,kw...)
        
    else

        offset = getOffset(samp,blank,pars,elements,internal,transformation;
                           num=num,den=den)

        p = plot(samp;
                 num=num,den=den,transformation=transformation,offset=offset,
                 seriestype=seriestype,ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                 i=i,legend=legend,show_title=show_title,
                 titlefontsize=titlefontsize,kw...)

        plotFitted!(p,samp,blank,pars,elements,internal;
                     num=num,den=den,transformation=transformation,
                     offset=offset,linecolor=linecol,linestyle=linestyle)
        
    end
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
              legend=:topleft,show_title=true,
              titlefontsize=10,kw...)
    elements = channels2elements(samp)
    return plot(samp,blank,pars,elements,internal;
                num=num,den=den,transformation=transformation,
                seriestype=seriestype,ms=ms,ma=ma,xlim=xlim,ylim=ylim,
                linecol=linecol,linestyle=linestyle,i=i,
                legend=legend,show_title=show_title,
                titlefontsize=titlefontsize,kw...)
end
function plot(samp::Sample;
              channels::AbstractVector=getChannels(samp),
              num=nothing,
              den=nothing,
              transformation=nothing,
              offset=nothing,
              seriestype=:scatter,ms=2,ma=0.5,
              xlim=:auto,ylim=:auto,
              i::Union{Nothing,Integer}=nothing,
              legend=:topleft,
              show_title=true,
              titlefontsize=10,
              kw...)
    xlab = names(samp.dat)[1]
    x = samp.dat[:,xlab]
    meas = samp.dat[:,channels]
    y = (isnothing(num) && isnothing(den)) ? meas : formRatios(meas,num,den)
    if isnothing(offset)
        offset = Dict(zip(names(y),fill(0.0,size(y,2))))
    end
    ty = transformeer(y;transformation=transformation,offset=offset)
    ratsig = isnothing(den) ? "signal" : "ratio"
    ylab = isnothing(transformation) ? ratsig : transformation*"("*ratsig*")"
    p = Plots.plot(x,Matrix(ty);
                   ms=ms,ma=ma,seriestype=seriestype,
                   label=permutedims(names(y)),
                   legend=legend,xlimits=xlim,ylimits=ylim,
                   kw...)
    Plots.xlabel!(xlab)
    Plots.ylabel!(ylab)
    if show_title
        title = samp.sname*" ["*samp.group*"]"
        if !isnothing(i)
            title = string(i) * ". " * title
        end
        Plots.title!(title;titlefontsize=titlefontsize)
    end
    dy = Plots.ylims(p)
    # plot t0:
    Plots.plot!(p,[samp.t0,samp.t0],collect(dy[[1,2]]);
                linecolor="grey",linestyle=:dot,label="")
    # plot selection windows:
    for win in [samp.bwin,samp.swin]
        for w in win
            from = x[w[1]]
            to = x[w[2]]
            Plots.plot!(p,[from,from,to,to,from],collect(dy[[1,2,2,1,1]]);
                        linecolor="black",linestyle=:dot,label="")
        end
    end
    return p
end
export plot

# minerals
function plotFitted!(p,
                     samp::Sample,
                     blank::AbstractDataFrame,
                     pars::NamedTuple,
                     channels::AbstractDict,
                     anchors::AbstractDict;
                     dt::Union{AbstractDict,Nothing}=nothing,
                     dead::AbstractFloat=0.0,
                     num=nothing,den=nothing,transformation=nothing,
                     offset::AbstractDict,linecolor="black",linestyle=:solid)
    pred = predict(samp,pars,blank,channels,anchors;
                   dt=dt,dead=dead)
    rename!(pred,[channels[i] for i in names(pred)])
    plotFitted!(p,samp,pred;
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,linestyle=linestyle)
end
# concentrations
function plotFitted!(p,
                     samp::Sample,
                     blank::AbstractDataFrame,
                     pars::AbstractVector,
                     elements::AbstractDataFrame,
                     internal::AbstractString;
                     num=nothing,den=nothing,transformation=nothing,
                     offset::AbstractDict,linecolor="black",linestyle=:solid)
    pred = predict(samp,pars,blank,elements,internal)
    plotFitted!(p,samp,pred;
                num=num,den=den,transformation=transformation,
                offset=offset,linecolor=linecolor,linestyle=linestyle)
end
# helper
function plotFitted!(p,
                     samp::Sample,
                     pred::AbstractDataFrame;
                     num=nothing,den=nothing,transformation=nothing,
                     offset::AbstractDict,linecolor="black",linestyle=:solid)
    x = windowData(samp,signal=true)[:,1]
    y = formRatios(pred,num,den)
    ty = transformeer(y;transformation=transformation,offset=offset)
    for tyi in eachcol(ty)
        Plots.plot!(p,x,tyi;linecolor=linecolor,linestyle=linestyle,label="")
    end
end
export plotFitted!
