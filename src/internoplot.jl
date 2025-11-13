"""
internoplot(samp::Sample,
            channels::AbstractDict,
            blank::AbstractDataFrame,
            pars::NamedTuple;
            method::Union{AbstractString,Nothing}=nothing,
            legend::Bool = false,
            nsigma::Integer = 2,
            i::Union{Integer,Nothing}=nothing,
            show_title::Bool=true,
            titlefontsize::Integer=10,
            plot_options...)

Plot an internal isochron
"""
function internoplot(samp::Sample,
                     channels::AbstractDict,
                     blank::AbstractDataFrame,
                     pars::NamedTuple;
                     method::Union{AbstractString,Nothing}=nothing,
                     legend::Bool = false,
                     nsigma::Integer = 2,
                     i::Union{Integer,Nothing}=nothing,
                     show_title::Bool=true,
                     titlefontsize::Integer=10,
                     plot_options...)
    Phat, Dhat, dhat = atomic(samp,channels,blank,pars)
    x0, sx0, y0, sy0, rx0y0 = internochron(samp,channels,blank,pars)
    E = [ [sx0^2 rx0y0*sx0*sy0 ]
          [rx0y0*sx0*sy0 sy0^2 ] ]
    if isnothing(method)
        p = internoplot(x0,y0,E,Phat,Dhat,dhat;
                        legend=legend,nsigma=nsigma,plot_options...)
    else
        Pname,Dname,dname = getPDd(method)
        xlab = Pname * "/" * Dname
        ylab = dname * "/" * Dname
        p = internoplot(x0,y0,E,Phat,Dhat,dhat;
                        legend=legend,nsigma=nsigma,
                        xlab=xlab,ylab=ylab,plot_options...)
        t, st, y0, sy0, rho = x0y02t(x0,y0,E,method)
        sdig = 2
        tdig = ceil(Int,log10(t/st)) + sdig
        tstring = "t = " *
            string(round(t,sigdigits=tdig)) * "±" *
            string(round(st,sigdigits=sdig)) * "Ma"
        Plots.annotate!(Plots.xlims(p)[2],
                        Plots.ylims(p)[2],
                        Plots.text(tstring,:right,titlefontsize))
        if method=="U-Pb"
            xmax = Plots.xlims(p)[2]
            ymax = Plots.ylims(p)[2]
            add_concordia_line(xmax,ymax)
        end
    end
    if show_title
        title = samp.sname*" ["*samp.group*"]"
        if !isnothing(i)
            title = string(i) * ". " * title
        end
        Plots.title!(title;titlefontsize=titlefontsize)
    end
    return p
end
"""
internoplot(x0::AbstractFloat,
            y0::AbstractFloat,
            E::Matrix,
            Phat::AbstractVector,
            Dhat::AbstractVector,
            dhat::AbstractVector;
            xlim::AbstractVector = [0,1.05*maximum([Phat./Dhat;x0+sqrt(E[1,1])])],
            ylim::AbstractVector = [0,1.05*maximum([dhat./Dhat;y0+sqrt(E[2,2])])],
            legend::Bool = false,
            nsigma::Integer = 2,
            xlab::AbstractString = "P/D",
            ylab::AbstractString = "d/D",
            plot_options...)
"""
function internoplot(x0::AbstractFloat,
                     y0::AbstractFloat,
                     E::Matrix,
                     Phat::AbstractVector,
                     Dhat::AbstractVector,
                     dhat::AbstractVector;
                     xlim::AbstractVector = [0,1.05*maximum([Phat./Dhat;x0+sqrt(E[1,1])])],
                     ylim::AbstractVector = [0,1.05*maximum([dhat./Dhat;y0+sqrt(E[2,2])])],
                     legend::Bool = false,
                     nsigma::Integer = 2,
                     xlab::AbstractString = "P/D",
                     ylab::AbstractString = "d/D",
                     plot_options...)
    nstep = 50
    x = range(xlim[1],xlim[2],nstep)
    y = @. y0 - x*y0/x0
    J = zeros(nstep,2)
    J[:,1] = x.*y0/x0^2
    J[:,2] = 1 .- x/x0
    covmat = J * E * transpose(J)
    sy = sqrt.(diag(covmat))
    p = Plots.plot(x,y,ribbon=nsigma*sy;legend=legend,xlim=xlim,ylim=ylim,plot_options...)
    Plots.plot!([0,x0],[y0,0];seriescolor=:black,legend=legend)
    Plots.plot!(Phat./Dhat,dhat./Dhat;seriestype=:scatter,legend=legend,plot_options...)
    Plots.xlabel!(xlab)
    Plots.ylabel!(ylab)
    return p
end
export internoplot

function add_concordia_line(xmax,ymax)
    L5, L8, U58 = UPb_helper()
    tmin = log(1+1/xmax)/L8
    function Pb76misfit(par)
        t = par[1]
        ypred = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
        return (ymax - ypred)^2
    end
    Pb76fit = Optim.optimize(Pb76misfit,[4000.0])
    tmax = Optim.minimizer(Pb76fit)[1]
    t = range(tmin[1],tmax;step=50)
    x = @. 1/(exp(L8*t)-1)
    y = @. U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    Plots.plot!(x,y,linewidth=1.5,linecolor=:black)
end
