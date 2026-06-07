"""
    KJfit(method::Gmethod)
    KJfit(method::Cmethod)

Create a default fit object for the supplied method type.
"""
function KJfit(method::Gmethod)
    return Gfit(method)
end

function KJfit(method::Cmethod)
    return Cfit()
end

"""
    Gfit(method::Gmethod; blank=DataFrame(), drift=zeros(method.ndrift), down=zeros(method.ndown), adrift=drift, covmat=..., bias=Dict())

Construct a geochronology fit object with optional initial parameters.
"""
function Gfit(method::Gmethod;
              blank::AbstractDataFrame = DataFrame(),
              drift::AbstractVector = zeros(method.ndrift),
              down::AbstractVector = zeros(method.ndown),
              adrift::AbstractVector = drift,
              covmat::AbstractMatrix = zeros(length([drift,down]),
                                             length([drift,down])),
              bias::AbstractDict = Dict())
    return Gfit(blank,drift,down,adrift,covmat,bias)
end

"""
    Cfit()

Construct an empty concentration fit container.
"""
function Cfit()
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end