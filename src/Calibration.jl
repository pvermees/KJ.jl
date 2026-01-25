function Calibration(;num::NamedTuple{(:ion,:channel),Tuple{String,String}}=NamedTuple(ion="",channel=""),
                      den::NamedTuple{(:ion,:channel),Tuple{String,String}}=NamedTuple(ion="",channel=""),
                      standards::AbstractSet=Set{String}())
    return Calibration(num,den,standards)
end

function REECalibration(;num::AbstractString,
                         den::AbstractString,
                         standards::AbstractSet=Set{String}())
    return REECalibration(num,den,standards)
end