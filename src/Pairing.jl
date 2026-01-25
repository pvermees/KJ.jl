function Pairing(;ion::AbstractString="",
                  proxy::AbstractString=ion,
                  channel::AbstractString=ion,
                  interferences::AbstractSet=Set{AbstractInterference}())
    return Pairing(ion,proxy,channel,interferences)
end