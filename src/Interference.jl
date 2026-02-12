import Base: ==, hash

==(a::Interference, b::Interference) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(Interference))
==(a::REEInterference, b::REEInterference) = all(getfield(a, f) == getfield(b, f) for f in fieldnames(REEinterference))

hash(a::Interference, h::UInt) = hash(a.ion, hash(a.proxy, hash(a.channel, h)))
hash(a::REEInterference, h::UInt) = hash(a.proxy, hash(a.REE, hash(a.REEO, h)))