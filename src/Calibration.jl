"""
    Calibration!(method::Gmethod; standards=Set{String}())

Initialize the calibration (bias correction) for a geochronology method.

This sets up the bias calibration structure using the sister/daughter isotope
ratio, which will later be fitted using reference materials.

# Arguments
- `method`: Geochronology method to set up calibration for
- `standards`: Set of standard names to use for calibration (default: empty set)
"""
function Calibration!(method::Gmethod;
                      standards::Set=Set{String}())
    method.bias = Calibration(num=(ion=method.d.proxy,channel=method.d.channel),
                              den=(ion=method.D.proxy,channel=method.D.channel),
                              standards=standards)
end
export Calibration!