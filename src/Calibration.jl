function Calibration!(method::Gmethod;
                      standards::Set=Set{String}())
    method.bias = Calibration(num=(ion=method.d.proxy,channel=method.d.channel),
                              den=(ion=method.D.proxy,channel=method.D.channel),
                              standards=standards)
end
export Calibration!