using Plasmatrace, PTpost, Plots, Infiltrator

myrun = load("/home/pvermees/Documents/Adelaide/own_data/Gt_Lu_Hf_lines/",
             instrument="Agilent")
method = "Lu-Hf"
channels = Dict("d"=>"Hf178 -> 260",
                "D"=>"Hf176 -> 258",
                "P"=>"Lu175 -> 175")
standards = Dict("GWA-2_gt" => "90667")
glass = Dict("NIST610" => "NIST610")
blk, fit = Plasmatrace.process!(myrun,method,channels,standards,glass,
                                nblank=2,ndrift=1,ndown=0)
istand = 15
p1 = Plasmatrace.plot(myrun[istand],method,channels,blk,fit,standards,glass,
                      den="Hf176 -> 258",transformation="log",ylim=[-4,13])

isamp = 7
P, D, d = Plasmatrace.atomic(myrun[isamp],channels,blk,fit)
x0, y0, E = PTpost.internochron(P,D,d;numerical=true)
p2 = PTpost.plot(x0,y0,E,P,D,d;markercolor=:white)#,ylim=[0,10])

p = Plots.plot(p1,p2,layout=(1,2))

display(p)
