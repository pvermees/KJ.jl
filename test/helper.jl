function synthetic()
    nb = 20
    ns = 100
    out = [
    Sample("BP_1", # sname
           DateTime("2025-01-01T08:00:00"), # datetime
           DataFrame("Lu175 -> 175"=P,
                     "Hf176 -> 258"=D,
                     "Hf178 -> 260"=d), # dat
           10.0, # t0
           [(0.0,9.0)], # bwin
           [(11.0,30.0)], # swin
           "sample"),
    Sample("HOG_1",
           DateTime("2025-01-01T08:00:00"),
           DataFrame("Lu175 -> 175"=P,
                     "Hf176 -> 258"=D,
                     "Hf178 -> 260"=d),
           10.0,
           [(0.0,9.0)],
           [(11.0,30.0)],
           "sample"),
    Sample("NIST612_1",
           DateTime("2025-01-01T08:00:00"),
           DataFrame("Lu175 -> 175"=P,
                     "Hf176 -> 258"=D,
                     "Hf178 -> 260"=d),
           10.0,
           [(0.0,9.0)],
           [(11.0,30.0)],
           "sample"),
    Sample("BP_2",
           DateTime("2025-01-01T08:00:00"),
           DataFrame("Lu175 -> 175"=P,
                     "Hf176 -> 258"=D,
                     "Hf178 -> 260"=d),
           10.0,
           [(0.0,9.0)],
           [(11.0,30.0)],
           "sample")
    ]
    
end
