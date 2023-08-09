include("../src/cobrexa.jl")

time_limit = 1800

organism = "iAF692"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iJR904"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iML1515"
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "e_coli_core"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iNF517"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iSB619"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iNJ661"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)

organism = "iCN900"
@show organism
# fba_data(organism)
loopless_fba_data(organism, time_limit=time_limit)
