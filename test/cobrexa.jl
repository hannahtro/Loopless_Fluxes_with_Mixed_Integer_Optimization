include("../src/cobrexa.jl")

organism = "iAF692"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iJR904"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iML1515"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "e_coli_core"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iNF517"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

#organism = "iSB619"
#@show organism
#fba_data(organism)
#loopless_fba_data(organism, time_limit=600)

organism = "iNJ661"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iCN900"
@show organism
fba_data(organism)
loopless_fba_data(organism, time_limit=600)
