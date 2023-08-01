include("../src/cobrexa.jl")

organism = "iAF692"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iJR904"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iML1515"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "e_coli_core"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iNF517"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iSB619"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iNJ661"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)

organism = "iCN900"
fba_data(organism)
loopless_fba_data(organism, time_limit=600)
