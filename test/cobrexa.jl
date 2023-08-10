include("../src/cobrexa.jl")

time_limit = 1800

organisms = [
    "iAF692", 
    "iJR904", 
    "iML1515", 
    "e_coli_core",
    "iNF517",
    "iSB619",
    "iNJ661",
    "iCN900",
    "iAF1260",
    "iEK1008",
    "iJO1366",
    "iMM904",
    "iSDY_1059",
    "iSFV_1184",
    "iSF_1195",
    "iS_1188",
    "iSbBS512_1146"
]

for organism in organisms
    @show organism
    loopless_fba_data(organism, time_limit=time_limit)
end 
