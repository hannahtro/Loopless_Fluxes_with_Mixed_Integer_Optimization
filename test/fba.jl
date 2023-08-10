include("../src/fba.jl")

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
    get_fba_data(organism, save_lp=false, json=true)
end 

# get_fba_data("iAF692", save_lp=false, json=true)