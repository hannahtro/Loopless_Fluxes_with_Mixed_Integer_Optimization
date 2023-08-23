include("../src/fba.jl")

organisms = [
    "iAF692", 
    # "iJR904", 
    # "iML1515", 
    "e_coli_core",
    # "iNF517",
    # "iSB619",
    # "iNJ661",
    # "iCN900",
    # "iAF1260",
    # "iEK1008",
    # "iJO1366",
    # "iMM904",
    # "iSDY_1059",
    # "iSFV_1184",
    # "iSF_1195",
    # "iS_1188",
    # "iSbBS512_1146"
]

for organism in organisms
    get_fba_data(organism, save_lp=false, json=true)
end 

# organism = "iML1515"
# get_fba_data(organism, save_lp=true, json=false)

# yeast model
# organisms = [
#     "Alloascoidea_hylecoeti",
#     "Ambrosiozyma_kashinagacola",
#     "Ambrosiozyma_monospora",
#     "Arthrobotrys_oligospora",
#     "Arxula_adeninivorans",
#     "Ascoidea_asiatica",
#     "Ascoidea_rubescens",
#     "Ashbya_aceri",
#     "Aspergillus_nidulans",
#     "Babjeviella_inositovora",
#     "Botrytis_cinerea"
# ]

# for organism in organisms
#     get_fba_data(organism, yeast=true)
# end 

