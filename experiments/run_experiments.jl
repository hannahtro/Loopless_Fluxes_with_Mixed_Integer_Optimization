using DelimitedFiles

# organisms = [
#     "iAF692", 
#     # "iJR904", 
#     # "iML1515", 
#     # "e_coli_core",
#     # "iNF517",
#     # "iSB619",
#     # "iNJ661",
#     # "iCN900",
#     # "iAF1260",
#     # "iEK1008",
#     # "iJO1366",
#     # "iMM904",
#     # "iSDY_1059",
#     # "iSFV_1184",
#     # "iSF_1195",
#     # "iS_1188",
#     # "iSbBS512_1146"
# ]

organisms = readdlm("../molecular_models/ecModel_small_model_names.txt", '\t', String, '\n')

time_limit = 3600*10
fast = true
json = true
yeast = true

for organism in organisms
    @show organism
    # run(`sbatch -A optimi batch.sh $organism $time_limit $fast $json $yeast`) # CB
    run(`sbatch -A optimi batch.sh $organism $time_limit $json $yeast`) # ll FBA
    # run(`sh batch.sh $organism $time_limit $fast $json $yeast`)
end