using DelimitedFiles

organisms = [
#    "iAF692", 
    "iJR904", 
#    "iML1515", 
#    "e_coli_core",
#    "iNF517",
#    "iSB619",
    "iNJ661",
#    "iCN900",
#    "iAF1260",
#    "iEK1008",
#    "iJO1366",
#    "iMM904",
#    "iSDY_1059",
#    "iSFV_1184",
#    "iSF_1195",
#    "iS_1188",
#    "iSbBS512_1146"
]

#organisms = readdlm("../molecular_models/ecModel_small_model_names.txt", '\t', String, '\n')
# organisms = [
#     "Hanseniaspora_uvarum",
#     "yHMPu5000035696_Hanseniaspora_singularis",
#     "yHMPu5000034963_Hanseniaspora_clermontiae",
#     "yHMPu5000035695_Hanseniaspora_pseudoguilliermondii",
#     "yHMPu5000035684_Kloeckera_hatyaiensis",
#     "Eremothecium_sinecaudum",
#     "yHMPu5000035659_Saturnispora_dispora",
#     "Tortispora_caseinolytica",
#     "Starmerella_bombicola_JCM9596",
#     "Eremothecium_gossypii",
#     "Ashbya_aceri"
# ]

time_limit = 1800
fast = true
json = true
yeast = false

# for organism in organisms
#     @show organism
#     run(`sbatch -A optimi batch.sh $organism $time_limit $fast $json $yeast`) # CB
#     # run(`sbatch -A optimi batch.sh $organism $time_limit $json $yeast`) # ll FBA
#     # run(`sh batch.sh $organism $time_limit $fast $json $yeast`)
# end

mis_numbers = [5, 10, 20, 30]
for organism in organisms
    for mis in mis_numbers
        @show organism, mis
        run(`sbatch -A optimi batch.sh $organism $time_limit $fast $json $yeast $mis`) # CB
        # run(`sbatch -A optimi batch.sh $organism $time_limit $json $yeast`) # ll FBA
        # run(`sh batch.sh $organism $time_limit $fast $json $yeast`)
    end
end

