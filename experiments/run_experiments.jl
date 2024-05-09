organisms = [
    "iAF692", 
    "iJR904", 
    "iML1515", 
    "e_coli_core",
    "iNF517",
    # "iSB619",
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

time_limit = 1800

seeds = 1:1#4
for organism in organisms
    for seed in seeds
        @show organism
        run(`sbatch -A optimi batch.sh $organism $seed`)
    end
end
