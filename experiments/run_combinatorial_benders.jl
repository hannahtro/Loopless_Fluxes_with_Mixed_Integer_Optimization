using DelimitedFiles

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
    "iSbBS512_1146",
    "RECON1",
    "Recon3D",
    "STM_v1_0",
    "iAB_RBC_283",
    "iAPECO1_1312",
    "iECB_1328",
    "iETEC_1333",
    "iHN637",
    "iIS312_Amastigote",
    "iJB785",
    "iJN746",
    "iLB1027_lipid",
    "iMM1415",
    "iND750",
    "iRC1080",
    "iSFxv_1172",
    "iSynCJ816",
    "iYO844",
    "iYS1720",
    "iZ_1308",
]

"""organisms = [
    "Hanseniaspora_uvarum",
    "yHMPu5000035696_Hanseniaspora_singularis",   
    "yHMPu5000034963_Hanseniaspora_clermontiae",
    "yHMPu5000035695_Hanseniaspora_pseudoguilliermondii",
    "yHMPu5000035684_Kloeckera_hatyaiensis",
    "Eremothecium_sinecaudum",
    "yHMPu5000035659_Saturnispora_dispora",
    "Tortispora_caseinolytica",
    "Starmerella_bombicola_JCM9596",
    "Eremothecium_gossypii",
    "Ashbya_aceri"
]
"""

time_limit = 1800
fast = true
json = true
yeast = false

mis_numbers = [3, 4, 5]# , 2] #[0, 5, 10, 20, 30]
#densities = [5, 10, 15, 20]
max_density = 1000
max_cuts = [2.0]
for organism in organisms
    for mis in mis_numbers
#        for max_density in densities
		for max_cut in max_cuts
	            @show organism, mis, max_density, max_cut
	            run(`sbatch -A optimi batch_cb.sh $organism $time_limit $fast $json $yeast $mis $max_density $max_cut`) # CB
		end
#        end
    end
end

# run(`sbatch -A optimi batch_cb.sh $"yHMPu5000034963_Hanseniaspora_clermontiae" $time_limit $fast $json $yeast $0`)

# run(`sbatch -A optimi batch_cb.sh $"Hanseniaspora_uvarum" $time_limit $fast $json $yeast $0.1`)
# run(`sbatch -A optimi batch_cb.sh $"Tortispora_caseinolytica" $time_limit $fast $json $yeast $0.1`)

# run(`sbatch -A optimi batch_cb.sh $"Hanseniaspora_uvarum" $time_limit $fast $json $yeast $0.2`)
