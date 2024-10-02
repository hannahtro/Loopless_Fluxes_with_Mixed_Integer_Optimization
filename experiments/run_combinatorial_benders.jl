using DelimitedFiles

organisms = [
    # "iAF692", 
    # "iJR904", 
    # "iML1515", 
    # "e_coli_core",
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
    # "iSbBS512_1146",
    # "RECON1",
    # "Recon3D",
    # "STM_v1_0",
    # "iAB_RBC_283",
    # "iAPECO1_1312",
    # "iECB_1328",
    # "iETEC_1333",
    # "iHN637",
    # "iIS312_Amastigote",
    # "iJB785",
    # "iJN746",
    # "iLB1027_lipid",
    # "iMM1415",
    # "iND750",
    # "iRC1080",
    # "iSFxv_1172",
    # "iSynCJ816",
    # "iYO844",
    # "iYS1720",
    # "iZ_1308",
     "iAF1260b",
     "iAF987",    
     "iAM_Pb448",
     "iAM_Pc455",
     "iAM_Pf480",
     "iAM_Pk459",
     "iAM_Pv461",
     "iAT_PLT_636",
     "iB21_1397",
     "iBWG_1329",
     "ic_1306",
     "iCHOv1",
     "iCHOv1_DG44",
     "iCN718",
     "iE2348C_1286",
     "iEC042_1314",
     "iEC1344_C",
     "iEC1349_Crooks",
     "iEC1356_Bl21DE3",
     "iEC1364_W",
     "iEC1368_DH5a",
     "iEC1372_W3110",
     "iEC55989_1330",
     "iECABU_c1320",
     "iECBD_1354",
     "iECD_1391",
     "iECDH10B_1368",
     "iEcDH1_1363",
     "iECDH1ME8569_1439",
     "iEcE24377_1341",
     "iECED1_1282",
     "iECH74115_1262",
     "iEcHS_1320",
     "iECIAI1_1343",
     "iECIAI39_1322",
     "iECNA114_1301",
     "iECO103_1326",
     "iECO111_1330",
     "iECO26_1355",
     "iECOK1_1307",
     "iEcolC_1368",
     "iECP_1309",
     "iECs_1301",
     "iECS88_1305",
     "iECSE_1348",
     "iECSF_1327",
     "iEcSMS35_1347",
     "iECSP_1301",
     "iECUMN_1333",
     "iECW_1372",
     "iEKO11_1354",
     "iG2583_1286",
     "iIS312",
     "iIS312_Epimastigote",
     "iIS312_Trypomastigote",
     "iIT341",
     "iJN1463",
     "iJN678",
     "iLF82_1304",
     "iLJ478",
     "iNRG857_1313",
     "iPC815",
     "iSBO_1134",
     "iSSON_1240",
     "iUMN146_1321",
     "iUMNK88_1353",
     "iUTI89_1310",
     "iWFL_1372",
     "iY75_1357",
     "iYL1228",
     "iYS854"
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

mis_numbers = [0.5]#1, 2, 5] #3, 4, 5]# , 2] #[0, 5, 10, 20, 30]
densities = [5, 10, 15, 20]
#max_density = 1000
# max_cuts = [0.2, 0.3, 0.4] #0.5]#0.5, 1.0, 1.5]
max_cut = 0.2
for organism in organisms
    for mis in mis_numbers
        for max_density in densities
	# for max_cut in max_cuts
	            @show organism, mis, max_density, max_cut
	            run(`sbatch -A optimi batch_cb.sh $organism $time_limit $fast $json $yeast $mis $max_density $max_cut`) # CB
		# end
        end
    end
end

# run(`sbatch -A optimi batch_cb.sh $"yHMPu5000034963_Hanseniaspora_clermontiae" $time_limit $fast $json $yeast $0`)

# run(`sbatch -A optimi batch_cb.sh $"Hanseniaspora_uvarum" $time_limit $fast $json $yeast $0.1`)
# run(`sbatch -A optimi batch_cb.sh $"Tortispora_caseinolytica" $time_limit $fast $json $yeast $0.1`)

# run(`sbatch -A optimi batch_cb.sh $"Hanseniaspora_uvarum" $time_limit $fast $json $yeast $0.2`)
