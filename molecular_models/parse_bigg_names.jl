organisms = []
open("bigg_models.txt") do io
    while !eof(io)
        line = readline(io)
        println(line)
        # print(split(line, "\t")[1])
        # println(" ")
        push!(organisms, first.(split.(line, "\t")))
    end
end


outfile = "bigg_model_names.txt"
f = open(outfile, "w")

used_organisms = [
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
    "iZ_1308"
]

for i in organisms
    if !(string(i) in used_organisms)
        println(i)
	    println(f, i)
    end
end
close(f) 