#using DelimitedFiles
using CSV
using COBREXA
using DataFrames

# organisms = readdlm("ecModel_names.txt", '\t', String, '\n')

# organism_data = []
# for organism in organisms
#     molecular_model = load_model("../molecular_models/ecModels/Classical/emodel_" * organism * "_classical.mat")
#     S = stoichiometry(molecular_model)
#     m, num_reactions = size(S)
#     push!(organism_data, (organism, m, num_reactions))
# end

# # num metabolites
# sort!(organism_data, by = organism_data -> organism_data[2])
# @show organism_data[1:10]
# min_m = organism_data[1:5]
# # sort!(organism_data, by = organism_data -> organism_data[2], rev=true)
# # @show organism_data[1:10]
# println("")

# # num reactions
# sort!(organism_data, by = organism_data -> organism_data[3])
# @show organism_data[1:10]
# min_r = organism_data[1:5]
# # sort!(organism_data, by = organism_data -> organism_data[3], rev=true)
# # @show organism_data[1:10]
# println("")

# sort!(organism_data, by = organism_data -> organism_data[2] + organism_data[3])
# @show organism_data[1:10]
# min_sum = organism_data[1:5]
# println("")

# @show unique(vcat(min_m, min_r, min_sum))
# smallest_organisms = [i[1] for i in unique(vcat(min_m, min_r, min_sum))]
# writedlm("ecModel_small_model_names.txt", smallest_organisms)

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
    "iZ_1308"
]

organism_data = []
for organism in organisms
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    push!(organism_data, (organism, m, num_reactions))
end

sort!(organism_data, by = organism_data -> organism_data[3], rev=false)
df = DataFrame(NamedTuple{(:organism, :metabolites, :reactions)}.(organism_data))

CSV.write("bigg_model_data.csv", df, append=false, writeheader=true)
