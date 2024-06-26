using DataFrames
using CSV
using Infiltrator

include("../src/cuts_decomposition.jl")

# organisms = [
#     "iAF692", 
#     "iJR904", 
#     "iML1515", 
#     "e_coli_core",
#     "iNF517",
#     "iSB619",
#     "iNJ661",
#     "iCN900",
#     "iAF1260",
#     "iEK1008",
#     "iJO1366",
#     "iMM904",
#     "iSDY_1059",
#     "iSFV_1184",
#     "iSF_1195",
#     "iS_1188",
#     "iSbBS512_1146"
# ]

# combinatorial_benders_data("iSbBS512_1146", time_limit=1800, fast=true)
# combinatorial_benders_data("iSFV_1184", time_limit=1800, fast=true)
# combinatorial_benders_data("iSB619", time_limit=1800, fast=true)
# combinatorial_benders_data("iAF692", time_limit=1800, fast=false, scip_tol=1e-5)

# for organism in organisms
#     type = "cb"
#     try 
#         combinatorial_benders_data(organism, time_limit=1800, fast=false)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end

#     type = "cb_fast"
#     try 
#         combinatorial_benders_data(organism, time_limit=1800, fast=true)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end
    
# #    type = "ch"
# #    try 
# #        constraint_handler_data(organism, time_limit=600)
# #    catch e 
# #        println(e)
# #        file = organism * "_" * type
# #        open(file * ".txt","a") do io
# #            println(io, e)
# #        end
# #    end
# end

# yeast model
organisms = [
    "Alloascoidea_hylecoeti",
#    "Ambrosiozyma_kashinagacola",
#    "Ambrosiozyma_monospora",
#    "Arthrobotrys_oligospora",
#    "Arxula_adeninivorans",
#    "Ascoidea_asiatica",
#    "Ascoidea_rubescens",
#    "Ashbya_aceri",
#    "Aspergillus_nidulans",
#    "Babjeviella_inositovora",
#    "Botrytis_cinerea"
]

# for organism in organisms
#     combinatorial_benders_data(organism, yeast=true, time_limit=1800*4, fast=false)
#     combinatorial_benders_data(organism, yeast=true, time_limit=1800*4, fast=true)
# end

for organism in organisms 
    type = "cb_fast_big_m"
    try 
        combinatorial_benders_data(organism, yeast=true, time_limit=60, fast=true, big_m=true)
    catch e 
        println(e)
        file = organism * "_" * type
        open(file * ".txt","a") do io
            println(io, e)
        end
    end
end
