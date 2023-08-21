include("../src/loopless_fba.jl")

# @show VERSION

# organism = "e_coli_core"
# # loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true, json=false)

# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)
  
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

# for organism in organisms
#     # loopless_relaxed_fba_data(organism, time_limit=1, nullspace_formulation=false, csv=false, save_lp=true)
#     # loopless_relaxed_fba_data(organism, time_limit=1, nullspace_formulation=true, csv=false, save_lp=true)
#     type = "thermo_feasible_fba"
#     try 
#         loopless_fba_data(organism, time_limit=1800)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end

#     type = "thermo_feasible_fba_nullspace"
#     try 
#         loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end
# end

# organism = "iAF692"
# loopless_fba_data(organism, time_limit=30, nullspace_formulation=false, json=true)

# # yeast model
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
#     loopless_fba_data(organism, yeast=true, time_limit=1800*4)
# end 


