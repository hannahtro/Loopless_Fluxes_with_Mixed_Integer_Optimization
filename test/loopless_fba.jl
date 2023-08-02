using Test 

include("../src/loopless_fba.jl")

@testset "simple model" begin
    println("--------------------------------------------------------")
    println("TEST SIMPLE MODEL")
    println("--------------------------------------------------------")
    S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
    lb = [0,-10,-10,-10,0,0,0]
    ub = [20,30,30,30,20,10,10]
    m, num_reactions = size(S)
    @show m, num_reactions

    # test thermodynamic feasible fba with nullspace formulation
    optimization_model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    internal_rxn_idxs = [2,3,4,6,7]

    add_loopless_constraints(optimization_model, S, internal_rxn_idxs)
    objective_value , _, solution, _, status = optimize_model(optimization_model)
    @test status == MOI.OPTIMAL
    @show objective_value

    # test thermodynamic feasible fba without nullspace formulation
    optimization_model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    internal_rxn_idxs = [2,3,4,6,7]

    add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
    objective_value , _, solution, _, status = optimize_model(optimization_model)
    @test status == MOI.OPTIMAL
    @show objective_value

    # check thermodynamic feasibility of solution
    flux_directions = solution[internal_rxn_idxs]
    @test thermo_feasible(internal_rxn_idxs, flux_directions, S)
    @test thermo_feasible_mu(internal_rxn_idxs, flux_directions, S)
end

# @show VERSION

# organism = "iAF692"
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=1000, block_limit=20, same_objective=false, nullspace_formulation=true, reduced=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=1000, block_limit=20, same_objective=false, nullspace_formulation=false, reduced=false)

# # iAF692_loopless_fba_blocked_100_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# # iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# # iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# # iAF692_loopless_fba_blocked_50_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)

# iAF692_loopless_fba_1800
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)
# # iAF692_loopless_fba_mu_1800
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)

# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# # iAF692_loopless_fba_blocked_shortest_cycles_100_1800_50
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


# organism = "iJR904"
# # iJR904_loopless_fba_blocked_100_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_shortest_cycles_100_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, shortest_cycles=true, block_limit=50, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=10000, same_objective=false, block_limit=50, nullspace_formulation=false)

# # # iJR904_loopless_fba_1800
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)

# # iJR904_loopless_indicator_fba_mu_1800
# loopless_indicator_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# # iJR904_loopless_indicator_fba_blocked_mu_1800_10_same_objective
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=50, nullspace_formulation=false)
# loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=100, nullspace_formulation=false)

# # iJR904_loopless_fba_blocked_1800_50
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_100
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, same_objective=false, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_200
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, same_objective=false, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_500
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, same_objective=false, nullspace_formulation=false)

# # iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50 TODO: same file twice? 
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, same_objective=false, vector_formulation=true, shortest_cycles=true, nullspace_formulation=false)

# # iJR904_loopless_fba_blocked_1800_50_same_objective
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=50, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_100_same_objective
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=100, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_200_same_objective
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=200, nullspace_formulation=false)
# # iJR904_loopless_fba_blocked_100_1800_500_same_objective
# loopless_fba_blocked_data(organism, time_limit=1800, ceiling=500, nullspace_formulation=false)


# organism = "iML1515"
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=false)
# loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)

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
  
# organisms = ["iAF692", "iJR904", "iML1515", "e_coli_core", "iNF517", "iSB619", "iNJ661", "iCN900"]

# for organism in organisms
#     type = "thermo_feasible_fba"
#     try 
#         loopless_fba_data(organism, time_limit=600)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end

#     type = "thermo_feasible_fba_nullspace"
#     try 
#         loopless_fba_data(organism, time_limit=600, nullspace_formulation=true)
#     catch e 
#         println(e)
#         file = organism * "_" * type
#         open(file * ".txt","a") do io
#             println(io, e)
#         end
#     end
# end

organism = "iAF692"
# loopless_fba_data(organism, time_limit=600, nullspace_formulation=true)
loopless_fba_bilinear_data(organism; time_limit=600, silent=false, type="loopless_bilinear_fba", csv=true)

