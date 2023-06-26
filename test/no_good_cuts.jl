using Test 
using DataFrames
using CSV

include("../src/cuts_decomposition.jl")

# @testset "simple model" begin
#     S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
#     lb = [0,-10,-10,-10,0,0,0]
#     ub = [20,30,30,30,20,10,10]
#     m, num_reactions = size(S)
#     @show m, num_reactions

#     model = build_fba_model(S, lb, ub)
#     internal_rxn_idxs = [2,3,4,6,7]

#     # no good cuts
#     objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S)
#     @test thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

#     # combinatorial Benders'
#     model = build_fba_model(S, lb, ub)
#     objective_value, dual_bounds, solution, time, termination, iter = combinatorial_benders(model, internal_rxn_idxs, S, fast=false)
#     @test termination == MOI.OPTIMAL 
#     feasible = thermo_feasible(internal_rxn_idxs, solution[internal_rxn_idxs], S)
#     @test feasible
#     println("--------------------------------------------------------")

#     # fast combinatorial Benders'
#     model = build_fba_model(S, lb, ub)
#     objective_value_fast, dual_bounds_fast, solution_fast, time_fast, termination_fast, iter_fast = combinatorial_benders(model, internal_rxn_idxs, S, fast=true)
#     @test termination_fast == MOI.OPTIMAL 
#     feasible = thermo_feasible(internal_rxn_idxs, solution[internal_rxn_idxs], S)
#     @test feasible

#     @test iter >= iter_fast
#     @test time >= time_fast
#     @test isapprox(objective_value,objective_value_fast)
#     @test solution == solution_fast
# end

# # TODO: no good cuts approach does not terminate in 200 iterations: verify that solution is eventually found
# @testset "iAF692" begin
#     organism = "iAF692"
#     model = deserialize("../data/" * organism * ".js")
#     print_model(model, "organism")

#     S = stoichiometry(model)
#     lb, ub = bounds(model)
#     internal_rxn_idxs = [
#         ridx for (ridx, rid) in enumerate(variables(model)) if
#         !is_boundary(reaction_stoichiometry(model, rid))
#     ]

    # model = build_fba_model(S, lb, ub)

    # time_limit = 2
    # objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

    # try 
    #     thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
    # catch 
    # else 
    #     @test time >= time_limit
    # end

#     # combinatorial Benders'
#     model = build_fba_model(S, lb, ub)
#     combinatorial_benders(model, internal_rxn_idxs, S, max_iter=10, fast=false)

#     # fast combinatorial Benders'
#     model = build_fba_model(S, lb, ub)
#     combinatorial_benders(model, internal_rxn_idxs, S, max_iter=10, fast=true)
# end

# no_good_cuts_data("iAF692", time_limit=3600)

# println("--------------------------------------------------------")
# combinatorial_benders_data("iAF692", max_iter=10, time_limit=30, csv=false, fast=false)
combinatorial_benders_data("iAF692", max_iter=10, time_limit=30, csv=false, fast=true)