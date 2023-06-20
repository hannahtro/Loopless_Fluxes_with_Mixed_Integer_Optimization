using Test 
using DataFrames
using CSV

include("../src/cuts_decomposition.jl")

@testset "simple model" begin
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,30,30,30,10]
    m, num_reactions = size(S)
    @show m, num_reactions

    model = build_fba_model(S, lb, ub)
    internal_rxn_idxs = [2,3,4]

    # no good cuts
    objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S)
    @test thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

    # combinatorial Benders'
    model = build_fba_model(S, lb, ub)
    combinatorial_benders(model, internal_rxn_idxs, S)
end

# TODO: does not terminate in 200 iterations: verify that solution is eventually found
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

#     model = build_fba_model(S, lb, ub)

#     time_limit = 2
#     objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

#     try 
#         thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
#     catch 
#     else 
#         @test time >= time_limit
#     end

#     # combinatorial Benders'
#     model = build_fba_model(S, lb, ub)
#     combinatorial_benders(model, internal_rxn_idxs, S)
# end

# no_good_cuts_data("iAF692", time_limit=3600)

