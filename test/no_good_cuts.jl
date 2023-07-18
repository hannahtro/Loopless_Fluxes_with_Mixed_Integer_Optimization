using Test 
using DataFrames
using CSV

include("../src/cuts_decomposition.jl")
include("../src/constraint_handler.jl")

@testset "simple model" begin
    S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
    lb = [0,-10,-10,-10,0,0,0]
    ub = [20,30,30,30,20,10,10]
    m, num_reactions = size(S)
    @show m, num_reactions

    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    internal_rxn_idxs = [2,3,4,6,7]

    println("### No good cuts")
    # no good cuts
    objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S)
    @test thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
    println("--------------------------------------------------------")

    println("### combinatorial Benders")
    # combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    objective_value, objective_values, dual_bounds, solution, time, termination, iter = combinatorial_benders(model, internal_rxn_idxs, S, fast=false)
    @show objective_value, solution
    @test termination == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution[internal_rxn_idxs], S)
    @test feasible
    @test objective_values == sort(objective_values, rev=true)
    println("--------------------------------------------------------")

    println("### fast combinatorial Benders")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, time_fast, termination_fast, iter_fast = combinatorial_benders(model, internal_rxn_idxs, S, fast=true)
    @test termination_fast == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution[internal_rxn_idxs], S)
    @test feasible
    @test objective_values_fast == sort(objective_values_fast, rev=true)

    @test iter >= iter_fast
    # @test time >= time_fast
    @test isapprox(objective_value,objective_value_fast)
    @test solution[1:num_reactions] == solution_fast[1:num_reactions]
    @show objective_value_fast, solution_fast
    println("--------------------------------------------------------")

    println("### constraint handler")
    # test constraint handler        
    scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true)
    # print(scip_model)
    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars)
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch")
    MOI.optimize!(scip_model)
    primal_objective_value = MOI.get(scip_model, MOI.ObjectiveValue())
    @show primal_objective_value
    # @show MOI.get(scip_model, MOI.VariablePrimal(), [x,a])
    solution = MOI.get(scip_model, MOI.VariablePrimal(), flux_vars)
    @show solution
    bin_vals = MOI.get(scip_model, MOI.VariablePrimal(), bin_vars)
    @show bin_vals
    feasible = thermo_feasible(internal_rxn_idxs, solution[internal_rxn_idxs], S)
    @test feasible
    @assert bin_vals[1:length(internal_rxn_idxs)] + bin_vals[length(internal_rxn_idxs)+1:end] == ones(length(internal_rxn_idxs))
end

# # TODO: no good cuts approach does not terminate in 200 iterations: verify that solution is eventually found
# @testset "iAF692" begin
#     organism = "iAF692"
#     molecular_model = deserialize("../data/" * organism * ".js")
#     print_model(molecular_model, "organism")

#     S = stoichiometry(molecular_model)
#     # lb, ub = bounds(molecular_model)
#     internal_rxn_idxs = [
#         ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
#         !is_boundary(reaction_stoichiometry(molecular_model, rid))
#     ]

#     # model = build_fba_model(S, lb, ub, optimizer=SCIP.Optimizer)
#     model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     time_limit = 2
#     objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

#     try 
#         thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
#     catch 
#     else 
#         @test time >= time_limit
#     end

#     # combinatorial Benders'
#     model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     combinatorial_benders(model, internal_rxn_idxs, S, max_iter=5, fast=false)

#     # fast combinatorial Benders'
#     model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     combinatorial_benders(model, internal_rxn_idxs, S, max_iter=5, fast=true)
# end

# no_good_cuts_data("iAF692", time_limit=3600)

# println("--------------------------------------------------------")
# combinatorial_benders_data("iAF692", time_limit=1800, csv=true, fast=false, silent=false)
# #combinatorial_benders_data("iAF692", time_limit=1800, csv=true, fast=true, silent=false)
# combinatorial_benders_data("iJR904", time_limit=1800, csv=true, fast=false, silent=false)
# #combinatorial_benders_data("iJR904", time_limit=1800, csv=true, fast=true, silent=false)
# combinatorial_benders_data("iML1515", time_limit=1800, csv=true, fast=false, silent=false)
# #combinatorial_benders_data("iML1515", time_limit=1800, csv=true, fast=true, silent=false)
