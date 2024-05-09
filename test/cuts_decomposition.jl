using Test 
using DataFrames
using CSV
using Infiltrator
using GLPK

include("../src/cuts_decomposition.jl")
include("../src/optimization_model.jl")

println("============================================================")
println("COMBINATORIAL BENDERS")
println("============================================================")

@testset "simple model" begin
    println("--------------------------------------------------------")
    println("TEST SIMPLE MODEL")
    println("--------------------------------------------------------")
    S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
    lb = [0,-10,-10,-20,0,0,0]
    ub = [20,30,30,30,20,10,10]
    m, num_reactions = size(S)
    internal_rxn_idxs = [2,3,4,6,7]
    @show m, num_reactions

    model = build_fba_model(S, lb, ub, set_objective=true)
    primal_objective_value, dual_objective_value, solution, time, status = optimize_model(model, print_objective=true, silent=true, mute=false)
    @show solution 

    println("--------------------------------------------------------")
    
    # println("no good cuts")
    # # no good cuts
    # objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S)
    # solution_flux = solution[1:num_reactions]
    # @show solution_flux
    # solution_direction = solution[num_reactions+1:end]
    # @show solution_direction
    # nonzero_flux_idxs = [idx for (idx,i) in enumerate(solution_flux[internal_rxn_idxs]) if !isapprox(0, i, atol=0.001)]
    # @show nonzero_flux_idxs
    # @test thermo_feasible_mu(internal_rxn_idxs[nonzero_flux_idxs], solution_direction[nonzero_flux_idxs], S)
    # println("--------------------------------------------------------")

    println("combinatorial Benders with no good cuts")
    # combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    objective_value_cb, objective_values_cb, dual_bounds_cb, solution_cb, x_cb, a_cb, G_cb, _, time_cb, termination_cb, iter_cb, cuts_cb, _, _, _ = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=false, indicator=true)
    @show objective_value_cb, solution_cb
    @test termination_cb == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_cb[internal_rxn_idxs], S)
    @test feasible
    @show objective_value_cb
    @test round.(objective_values_cb, digits=5) == round.(sort(objective_values_cb, rev=true), digits=5)
    # @test isapprox(objective_value, objective_value_cb, atol=0.0001)
    println("--------------------------------------------------------")

    println("fast combinatorial Benders with big M")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_big_m, objective_values_big_m, dual_bounds_big_m, solution_big_m, _, _, _, _, time_big_m, termination_big_m, iter_big_m = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, big_m=true)
    @test termination_big_m == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_big_m[internal_rxn_idxs], S)
    @test feasible
    @test isapprox(objective_value_big_m, objective_value_cb)
    # @test isapprox(solution_big_m[1:num_reactions], solution_fast[1:num_reactions], atol=0.00001)
    @show objective_value_big_m, objective_value_cb
    println("--------------------------------------------------------")

    println("fast combinatorial Benders with indicator")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, x_fast, a_fast, G_fast, _, time_fast, termination_fast, iter_fast, cuts_fast, _, _, _ = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, indicator=true)
    @test termination_fast == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_fast[internal_rxn_idxs], S)
    @test feasible
    @test round.(objective_values_fast, digits=5) == round.(sort(objective_values_fast, rev=true), digits=5)

    @test iter_cb >= iter_fast
    # @test time >= time_fast
    @test isapprox(objective_value_cb, objective_value_fast)
    @test isapprox(solution_cb[1:num_reactions], solution_fast[1:num_reactions], atol=0.00001)
    @show objective_value_fast, solution_fast
    println("--------------------------------------------------------")

    println("fast combinatorial Benders with indicator AND big M")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, x_fast, a_fast, G_fast, _, time_fast, termination_fast, iter_fast, cuts_fast, _, _, _ = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, indicator=true, big_m=true)
    @test termination_fast == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_fast[internal_rxn_idxs], S)
    @test feasible
    @test round.(objective_values_fast, digits=5) == round.(sort(objective_values_fast, rev=true), digits=5)

    @test iter_cb >= iter_fast
    # @test time >= time_fast
    @test isapprox(objective_value_cb, objective_value_fast)
    @test isapprox(solution_cb[1:num_reactions], solution_fast[1:num_reactions], atol=0.00001)
    @show objective_value_fast, solution_fast
    println("--------------------------------------------------------")

    println("fast combinatorial Benders with multiple MISs")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_multiple_mis, objective_values_multiple_mis, dual_bounds_multiple_mis, solution_multiple_mis, x_multiple_mis, a_multiple_mis, G_multiple_mis, _, time_multiple_mis, termination_multiple_mis, iter_multiple_mis, cuts_multiple_mis, _, _, _ = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, multiple_mis=10, indicator=true)
    @test termination_multiple_mis == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_multiple_mis[internal_rxn_idxs], S)
    @test feasible
    @test round.(objective_values_multiple_mis, digits=5) == round.(sort(objective_values_multiple_mis, rev=true), digits=5)

    @test iter_fast >= iter_multiple_mis
    # @test time >= time_fast
    @test isapprox(objective_value_multiple_mis, objective_value_fast)
    @test isapprox(solution_cb[1:num_reactions], solution_multiple_mis[1:num_reactions], atol=0.00001)
    @show objective_value_multiple_mis, solution_multiple_mis
    println("--------------------------------------------------------")

    println("test loopless violation")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # run 2 iterations where solution is not loopless
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, x_fast, a_fast, G_fast, _, time_fast, termination_fast, iter_fast, cuts_fast, _, _, _  = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, max_iter=2, indicator=true)
    flux = solution_fast[1:num_reactions]
    flux_directions = solution_fast[num_reactions+1:num_reactions+length(internal_rxn_idxs)]
    @show flux, flux_directions
    diff, solution = check_loopless_violation(flux, flux_directions, S, internal_rxn_idxs)
    @show diff, solution
    
    # solve till optimal solution found
    model = build_fba_model(S, lb, ub, set_objective=true)
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, x_fast, a_fast, G_fast, _, time_fast, termination_fast, iter_fast, cuts_fast, _, _, _  = combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, fast=true, indicator=true)
    flux = solution_fast[1:num_reactions]
    flux_directions = solution_fast[num_reactions+1:num_reactions+length(internal_rxn_idxs)]
    @show flux, flux_directions
    diff, solution = check_loopless_violation(flux, flux_directions, S, internal_rxn_idxs)
    @show diff, solution
    println("--------------------------------------------------------")

#     # println("constraint handler")
#     # # test constraint handler    
#     # scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, silent=true)
#     # # print(scip_model)
#     # ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [])
#     # SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch", enforce_priority=-99999999, check_priority=-99999999)

#     # MOI.optimize!(scip_model)
#     # @test MOI.get(scip_model, MOI.TerminationStatus()) == MOI.OPTIMAL
#     # objective_value_ch = MOI.get(scip_model, MOI.ObjectiveValue())
#     # solution_ch = [MOI.get(ch.o, MOI.VariablePrimal(1), MOI.VariableIndex(i)) for i in 1:length(internal_rxn_idxs) + num_reactions]
#     # @show objective_value_ch
#     # @show solution_ch 
#     # feasible = thermo_feasible(internal_rxn_idxs, solution_ch[internal_rxn_idxs], S)
#     # @test feasible
#     # @test isapprox(objective_value_ch,objective_value_fast)
    
#     # print SCIP solution
#     # SCIP.SCIPprintSol(ch.o, SCIP.SCIPgetBestSol(ch.o), C_NULL, SCIP.TRUE)
end

# @testset "iAF692" begin
#     println("")
#     println("--------------------------------------------------------")
#     println("TEST iAF692")
#     println("--------------------------------------------------------")
#     organism = "iAF692"
#     molecular_model = deserialize("../molecular_models/" * organism * ".js")
#     # print_model(molecular_model, "organism")

#     S = stoichiometry(molecular_model)
#     lb, ub = bounds(molecular_model)
#     internal_rxn_idxs = [
#         ridx for (ridx, rid) in enumerate(reactions(molecular_model)) if
#         !is_boundary(reaction_stoichiometry(molecular_model, rid))
#     ]

#     # # model = build_fba_model(S, lb, ub, optimizer=SCIP.Optimizer)
#     # model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     # time_limit = 2
#     # objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

#     # try 
#     #     thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
#     # catch 
#     # else 
#     #     @test time >= time_limit
#     # end

#     println("combinatorial Benders")
#     # combinatorial Benders'
#     model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, max_iter=5, fast=false)
#     println("--------------------------------------------------------")

#     println("fast combinatorial Benders")
#     # fast combinatorial Benders'
#     model = make_optimization_model(molecular_model, SCIP.Optimizer)
#     combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, max_iter=5, fast=true)
#     println("--------------------------------------------------------")
    
#     println("fast combinatorial Benders with multiple MISs")
#     # fast combinatorial Benders'
#     model = build_fba_model(S, lb, ub, set_objective=true)
#     # print(model)
#     combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, max_iter=5, fast=true, multiple_mis=10)
    
#     println("fast combinatorial Benders with big M")
#     # fast combinatorial Benders'
#     model = build_fba_model(S, lb, ub, set_objective=true)
#     combinatorial_benders(model, internal_rxn_idxs, S, lb, ub, max_iter=5, fast=true, big_m=true)

#     # println("constraint handler")
#     # # test constraint handler
#     # # extract objective
#     # objective_func = objective_function(model)
#     # objective_func_vars = [i.index for i in objective_func.terms.keys]
           
#     # scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, time_limit=10, objective_func_vars=objective_func_vars, objective_func_coeffs=objective_func.terms.vals)
#     # # print(scip_model)
#     # ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [])
#     # SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch", enforce_priority=-7000000, check_priority=-7000000)
#     # MOI.optimize!(scip_model)
#     # # @test MOI.get(scip_model, MOI.TerminationStatus()) == MOI.TIME_LIMIT
# end
