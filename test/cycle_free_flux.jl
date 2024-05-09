using Test 

include("../src/fba.jl")
include("../src/cycle_free_flux.jl")
include("../src/loopless_fba.jl")

println("============================================================")
println("CycleFreeFlux")
println("============================================================")

# @testset "simple model with 1 loop" begin
#     println("--------------------------------------------------------")
#     println("TEST SIMPLE MODEL (1 loop)")
#     println("--------------------------------------------------------")
#     S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
#     lb = [0,-30,-30,-30,0]
#     ub = [10,30,30,30,10]
#     internal_rxn_idxs = [2,3,4]
#     m, num_reactions = size(S)
#     @show m, num_reactions
#     println("-----------------------------------")

#     # infeasible loop
#     model = build_fba_model(S, lb, ub)
#     x = model[:x]
#     @objective(model, Max, x[2]+x[3]+x[4])
#     objective_fba, _, solution, _, _ = optimize_model(model)
#     @show objective_fba, solution

#     # test feasibility, filter non-zero fluxes, set binaries accordingly
#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     @test !feasible

#     # test CycleFreeFlux algorithm
#     # extract objective function 
#     objective_func = objective_function(model)
#     objective_func_vars = [i.index.value for i in objective_func.terms.keys]
#     objective_func_coeffs = objective_func.terms.vals

#     # example for loop in solution
#     objective_val, vars = cycle_free_flux(solution, objective_fba, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs)
#     @show objective_val, vars
#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(vars) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [vars[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     @show thermo_feasible
#     @test !thermo_feasible

#     # ll FBA
#     model = build_fba_model(S, lb, ub)
#     x = model[:x]
#     @objective(model, Max, x[2]+x[3]+x[4])
#     add_loopless_constraints(model, S, internal_rxn_idxs)
#     objective_ll_fba, _, solution, _, _ = optimize_model(model)
#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     @show objective_ll_fba, solution
#     @test feasible
# end

# @testset "simple model" begin
#     println("--------------------------------------------------------")
#     println("TEST SIMPLE MODEL (2 loops)")
#     println("--------------------------------------------------------")
#     S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
#     lb = [0,-10,-10,-10,0,0,0]
#     ub = [20,30,30,30,20,10,10]
#     m, num_reactions = size(S)
#     @show m, num_reactions

#     # solve fba
#     optimization_model = build_fba_model(S, lb, ub, set_objective=true)
#     internal_rxn_idxs = [2,3,4,6,7]
#     objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(optimization_model, print_objective=true)
#     # solution contains loop
#     @show vars_fba

#     # test ll-ness through internal reactions
#     internal_flux_directions = [vars_fba[idx] >= 1e-5 ? 1 : 0 for idx in internal_rxn_idxs]
#     thermo_feasible = thermo_feasible_mu(internal_rxn_idxs, internal_flux_directions, S)
#     @show thermo_feasible
#     @test !thermo_feasible

#     # test ll-ness through internal, non-zero reactions
#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(vars_fba) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [vars_fba[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     @show thermo_feasible
#     @test !thermo_feasible

#     # test CycleFreeFlux algorithm
#     # extract objective function 
#     objective_func = objective_function(optimization_model)
#     objective_func_vars = [i.index.value for i in objective_func.terms.keys]
#     objective_func_coeffs = objective_func.terms.vals

#     # example for loop in solution
#     @warn "solution of algorithm is not loopless"
#     objective_val, vars = cycle_free_flux(vars_fba, objective_fba, S, internal_rxn_idxs, objective_func_vars, objective_func_coeffs)
#     @show vars
#     non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(vars) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
#     non_zero_flux_directions = [vars[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
#     thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
#     @show thermo_feasible
#     @test !thermo_feasible
# end 

@testset "FVA of simple model with 1 loop" begin
    println("--------------------------------------------------------")
    println("TEST SIMPLE MODEL (1 loop)")
    println("--------------------------------------------------------")
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-30,-30,-30,0]
    ub = [10,30,30,30,10]
    internal_rxn_idxs = [2,3,4]
    m, num_reactions = size(S)
    @show m, num_reactions
    println("-----------------------------------")
    cycles = cycle_free_fva(S, internal_rxn_idxs, lb, ub, optimizer=HiGHS.Optimizer)
    @show cycles
    @test length(cycles) == 1
    thermo_feasible = thermo_feasible_mu(cycles[1][1], cycles[1][2], S)
    @test !thermo_feasible
end

@testset "FVA of e_coli_core" begin
    println("--------------------------------------------------------")
    println("TEST e coli core")
    println("--------------------------------------------------------")
    organism = "e_coli_core"
    molecular_model = load_model("../molecular_models/" * organism * ".json")

    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)   
    println("-----------------------------------")
    # cycles = cycle_free_fva(S, internal_rxn_idxs, lb, ub, optimizer=HiGHS.Optimizer)
    # # @show cycles
    # @test length(cycles) == 50
    # cycles_to_block = [cycle for cycle in cycles if !thermo_feasible_mu(cycle[1], cycle[2], S)]
    # @show length(cycles_to_block)

    loopless_fba_blocked_data(organism; time_limit=10,shortest_cycles=false, block_limit=100, type="loopless_fba_blocked", nullspace_formulation=false, reduced=false, json=false, cff=true, cff_time_limit=10)
end