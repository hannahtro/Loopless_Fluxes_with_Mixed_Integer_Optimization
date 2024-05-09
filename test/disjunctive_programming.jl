using Test 

include("../src/disjunctive_programming.jl")
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
    internal_rxn_idxs = [2,3,4,6,7]

    println("ll FBA")
    optimization_model = build_fba_model(S, lb, ub, set_objective=true)
    add_loopless_constraints_mu(optimization_model, S, internal_rxn_idxs)
    print(optimization_model)
    objective_value , _, solution, _, status = optimize_model(optimization_model)
    @test status == MOI.OPTIMAL
    @show objective_value

    println("DP")
    model = build_disjunctive_fba_model(S, lb, ub, set_objective=true)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))
    @show max_flux_bound
    add_disjunctive_loopless_constraints(model, S, internal_rxn_idxs, max_flux_bound)
    print(model)

    println("")
    println(" Big M ")
    # print(model)
    set_attribute(model, MOI.Silent(), true)
    primal_objective_value, _, solution, time_taken, status = optimize_gdp_model(model)
    @test status == MOI.OPTIMAL
    @test isapprox(primal_objective_value, 80, atol=1e-6) 
    @show primal_objective_value, solution
    print(model)
    flux = value.(model[:x])
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test feasible

    println("")
    println(" Indicator ")
    primal_objective_value, _, solution, time_taken, status = optimize_gdp_model(model, gdp_method="Indicator")
    @test status == MOI.OPTIMAL
    @test isapprox(primal_objective_value, 80, atol=1e-6) 
    @show primal_objective_value, solution
    # print(model)
    flux = value.(model[:x])
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test feasible

    println("")
    println(" Hull ")
    primal_objective_value, _, solution, time_taken, status = optimize_gdp_model(model, gdp_method="Hull")
    @test status == MOI.OPTIMAL
    @test isapprox(primal_objective_value, 80, atol=1e-6) 
    @show primal_objective_value, solution
    print(model)
    flux = value.(model[:x])
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test feasible
end