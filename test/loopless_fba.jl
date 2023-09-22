using Test 

include("../src/loopless_fba.jl")

println("============================================================")
println("LOOPLESS FBA")
println("============================================================")

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

    # check thermodynamic feasibility of solution through non zero flux
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test feasible
end

@testset "iAF692" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST iAF692")
    println("--------------------------------------------------------")
    organism = "iAF692"
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, "organism")

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    # ll-FBA
    primal_objective_value, solution, status = loopless_fba_data(organism, time_limit=1800, json=false)

    # test feasibility, filter non-zero fluxes, set binaries accordingly
    solution = solution[1:num_reactions]
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test feasible
end 
