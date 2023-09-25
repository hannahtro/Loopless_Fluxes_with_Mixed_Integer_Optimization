using Test 

include("../src/cobrexa.jl")

println("============================================================")
println("COBREXA")
println("============================================================")

@testset "ecoli core" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST ecoli core")
    println("--------------------------------------------------------")
    # load data 
    organism = "e_coli_core"
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    S = stoichiometry(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)

    # FBA 
    primal_objective_value, solution, status = cobrexa_fba_data(organism, time_limit=1800, json=false)
    @test status == MOI.OPTIMAL

    # assert that there is no loop in solution
    internal_rxn_idxs = []
    flux_directions = solution[internal_rxn_idxs]
    @test thermo_feasible(internal_rxn_idxs, flux_directions, S)
    @test thermo_feasible_mu(internal_rxn_idxs, flux_directions, S)
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
    @show size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)

    # FBA
    primal_objective_value, solution, status = cobrexa_fba_data(organism, time_limit=1800, json=false)

    # test feasibility, filter non-zero fluxes, set binaries accordingly
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for (idx,val) in enumerate(non_zero_flux_indices)]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @show feasible
    @test !feasible

    # ll-FBA
    primal_objective_value, solution, status = cobrexa_loopless_fba_data(organism, time_limit=1800, json=false)

    # test feasibility, filter non-zero fluxes, set binaries accordingly
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for (idx,val) in enumerate(non_zero_flux_indices)]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @show feasible
    @test thermo_feasible(non_zero_flux_indices, non_zero_flux_directions, S)
    @test thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
end 