using Test 

include("../src/fba.jl")

println("============================================================")
println("FBA")
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

    # solve fba
    optimization_model = build_fba_model(S, lb, ub, set_objective=true)
    internal_rxn_idxs = [2,3,4,6,7]
    objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(optimization_model, print_objective=true)
    # solution contains loop
    @show vars_fba

    # test ll-ness through internal reactions
    internal_flux_directions = [vars_fba[idx] >= 1e-5 ? 1 : 0 for idx in internal_rxn_idxs]
    thermo_feasible = thermo_feasible_mu(internal_rxn_idxs, internal_flux_directions, S)
    @show thermo_feasible
    @test !thermo_feasible

    # test ll-ness through internal, non-zero reactions
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(vars_fba) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [vars_fba[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @show thermo_feasible
    @test !thermo_feasible
end 

@testset "iAF692" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST iAF692")
    println("--------------------------------------------------------")
    organism = "iAF692"
    molecular_model = deserialize("../molecular_models/" * organism * ".js")
    # print_model(molecular_model, "organism")

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)

    # test ll-ness through internal reactions
    internal_flux_directions = [vars_fba[idx] >= 1e-5 ? 1 : 0 for idx in internal_rxn_idxs]
    feasible = thermo_feasible_mu(internal_rxn_idxs, internal_flux_directions, S)
    @show feasible
    @test !feasible

    # test feasibility, filter non-zero fluxes, set binaries accordingly
    solution = vars_fba[1:num_reactions]
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @test !feasible
end 