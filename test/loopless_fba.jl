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
