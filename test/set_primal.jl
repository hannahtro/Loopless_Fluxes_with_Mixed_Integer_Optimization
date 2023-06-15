
using Test 

include("../src/set_primal.jl")


@testset "simple model with feasible loop" begin
    S = [[1,0,0,0,0,0,0,0,0] [0,0,0,0,0,1,0,0,0] [0,1,0,0,0,-1,0,0,0] [0,-1,0,0,0,0,1,0,0] [0,0,0,0,0,0,-1,0,0] [0,0,0,0,-1,0,0,0,0] [-1,-1,1,1,0,0,0,0,0] [0,1,0,-1,1,0,0,0,0] [0,0,-1,0,0,0,0,0,0] [0,0,0,0,0,0,0,1,0] [0,0,0,1,0,0,0,-1,0] [0,0,0,-1,0,0,0,0,1] [0,0,0,0,0,0,0,0,-1]]
    m, num_reactions = size(S)
    lb = [0,0,0,0,0,0,0,0,0,0,0,0,0]
    ub = [10,10,10,10,10,10,10,10,10,10,10,10,10]
    model = build_model(S, lb, ub)
    internal_rxn_idxs = [3,4,7,8,11,12]

    add_loopless_constraints(model, S, internal_rxn_idxs)
    objective, _, solution, time, _ = optimize_model(model, time_limit=1800)
    nodes = MOI.get(model, MOI.NodeCount())

    model = build_model(S, lb, ub)
    objective_value_primal, time_primal, nodes_primal = loopless_fba_set_primal("simple_model", model, S, internal_rxn_idxs, nullspace_formulation=true, flux=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,])

    model = build_model(S, lb, ub)
    objective_value_primal_mu, time_primal_mu, nodes_primal_mu = loopless_fba_set_primal("simple_model", model, S, internal_rxn_idxs, nullspace_formulation=false, flux=[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,])

    @test isapprox(objective, objective_value_primal)
    @test isapprox(objective_value_primal, objective_value_primal_mu)
    @test nodes_primal == nodes_primal_mu == 0
    @test nodes >= nodes_primal
end

@testset "set optimal solution as primal in loopless FBA without nullspace formulation" begin
    organism = "iAF692"

    objective_value, solution, time, nodes = loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true, csv=false)
    loopless_fba_set_primal(organism, nullspace_formulation=true, load=false, time_limit=1800, csv=false)

    # objective_value_primal, time_primal, nodes_primal = loopless_fba_set_primal(organism, nullspace_formulation=true, flux=solution[1:690], load=false, time_limit=1800)

    # @test isapprox(objective_value_primal,objective_value, atol=0.001)
    # @test nodes_primal < nodes
end

# TODO: no assignment for G found for given flux
# @testset "set optimal solution as primal in loopless FBA with nullspace formulation" begin
#     organism = "iAF692"

#     objective_value, solution, time, nodes = loopless_fba_data(organism, time_limit=1800, nullspace_formulation=true)
#     objective_value_primal, time_primal, nodes_primal = loopless_fba_set_primal(organism, mu=false, flux=solution[1:690], load=false, time_limit=1800)

#     @test isapprox(objective_value_primal,objective_value, atol=0.001)
#     @test time_primal < time
#     @test nodes_primal < nodes
# end