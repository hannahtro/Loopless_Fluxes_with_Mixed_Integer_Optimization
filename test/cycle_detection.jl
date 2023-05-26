using Test 

include("../src/cycle_detection.jl")
include("../src/loopless_constraints.jl")

S = [[0,1,1,-1,0] [-1,1,1,0,0] [0,0,-1,0,1] [0,0,0,1,-1] [1,0,0,0,0] [0,0,0,-1,0] [0,-1,0,0,0]]
# @show S
# @show size(S)
_, num_reactions = size(S)
lb = [-10,-10,-10,-10,0,0,0]
ub = [10,10,10,10,10,10,10]

@testset "block cycles in transformed S of simple model" begin
    S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
    # @show size(S_transform)
    # @show S_transform'
    _, num_reactions_transform = size(S_transform)
    @test size(S_transform)[1] == size(S)[1]

    model = build_model(S_transform, lb_transform, ub_transform)
    x = model[:x]
    @objective(model, Max, x[2]+x[5]+x[6])
    _, solution_loop, _, _ = optimize_model(model)
    @show solution_loop
    @test length(solution_loop) == size(S_transform)[2]

    # get original reactions
    cycles, edge_mapping = ubounded_cycles(S_transform, solution_loop)
    # @show cycles
    @test length(cycles) == 1
    # @show edge_mapping

    unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_loop, reaction_mapping)
    @show unbounded_cycles
    @show unbounded_cycles_original
    # @show flux_directions

    # block cycle in transformed S
    optimization_model = build_model(S_transform, lb_transform, ub_transform)
    internal_rxn_idxs = [1,2,3,4,5,6]
    add_loopless_constraints(optimization_model, S_transform, internal_rxn_idxs)

    # @show optimization_model
    block_cycle_constraint(optimization_model, unbounded_cycles_original, flux_directions)
    # @show optimization_model

    x = optimization_model[:x]
    @objective(optimization_model, Max, x[2]+x[5]+x[6])
    _, solution, _, _ = optimize_model(optimization_model)

    @show solution[1:num_reactions_transform] # x, a, G
    # @show solution[num_reactions_transform+1:num_reactions_transform+length(internal_rxn_idxs)]
    # @show solution[num_reactions_transform+1+length(internal_rxn_idxs):end]
    @test solution[1:num_reactions_transform] != solution_loop
end

@testset "block cycles in S of simple model" begin
    model = build_model(S, lb, ub)
    x = model[:x]
    @objective(model, Max, x[1]+x[3]+x[4])
    _, solution_loop, _, _ = optimize_model(model)
    @show solution_loop

    # check if cycle exists
    # get original reactions
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution_loop)
    # @show solution_transform
    cycles, edge_mapping = ubounded_cycles(S_transform, solution_transform)
    @test length(cycles) == 1
    # @show edge_mapping

    unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)
    @show unbounded_cycles
    @show unbounded_cycles_original
    # @show flux_directions

    # block cycle in original S
    optimization_model = build_model(S, lb, ub)
    x = optimization_model[:x]
    @objective(optimization_model, Max, x[1]+x[3]+x[4])
    internal_rxn_idxs = [1,2,3,4]
    add_loopless_constraints(optimization_model, S, internal_rxn_idxs)

    # @show optimization_model
    block_cycle_constraint(optimization_model, unbounded_cycles_original, flux_directions)

    # @show optimization_model
    _, solution, _, _ = optimize_model(optimization_model)

    @test solution[1:num_reactions] != solution_loop
    @show solution[1:num_reactions] # x
    # @show solution[num_reactions+1:num_reactions+length(internal_rxn_idxs)] # a
    # @show solution[num_reactions+1+length(internal_rxn_idxs):end] # G
end