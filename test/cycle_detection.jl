using Test 
using GraphPlot, Cairo, Compose

include("../src/cycle_detection.jl")
include("../src/loopless_constraints.jl")

# S = [[0,1,1,-1,0] [-1,1,1,0,0] [0,0,-1,0,1] [0,0,0,1,-1] [1,0,0,0,0] [0,0,0,-1,0] [0,-1,0,0,0]]
# # @show S
# # @show size(S)
# m, num_reactions = size(S)
# lb = [-10,-10,-10,-10,0,0,0]
# ub = [10,10,10,10,10,10,10]

# @testset "block cycle in transformed S of simple model" begin
#     S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
#     # @show size(S_transform)
#     # @show S_transform'
#     _, num_reactions_transform = size(S_transform)
#     @test size(S_transform)[1] == size(S)[1]

#     model = build_fba_model(S_transform, lb_transform, ub_transform)
#     x = model[:x]
#     @objective(model, Max, x[2]+x[5]+x[6])
#     _, _, solution_loop, _, _ = optimize_model(model)
#     @show solution_loop
#     @test length(solution_loop) == size(S_transform)[2]

#     # get original reactions
#     cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_loop)
#     nodelabel = ["A", "B", "C", "D", "E"]
#     # gplothtml(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout)

#     # @show cycles
#     @test length(cycles) == 1
#     # @show edge_mapping

#     unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_loop, reaction_mapping)
#     @show unbounded_cycles
#     @show unbounded_cycles_original
#     # @show flux_directions

#     # block cycle in transformed S
#     optimization_model = build_fba_model(S_transform, lb_transform, ub_transform)
#     internal_rxn_idxs = [1,2,3,4,5,6]
#     add_loopless_constraints(optimization_model, S_transform, internal_rxn_idxs)

#     # @show all_constraints(optimization_model, include_variable_in_set_constraints=true)
#     block_cycle_constraint(optimization_model, unbounded_cycles, flux_directions, internal_rxn_idxs, S_transform)
#     # @show all_constraints(optimization_model, include_variable_in_set_constraints=true)

#     x = optimization_model[:x]
#     @objective(optimization_model, Max, x[2]+x[5]+x[6])
#     _, _, solution, _, _ = optimize_model(optimization_model)

#     @show solution[1:num_reactions_transform] # x, a, G
#     # @show solution[num_reactions_transform+1:num_reactions_transform+length(internal_rxn_idxs)]
#     # @show solution[num_reactions_transform+1+length(internal_rxn_idxs):end]
#     @test solution[1:num_reactions_transform] != solution_loop
# end

# TODO: verify thermo feasibilty of cycle
# thermodynamically feasible but not biological valid as no external source is take up
# @testset "block cycle in S of simple model" begin
#     model = build_fba_model(S, lb, ub)
#     x = model[:x]
#     @objective(model, Max, x[1]+x[3]+x[4])
#     _, _, solution_loop, _, _ = optimize_model(model)
#     @show solution_loop

#     # check if cycle exists
#     # get original reactions
#     S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution_loop)
#     # @show solution_transform
#     cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
#     @test length(cycles) == 1
#     @show cycles
#     # @show edge_mapping

#     # nodelabel = ["A", "B", "C", "D", "E"]
#     # gplothtml(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout)

#     unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)
#     @show unbounded_cycles
#     @show unbounded_cycles_original
#     # @show flux_directions

#     # block cycle in original S
#     optimization_model = build_fba_model(S, lb, ub)
#     x = optimization_model[:x]
#     @objective(optimization_model, Max, x[1]+x[3]+x[4])
#     internal_rxn_idxs = [1,2,3,4]
#     add_loopless_constraints(optimization_model, S, internal_rxn_idxs)

#     # @show optimization_model
#     num_blocked_cycles = block_cycle_constraint(optimization_model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S)
#     @show num_blocked_cycles

#     # @show optimization_model
#     _, _, solution, _, _ = optimize_model(optimization_model)

#     @show solution[1:num_reactions] # x
#     @test solution[1:num_reactions] != solution_loop
#     # @show solution[num_reactions+1:num_reactions+length(internal_rxn_idxs)] # a
#     # @show solution[num_reactions+1+length(internal_rxn_idxs):end] # G

#     S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution)
#     cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
#     #TODO: add node position
#     # locs_x = [2,1,3,2,4]
#     # locs_y = [1,2,2,3,3]
#     # locs_x, locs_y = shell_layout(G)
#     # @show locs_x, locs_y
#     # gplothtml(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout)
#     # gplothtml(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout, locs_x_in=locs_x, locs_y_in=locs_y)
# end

# @testset "block cycle in simple graph with feasible cycle" begin
#     S = [[1,0,0,0,0,0,0,0,0] [0,0,0,0,0,1,0,0,0] [0,1,0,0,0,-1,0,0,0] [0,-1,0,0,0,0,1,0,0] [0,0,0,0,0,0,-1,0,0] [0,0,0,0,-1,0,0,0,0] [-1,-1,1,1,0,0,0,0,0] [0,1,0,-1,1,0,0,0,0] [0,0,-1,0,0,0,0,0,0] [0,0,0,0,0,0,0,1,0] [0,0,0,1,0,0,0,-1,0] [0,0,0,-1,0,0,0,0,1] [0,0,0,0,0,0,0,0,-1]]
#     m, num_reactions = size(S)
#     @show m, num_reactions
#     lb = [0,0,0,0,0,0,0,0,0,0,0,0,0]
#     ub = [10,10,10,10,10,10,10,10,10,10,10,10,10]
#     # @show length(lb), length(ub)
#     # println("-----------------------------------")

#     # feasible loop
#     model = build_fba_model(S, lb, ub)
#     x = model[:x]
#     @objective(model, Max, x[1]+x[3]+x[4])
#     add_loopless_constraints(model, S, [3,4,7,8,11,12])
#     _, _, solution, _, _ = optimize_model(model)
#     @show solution
#     S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution)
#     cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
#     @test length(cycles) == 1
#     println("# thermo feasible cylces: ", length(cycles))
#     println("-----------------------------------")

#     # verify that found loop is thermodynamically feasible
#     model = build_fba_model(S, lb, ub)
#     x = model[:x]
#     @objective(model, Max, x[1]+x[3]+x[4])
#     # add_loopless_constraints(model, S, [3,4,7,8,11,12])
#     _, _, solution_loop, _, _ = optimize_model(model)
#     # @show solution_loop
#     # check if cycle exists
#     # get original reactions
#     S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution_loop)
#     # @show solution_transform
#     cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
#     @test length(cycles) == 1
#     # @show edge_mapping
#     # nodelabel = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
#     # gplothtml(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout)
#     unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)
#     # @show unbounded_cycles
#     @show unbounded_cycles_original
#     internal_rxn_idxs = [3,4,7,8,11,12]
#     feasible = thermo_feasible(unbounded_cycles_original[1], flux_directions[1], S)
#     @test feasible
# end

@testset "block cycle in simple graph with infeasible cycle" begin
    S = [[1,0,0] [-1,1,0] [0,-1,1] [-1,0,1] [0,0,-1]]
    lb = [0,-10,-10,-10,0]
    ub = [10,10,10,10,10]
    m, num_reactions = size(S)
    @show m, num_reactions
    println("-----------------------------------")

    # infeasible loop
    model = build_fba_model(S, lb, ub)
    x = model[:x]
    @objective(model, Max, x[2]+x[3]+x[4])
    add_loopless_constraints(model, S, [2,3,4])
    print(model)
    _, _, solution, _, _ = optimize_model(model)
    @show solution[1:num_reactions]
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution)
    cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
    @test length(cycles) == 0
    println("thermo feasible cylce: ", cycles)
    println("-----------------------------------")

    # infeasible loop using add_loopless_constraints_mu
    model = build_fba_model(S, lb, ub)
    x = model[:x]
    @objective(model, Max, x[2]+x[3]+x[4])
    add_loopless_constraints_mu(model, S, [2,3,4])
    print(model)
    _, _, solution, _, _ = optimize_model(model)
    @show solution[1:num_reactions]
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution)
    cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
    @test length(cycles) == 0
    println("thermo feasible cylce: ", cycles)
    println("-----------------------------------")

    # infeasible loop using add_loopless_constraints_mu_reduced
    model = build_fba_model(S, lb, ub)
    x = model[:x]
    @objective(model, Max, x[2]+x[3]+x[4])
    add_loopless_constraints_mu_reduced(model, S, [2,3,4])
    print(model)
    _, _, solution, _, _ = optimize_model(model)
    @show solution[1:num_reactions]
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution)
    cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
    @test length(cycles) == 0
    println("thermo feasible cylce: ", cycles)
    println("-----------------------------------")

    # verify that found loop is thermodynamically feasible
    model = build_fba_model(S, lb, ub)
    x = model[:x]
    @objective(model, Max, x[2]-x[4])
    _, _, solution_loop, _, _ = optimize_model(model)
    @show solution_loop
    # check if cycle exists
    # get original reactions
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform = split_hyperarcs(S, lb, ub, solution_loop)
    @test S == S_transform
    @test solution_loop == solution_transform
    cycles, edge_mapping, G = ubounded_cycles(S_transform, solution_transform)
    @test length(cycles) == 1
    nodelabel = ["A", "B", "C"]
    # draw(PDF("cycle.pdf", 16cm, 16cm), gplot(G, nodesize=0.1, nodelabel=nodelabel, layout=circular_layout))

    # @show edge_mapping
    unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)
    @show unbounded_cycles_original
    internal_rxn_idxs = [2,3,4]
    feasible = thermo_feasible(unbounded_cycles_original[1], flux_directions[1], S)
    @test !feasible
end

# @testset "block cycle in iAF692" begin
#     # test organism
#     organism = "iAF692"

#     # transform S
#     molecular_model = deserialize("../data/" * organism * ".js")
#     print_model(molecular_model, organism)

#     # split hyperarcs
#     S = stoichiometry(molecular_model)
#     # @show size(S)
#     lb, ub = bounds(molecular_model)
#     S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
#     @test size(S_transform)[2] == length(lb_transform) == length(ub_transform)
#     m, n = size(S_transform)

#     optimization_model = build_fba_model(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
#     _, _, solution, _, _ = optimize_model(optimization_model)
#     # @show size(solution)

#     # get original reactions
#     cycles, edge_mapping, _ = ubounded_cycles(S_transform, solution, ceiling=100)
#     @test length(cycles) == 100
#     # @show cycles
#     # @show edge_mapping
#     unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
#     # @show unbounded_cycles
#     # @show flux_directions
#     add_loopless_constraints(molecular_model, optimization_model)

#     # @show optimization_model
#     # @show unbounded_cycles_original
#     internal_rxn_idxs = [
#         ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
#         !is_boundary(reaction_stoichiometry(molecular_model, rid))
#     ]
#     num_blocked_cycles = block_cycle_constraint(optimization_model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S)
# end


# TODO: test feasibility mu