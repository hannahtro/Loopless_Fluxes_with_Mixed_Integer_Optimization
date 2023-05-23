using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Graphs

include("functions.jl")
include("split_hyperarcs.jl")

function max_flow(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
    # make optimization model
    optimization_model = Model(optimizer)
    _, n = size(S_transform)
    @variable(optimization_model, x[1:n])
    @constraint(optimization_model, mb, S_transform * x .== 0) # mass balance #TODO set coefficients to -1/1?
    @constraint(optimization_model, lbs, lb_transform .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= ub_transform) # upper bounds
    # @show optimization_model

    # perform max flow as MIP
    @objective(optimization_model, Max, sum(x)) #TODO use original objective
    set_attribute(optimization_model, MOI.Silent(), true)
    optimize!(optimization_model)
    solution = [value(var) for var in all_variables(optimization_model)]
    # @show solution
end

function ubounded_cycles(S_transform, solution)
    # delete non used reactions
    println("")
    # @show size(S_transform)
    # @show S_transform
    non_zero_reactions = findall(!iszero,solution)
    # @show non_zero_reactions
    S_transform_reduced = []
    for (id,row) in enumerate(eachcol(S_transform))
        if id in non_zero_reactions
            push!(S_transform_reduced, row)
        end
    end
    # list to array
    S_transform_reduced = mapreduce(permutedims, vcat, S_transform_reduced)'
    @show size(S_transform_reduced)

    edge_mapping = Dict()
    # build graph
    G = DiGraph()
    add_vertices!(G, size(S_transform)[1])
    for (idx,col) in enumerate(eachcol(S_transform))
        if idx in non_zero_reactions
            edge_mapping[idx] = []
            # @show idx
            metabolite_indices = findall(!iszero, col)
            # @show metabolite_indices
            # @show col[metabolite_indices[1]]
            if col[metabolite_indices[1]] > 0 #TODO verify direction
                add_edge!(G, metabolite_indices[1], metabolite_indices[2])
                push!(edge_mapping[idx], metabolite_indices[1]) 
                push!(edge_mapping[idx], metabolite_indices[2])
            else
                add_edge!(G, metabolite_indices[2], metabolite_indices[1])
                push!(edge_mapping[idx], metabolite_indices[2]) 
                push!(edge_mapping[idx], metabolite_indices[1])        
            end
        end    
    end

    # @show edges(G)
    # @show weights(G)
    # @show vertices(G)

    # compute cycles
    cycles = simplecycles(G)
    # @show length(cycles)
    return cycles, edge_mapping
end 

function unbounded_cycles_S(cycles, edge_mapping)
    edge_mapping_reverse = Dict(value => key for (key, value) in edge_mapping)
    # @show edge_mapping_reverse

    unbounded_cycles = []
    for cycle in cycles
        unbounded_cycle = []
        temp_cycle = deepcopy(cycle)
        append!(temp_cycle, cycle[1])
        # @show temp_cycle
        # @show cycle
        for (idx,m) in enumerate(cycle)
            key = [m,temp_cycle[idx+1]]
            push!(unbounded_cycle,(edge_mapping_reverse[key]))
        end
        push!(unbounded_cycles, unbounded_cycle)
    end
    return unbounded_cycles
end

# test network
S = [[0,1,1,-1,0] [-1,1,1,0,0] [0,0,-1,0,1] [0,0,0,1,-1]]
# @show S
# @show size(S)
lb = [-10,-10,-10,-10]
ub = [10,10,10,10,10]
S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
@show S_transform'
@assert size(S_transform)[1] == size(S)[1]

solution = max_flow(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
@show solution
@assert length(solution) == size(S_transform)[2]

# get original reactions
cycles, edge_mapping = ubounded_cycles(S_transform, solution)
@show cycles
@show edge_mapping

unbounded_cycles = unbounded_cycles_S(cycles, edge_mapping)
@show unbounded_cycles

# # test organism
# organism = "iAF692"

# # transform S
# molecular_model = deserialize("data/" * organism * ".js")
# print_model(molecular_model)

# # split hyperarcs
# S = stoichiometry(molecular_model)
# lb, ub = bounds(molecular_model)
# S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
# @show size(S_transform)
# @show size(lb_transform), size(ub_transform)
# m, n = size(S_transform)

# solution = max_flow(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
# # get original reactions
# cycles = ubound_cycles(S_transform, solution)
# @show cycles[1]
# # @show edge_mapping