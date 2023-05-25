using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Graphs

include("functions.jl")
include("split_hyperarcs.jl")

# function constraints_to_matrix(model)
#     # @show model
#     @show list_of_constraint_types(model)
#     @show MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}
#     # @show MOI.get(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}(1))
#     # @show MOI.get(model, MOI.ConstraintSet(), all_constraints(model; include_variable_in_set_constraints = false)[1].index)

#     # @show model[:ubs][1]
# end 

#TODO: why do we get loops in max_flow?
function build_model(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer, cycles=[], flux_values=[])
    # make optimization model
    optimization_model = Model(optimizer)
    _, n = size(S_transform)
    @show size(S_transform)
    @show size(lb_transform)
    @show size(ub_transform)

    @variable(optimization_model, x[1:n])
    @constraint(optimization_model, mb, S_transform * x .== 0) # mass balance #TODO set coefficients to -1/1?
    @constraint(optimization_model, lbs, lb_transform .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= ub_transform) # upper bounds
    # @show optimization_model

    #TODO: remove this part
    if !isempty(cycles)
        for (idx,cycle) in enumerate(cycles)
            cycle_vars = [x[i] for i in cycle]
            @show cycle_vars
            @constraint(optimization_model, sum(cycle_vars) <= flux_values[idx])
        end
    end
    @objective(optimization_model, Max, sum(x)) #TODO use original objective

    return optimization_model
end

function ubounded_cycles(S_transform, solution)
    # filter non used reactions
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

    # map edges to reaction in transformed S_transform
    # build graph of used edges in solution
    edge_mapping = Dict()
    # build graph
    G = DiGraph()
    add_vertices!(G, size(S_transform)[1])
    @show non_zero_reactions
    for (idx,col) in enumerate(eachcol(S_transform))
        if idx in non_zero_reactions
            edge_mapping[idx] = []
            metabolite_indices = findall(!iszero, col)
            # for internal reactions only
            if length(metabolite_indices) > 1
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
    end

    @show length(edges(G))
    # @show weights(G)
    # @show vertices(G)

    # compute cycles
    cycles = simplecycles_iter(G, 10^5)
    @show length(cycles)
    return cycles, edge_mapping
end 

function unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
    edge_mapping_reverse = Dict(value => key for (key, value) in edge_mapping)
    # @show edge_mapping_reverse

    unbounded_cycles = []
    for cycle in cycles
        unbounded_cycle = []
        temp_cycle = deepcopy(cycle)
        push!(temp_cycle, cycle[1])
        # @show temp_cycle
        # @show cycle
        for (idx,r) in enumerate(cycle)
            key = [r,temp_cycle[idx+1]]
            push!(unbounded_cycle,(edge_mapping_reverse[key]))
        end
        push!(unbounded_cycles, unbounded_cycle)
    end

    # map reactions participating in cycles to reactions in S
    # @show reaction_mapping
    reaction_mapping_reverse = Dict()
    for (key,value) in reaction_mapping
        for val in value
            reaction_mapping_reverse[val]=key
        end
    end
    # @show reaction_mapping_reverse

    # for each cycle, compute the current flux through reactions
    unbounded_cycles_original = [] 
    flux_values = []
    for cycle in unbounded_cycles
        unbounded_cycle_original = []
        flux_value = []
        for r in cycle 
            push!(unbounded_cycle_original, reaction_mapping_reverse[r])
            push!(flux_value,solution[r])
        end
        push!(unbounded_cycles_original, unbounded_cycle_original)
        push!(flux_values, sum(flux_value) - minimum(flux_value)) #TODO: check when fluxes are negative
    end
    unique!(unbounded_cycles_original)

    return unbounded_cycles, unbounded_cycles_original, flux_values
end

# test network
println("FBA of transformed S with loop")
S = [[0,1,1,-1,0] [-1,1,1,0,0] [0,0,-1,0,1] [0,0,0,1,-1] [1,0,0,0,0] [0,0,0,-1,0] [0,-1,0,0,0]]
# @show S
@show size(S)
_, num_reactions = size(S)
lb = [-10,-10,-10,-10,0,0,0]
ub = [10,10,10,10,10,10,10]
S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
@show size(S_transform)
@show S_transform'
_, num_reactions_transform = size(S_transform)
@assert size(S_transform)[1] == size(S)[1]

model = build_model(S_transform, lb_transform, ub_transform)
x = model[:x]
@objective(model, Max, x[2]+x[5]+x[6])
_, solution, _, _ = optimize_model(model)
@show solution
@assert length(solution) == size(S_transform)[2]

# get original reactions
cycles, edge_mapping = ubounded_cycles(S_transform, solution)
@show cycles
@show edge_mapping

unbounded_cycles, unbounded_cycles_original, flux_values = unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
@show unbounded_cycles
@show flux_values

println("")
println("loopless FBA of transformed S with blocked cycle")
# block cycle in transformed S
optimization_model = build_model(S_transform, lb_transform, ub_transform)
internal_rxn_idxs = [1,2,3,4,5,6]
add_loopless_constraints(optimization_model, S_transform, internal_rxn_idxs)

# @show optimization_model
a = optimization_model[:a]
x = optimization_model[:x]
for cycle in unbounded_cycles_original
    cycle_vars = [a[i] for i in cycle]
    @show cycle_vars
    @constraint(optimization_model, sum(cycle_vars) >= 1)
end
# @show optimization_model
x = optimization_model[:x]
@objective(optimization_model, Max, x[2]+x[5]+x[6])
_, solution, _, _ = optimize_model(optimization_model)

@show solution[1:num_reactions_transform] # x, a, G
@show solution[num_reactions_transform+1:num_reactions_transform+length(internal_rxn_idxs)]
@show solution[num_reactions_transform+1+length(internal_rxn_idxs):end]


println("")
println("FBA of S with loop")
model = build_model(S, lb, ub)
x = model[:x]
@objective(model, Max, x[1]+x[3]+x[4])
_, solution, _, _ = optimize_model(model)
@show solution

println("")
println("loopless FBA of S with blocked cycle")
# # block cycle in original S
optimization_model = build_model(S, lb, ub)
x = optimization_model[:x]
@objective(optimization_model, Max, x[1]+x[3]+x[4])
internal_rxn_idxs = [1,2,3,4]
add_loopless_constraints(optimization_model, S, internal_rxn_idxs)

# @show optimization_model
a = optimization_model[:a] #TODO: ensure direction of flux
x = optimization_model[:x]
for cycle in unbounded_cycles_original
    cycle_vars = [a[i] for i in cycle]
    @show cycle_vars
    @constraint(optimization_model, sum(cycle_vars) >= 1)
end
# @show optimization_model
_, solution, _, _ = optimize_model(optimization_model)

@show solution[1:num_reactions] # x
@show solution[num_reactions+1:num_reactions+length(internal_rxn_idxs)] # a
@show solution[num_reactions+1+length(internal_rxn_idxs):end] # G

# # test organism
# organism = "iAF692"

# # transform S
# molecular_model = deserialize("data/" * organism * ".js")
# print_model(molecular_model)

# # split hyperarcs
# S = stoichiometry(molecular_model)
# lb, ub = bounds(molecular_model)
# S_transform, lb_transform, ub_transform, reaction_mapping = split_hyperarcs(S, lb, ub)
# # @show size(S_transform)
# # @show size(lb_transform), size(ub_transform)
# m, n = size(S_transform)

# solution = max_flow(S_transform, lb_transform, ub_transform; optimizer=SCIP.Optimizer)
# # @show solution
# # get original reactions
# cycles, edge_mapping = ubounded_cycles(S_transform, solution)
# # @show cycles[1]
# # @show edge_mapping
# unbounded_cycles, flux_values = unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
# # @show unbounded_cycles
# # @show flux_values
