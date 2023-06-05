using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Graphs

include("optimization_model.jl")
include("split_hyperarcs.jl")

# function constraints_to_matrix(model)
#     # @show model
#     @show list_of_constraint_types(model)
#     @show MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}
#     # @show MOI.get(model, MOI.ConstraintSet(), MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}(1))
#     # @show MOI.get(model, MOI.ConstraintSet(), all_constraints(model; include_variable_in_set_constraints = false)[1].index)

#     # @show model[:ubs][1]
# end 

"""
returns cycles as list of nodes found in graph G, where G is constructed using 
the reactions in transformed S that are non zero in the solution,
maps edges to reactions in transformed S
"""
function ubounded_cycles(S_transform, solution; ceiling=10^5, smallest_cycles=false)
    # filter non used reactions
    m, n = size(S_transform)
    # @show length(solution)
    @assert size(S_transform)[2] == length(solution)
    # ignore too small fluxes
    solution = [isapprox(0,i) ? 0 : i for i in solution]
    non_zero_reactions = findall(!iszero,solution)
    
    # map edges to reaction in transformed S_transform
    # build graph of used edges in solution
    edge_mapping = Dict()
    # build graph
    G = DiGraph()
    # @show size(S_transform)[1]
    add_vertices!(G, size(S_transform)[1])
    # @show non_zero_reactions
    for (idx,col) in enumerate(eachcol(S_transform))
        if idx in non_zero_reactions
            edge_mapping[idx] = [] # reaction is key
            metabolite_indices = findall(!iszero, col)
            @assert maximum(metabolite_indices) <= m
            # for internal reactions only
            if length(metabolite_indices) > 1
                @assert length(metabolite_indices) == 2
                if col[metabolite_indices[1]] * solution[idx] < 0
                    add_edge!(G, metabolite_indices[1], metabolite_indices[2])
                    push!(edge_mapping[idx], metabolite_indices[1]) 
                    push!(edge_mapping[idx], metabolite_indices[2])
                else
                    add_edge!(G, metabolite_indices[2], metabolite_indices[1])
                    push!(edge_mapping[idx], metabolite_indices[2]) 
                    push!(edge_mapping[idx], metabolite_indices[1])        
                end
            end
            # @show edge_mapping
        end    
    end

    # @show nv(G)
    # @show neighbors(G,1)
    # @show neighbors(G,2)
    # @show neighbors(G,3)
    # @show neighbors(G,4)
    # @show neighbors(G,5)
    
    # compute cycles
    # cycles are nodes in the network
    # first element ≠ last element
    if smallest_cycles 
        cycles = simplecycles_iter(G) 
        cycles_length = [length(c) for c in cycles]
        # @show issorted(cycles_length, rev=false)
        cycles_length, cycles = getindex.((cycles_length, cycles), (sortperm(cycles_length),))    
        # @show issorted(cycles_length, rev=false)
        cycles = [c for c in cycles if length(c) > 2]
        cycles = cycles[1:ceiling]
    else 
        cycles = simplecycles_iter(G, ceiling)
    end
    return cycles, edge_mapping, G
end 

"""
maps cycles found in G back to transformed S and S,
returns direction of reactions in cycle
"""
function unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
    edge_mapping_reverse = Dict(value => key for (key, value) in edge_mapping)
    # @show collect(keys(edge_mapping)) # reaction in transformed S to pair of nodes in G
    # @show collect(keys(edge_mapping_reverse)) # pair of nodes to reaction

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
    flux_directions = []
    for cycle in unbounded_cycles
        unbounded_cycle_original = []
        flux_directions_cycle = []
        for r in cycle 
            push!(unbounded_cycle_original, reaction_mapping_reverse[r])
            push!(flux_directions_cycle,solution[r])
        end
        push!(unbounded_cycles_original, unbounded_cycle_original)
        push!(flux_directions, [sign(i) for i in flux_directions_cycle]) #TODO: check when fluxes are negative
    end
    duplicates = [(i, count(==(i), unbounded_cycles_original)) for i in unique(unbounded_cycles_original)]
    
    # unique!(unbounded_cycles_original)

    return unbounded_cycles, unbounded_cycles_original, flux_directions
end
