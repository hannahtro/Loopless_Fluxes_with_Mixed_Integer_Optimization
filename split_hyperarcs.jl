using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra

include("functions.jl") 

function split_hyperarcs(S, lb, ub)
    S_transpose = S'
    n,m = size(S_transpose)
    @show n,m # n reactions, m metabolites
    S_transform = []
    lb_transform = []
    ub_transform = []

    for (idx,row) in enumerate(eachrow(S_transpose))
        # @show row # reaction
        
        forward_arc_m = [] # metabolite
        backward_arc_m = []
        forward_coef = []
        backward_coef = []

        # separate forward and backward arcs
        for (id,coef) in enumerate(row)
            if coef > 0
                push!(forward_coef, coef)
                push!(forward_arc_m, id)
            elseif coef < 0
                push!(backward_coef, coef)
                push!(backward_arc_m, id)  
            end
        end
        # @show backward_coef, backward_arc_m
        # @show forward_coef, forward_arc_m
        
        # add splitted reactions
        for (forward_idx,forward_m) in enumerate(forward_arc_m)
            for (backard_idx,backward_m) in enumerate(backward_arc_m)
                push!(lb_transform,lb[idx])
                push!(ub_transform,ub[idx])
                split_arc = zeros(m)
                split_arc[forward_m] = forward_coef[forward_idx]
                split_arc[backward_m] = backward_coef[backard_idx]
                push!(S_transform, split_arc)
                # @show split_arc
            end
        end  
    end

    S_transform = mapreduce(permutedims, vcat, S_transform)
    S_transform = S_transform'
    return S_transform, lb_transform, ub_transform
end

# test network
# S = [[0,1,1,-1] [-1,1,1,0]]
# @show S
# lb = [-5,-10]
# ub = [5,10]
# S_transform, lb_transform, ub_transform = split_hyperarcs(S, lb, ub)
# @show S_transform
# @show lb_transform, ub_transform

# # bigger model
# organism = "iJR904"

# # build model
# optimizer = SCIP.Optimizer

# molecular_model = deserialize("data/" * organism * ".js")
# print_model(molecular_model, organism)

# S = stoichiometry(molecular_model)
# lb,ub = bounds(molecular_model)
# S_transform, lb_transform, ub_transform = split_hyperarcs(S, lb, ub)
# @show size(S_transform)
# @show size(lb_transform), size(ub_transform)

# [`variables`](@ref)       # number of reactions
# [`metabolites`](@ref)     # same
# [`stoichiometry`](@ref)   # computed
# [`bounds`](@ref)          # upper and lower_bound on each reaction
# [`objective`](@ref)       # set correct reaction

# struct SplitModel <: AbstractMetabolicModel
#     stoichiometry
#     original_model::AbstractMetabolicModel
# end

# Everything.n_reactions(m::SplitModel) = size(m.stoichiometry)[2]
# Everything.n_metabolites(m::SplitModel) = size(m.stoichiometry)[1]

# Everything.objective(m::SplitModel) = m.original_model.objective # TODO update
# Everything.metabolites(m::SplitModel) = m.original_model.metabolites
# Everything.stoichiometry(m::SplitModel) = m.stoichiometry
# Everything.reactions(m::SplitModel) = ["rxn$i" for i in 1:n_reactions(m)]

# function Everything.bounds(m::SplitModel) # TODO update
#     m = size(m.stoichiometry)[2]
#     zeros(m) # lower bounds
#     ones(m) # upper bounds
# end

# transformed_model = SplitModel(S_transform, molecular_model)
# @show objective(transformed_model)
# # @show length(metabolites(transformed_model))
# # @show length(reactions(transformed_model))
# # @show length(Everything.bounds(transformed_model))

# make_optimization_model(transformed_model, optimizer)