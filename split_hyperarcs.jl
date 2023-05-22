using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra

include("functions.jl")


organism = "iJR904"

# build model
# optimizer = SCIP.Optimizer

# molecular_model = deserialize("data/" * organism * ".js")
# print_model(molecular_model, organism)

# S = stoichiometry(molecular_model)
# S_transpose = S'

# @show size(S)
# @show size(S_transpose)

# S_transform = [;;]

# for row in eachrow(S_transpose)
#     @show row[1:10]
# end    

S = [[0,1,1,-1] [-1,1,1,0]]
S_transpose = S'
n,m = size(S_transpose)
@show n,m # n reactions, m metabolites
S_transform = [;;]

for row in eachrow(S_transpose)
    @show row
    
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
    @show backward_coef, backward_arc_m
    @show forward_coef, forward_arc_m
    
    # add splitted reactions
    for (forward_idx,forward_m) in enumerate(forward_arc_m)
        for (backard_idx,backward_m) in enumerate(backward_arc_m)
            split_arc = zeros(m)
            split_arc[forward_m] = forward_coef[forward_idx]
            split_arc[backward_m] = backward_coef[backard_idx]
            @show split_arc
            vcat(S_transform, split_arc)
            @show S_transform
        end
    end  
end

@show S_transform

# model = make_optimization_model(molecular_model, optimizer)
# # @show model

# list_of_constraints = all_constraints(model; include_variable_in_set_constraints = true)[1:10]

# for reaction in list_of_constraints
#     constraint = constraint_object(reaction)
#     vars = constraint.func.terms.keys
#     vals = constraint.func.terms.vals

#     pos_vars = []
#     neg_vars = []
#     pos_vals = []
#     neg_vals = []

#     # separate forward and backward arcs
#     for (id,v) in enumerate(vals)
#         if v > 0
#             push!(pos_vals, v)
#             push!(pos_vars, vars[id])
#         else
#             push!(neg_vals, v)

#             push!(neg_vars, vars[id])  
#         end
#     end

#     # add splitted reactions
#     for (pos_id,pos_reaction) in enumerate(pos_vars)
#         for (neg_id,neg_reaction) in enumerate(neg_vars)
#             @constraint(transformed_model, )
#         end
#     end   

#     # fix arc direction
# end
