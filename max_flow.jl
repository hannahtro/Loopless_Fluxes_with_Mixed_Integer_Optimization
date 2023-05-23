using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Graphs

include("functions.jl")
include("split_hyperarcs.jl")

organism = "iAF692"

# transform S
molecular_model = deserialize("data/" * organism * ".js")
print_model(molecular_model)

# split hyperarcs
S = stoichiometry(molecular_model)
lb, ub = bounds(molecular_model)
S_transform, lb_transform, ub_transform = split_hyperarcs(S, lb, ub)
@show size(S_transform)
@show size(lb_transform), size(ub_transform)
m, n = size(S_transform)

# make optimization model
optimizer = SCIP.Optimizer
optimization_model = Model(optimizer)
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

# delete non used reactions
non_zero_reactions = findall(!iszero,solution)
S_transform = S_transform[1:end, 1:end .∉ [non_zero_reactions]] 
@show size(S_transform)

# build graph
G = DiGraph()
add_vertices!(G, size(S_transform)[1])
for (idx,col) in enumerate(eachcol(S_transform))
    # @show idx
    metabolite_indices = findall(!iszero, col)
    # @show metabolite_indices
    # @show col[metabolite_indices[1]]
    if col[metabolite_indices[1]] > 0 #TODO verify direction
        add_edge!(G, metabolite_indices[1], metabolite_indices[2])
    else
        add_edge!(G, metabolite_indices[2], metabolite_indices[1])
    end
end

@show edges(G)
@show weights(G)
@show vertices(G)

# compute cycles
cycles = [simplecycles(G)]
@show length(cycles[1])

# get original reactions

