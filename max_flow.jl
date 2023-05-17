using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra

include("functions.jl")

organism = "iAF692"

# transform S
# build model
optimizer = SCIP.Optimizer

molecular_model = deserialize("data/" * organism * ".js")
print_model(molecular_model)

model = make_optimization_model(molecular_model, optimizer)
@show model

constraints = all_constraints(model, include_variable_in_set_constraints = false)
@show constraints[1:10]

# perform max flow as MIP