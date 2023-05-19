using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra

include("functions.jl")


organism = "iJR904"

# build model
optimizer = SCIP.Optimizer

molecular_model = deserialize("data/" * organism * ".js")
print_model(molecular_model, organism)

model = make_optimization_model(molecular_model, optimizer)
# @show model

list_of_constraints = all_constraints(model; include_variable_in_set_constraints = true)[1:10]
for reaction in list_of_constraints
    println(reaction)
    println(constraint_object(reaction))
end