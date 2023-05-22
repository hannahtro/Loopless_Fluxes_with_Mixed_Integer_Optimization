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

# use split_hyperarcs.jl

# perform max flow as MIP