using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra
using Boscia, FrankWolfe
include("functions.jl")

# could not extract json.gz
# molecular_model = load_model(ObjectModel, "data/iML1515.json")
# serialize("data/iML1515.js", molecular_model)

# build model
optimizer = SCIP.Optimizer

molecular_model = deserialize("data/iML1515.js")
print_model(molecular_model, "Escherichia coli str. K-12 substr. MG1655")

model = make_optimization_model(molecular_model, optimizer)
@show model
set_attribute(model, MOI.Silent(), true)

# FBA
optimize_model(model, print_objective=true)

# loopless FBA
add_loopless_constraints(molecular_model, model)
@show model
optimize_model(model, "loopless FBA")