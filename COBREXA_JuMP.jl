
using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP

optimizer = SCIP.Optimizer

molecular_model = deserialize("data/ec_e_coli_core.js")
println("E COLI COBRA MODEL")
println("number of metabolites : ", length(molecular_model.metabolites))
println("number of reactions : ", length(molecular_model.reactions))
println("number of genes : ", length(molecular_model.genes))
# @show model.annotations
# @show model.notes
println("objective : ", molecular_model.objective)
println("")

model = make_optimization_model(molecular_model, optimizer)
@show model
@show objective_function(model)

set_objective_coefficient(model, model.x[1], 0.25)