using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe
include("functions.jl")

# build model
optimizer = SCIP.Optimizer

molecular_model = deserialize("data/ec_e_coli_core.js")

print_model(molecular_model, "E COLI COBRA MODEL")

# add loop
r_1 = ReactionBidirectional("r_1", Dict("glc__D_e" => -1.0, "gln__L_c" => 1.0))
r_2 = ReactionBidirectional("r_2", Dict("gln__L_c" => -1.0, "gln__L_e" => 1.0))
r_3 = ReactionBidirectional("r_3", Dict("gln__L_e" => -1.0, "glc__D_e" => 1.0))

add_reactions!(molecular_model, [r_1, r_2, r_3])

model = make_optimization_model(molecular_model, optimizer)
@show model
set_attribute(model, MOI.Silent(), true)

# FBA
optimize_model(model, print_objective=true)

# FBA with modified objective
x = model[:x]
@objective(model, Max, sum(x))
optimize_model(model)

# loopless FBA
add_loopless_constraints(molecular_model, model)
optimize_model(model, "loopless FBA")
