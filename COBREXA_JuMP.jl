
using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe
include("functions.jl")

optimizer = SCIP.Optimizer

molecular_model = deserialize("data/ec_e_coli_core.js")
print_model(molecular_model, "E COLI COBRA MODEL")

model = make_optimization_model(molecular_model, optimizer)
set_attribute(model, MOI.Silent(), true)
# FBA
optimize_model(model, print_objective=true)

# get reference to variables
x = model[:x]
reference_flux =  MOI.get.(model, MOI.VariablePrimal(), x)
# @show reference_flux

# loopless FBA
add_loopless_constraints(molecular_model, model)
optimize_model(model, "loopless FBA")

# block reaction
set_lower_bound(x[1], 0)
set_upper_bound(x[1], 0)

optimize_model(model, "loopless FBA with blocked reaction")
knock_out_flux =  MOI.get.(model, MOI.VariablePrimal(), x)
# @show knock_out_flux

# MOMA
# modify objective
L = - reference_flux
Q = I(length(x))
@objective(model, Min, 1/2 * x' * Q * x + L' * x)
# @show objective_function(model)
@show typeof(objective_function(model))
optimize_model(model, "loopless MOMA")

# Boscia
println("")
println("loopless MOMA in Boscia")
println("----------------------------------")

function f(x)
    length_var = size(Q)[1]
    f_x = 1/2 * x[1:length_var]' * Q * x[1:length_var] + L' * x[1:length_var]
end
# @show f(ones(length(x)))

function grad!(storage, x)
    length_var = size(Q)[1]
    storage[1:length_var] = Q * x[1:length_var] + L
    storage[length_var+1:length(x)] .= 0
end

set_objective_sense(model, FEASIBILITY_SENSE)
moi_model = backend(model)
lmo = FrankWolfe.MathOptLMO(moi_model)
x, _, result = Boscia.solve(f, grad!, lmo, verbose=true)
