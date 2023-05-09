
using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra

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
println("FBA")
println("----------------------------------")
model = make_optimization_model(molecular_model, optimizer)
set_attribute(model, MOI.Silent(), true)
@show model
@show objective_function(model)

# get reference to variables
x = model[:x]

optimize!(model)
reference_flux =  MOI.get.(model, MOI.VariablePrimal(), x)
# @show reference_flux

# @show MOI.get(model, MOI.TerminationStatus())
@show MOI.get(model, MOI.ObjectiveValue())
# @show  MOI.get.(model, MOI.VariablePrimal(), x)

println("")
println("loopless FBA")
println("----------------------------------")
# loopless model
internal_rxn_idxs = [
    ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
    !is_boundary(reaction_stoichiometry(molecular_model, rid))
]

@show length(internal_rxn_idxs)

N_int = nullspace(Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])) # no sparse nullspace function

a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

for (cidx, ridx) in enumerate(internal_rxn_idxs)
    @constraint(model, -1000 * (1 - a[cidx]) <= x[ridx])
    @constraint(model, x[ridx] <= 1000 * a[cidx])

    @constraint(
        model,
        -1000 * a[cidx] + (1 - a[cidx]) <= G[cidx]
    )
    @constraint(
        model,
        G[cidx] <=
        -a[cidx] + 1000 * (1 - a[cidx])
    )
end

@constraint(model, N_int' * G .== 0)

@show model
optimize!(model)
@show MOI.get(model, MOI.ObjectiveValue())

println("")
println("loopless FBA with blocked reaction")
println("----------------------------------")
# block reaction
set_lower_bound(x[1], 0)
set_upper_bound(x[1], 0)

@show model

optimize!(model)
knock_out_flux =  MOI.get.(model, MOI.VariablePrimal(), x)
# @show knock_out_flux
@show MOI.get(model, MOI.ObjectiveValue())

println("")
println("loopless MOMA")
println("----------------------------------")
# MOMA
# modify objective
L = - reference_flux
Q = I(length(x))
@objective(model, Min, 1/2 * x' * Q * x + L' * x)
# @show objective_function(model)
optimize!(model)
@show MOI.get(model, MOI.ObjectiveValue())
