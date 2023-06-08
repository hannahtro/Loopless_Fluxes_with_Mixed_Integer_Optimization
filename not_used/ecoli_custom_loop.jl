using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra
using Boscia, FrankWolfe
include("functions.jl")
include("max_flow.jl")

# build model
optimizer = SCIP.Optimizer
# optimizer = Tulip.Optimizer

molecular_model = deserialize("data/ec_e_coli_core.js")
S = stoichiometry(molecular_model)
@show size(stoichiometry(molecular_model))

print_model(molecular_model, "E COLI COBRA MODEL")

# add loop
r_1 = ReactionBidirectional("r_1", Dict("glc__D_e" => -1.0, "gln__L_c" => 1.0))
r_2 = ReactionBidirectional("r_2", Dict("gln__L_c" => -1.0, "gln__L_e" => 1.0))
r_3 = ReactionBidirectional("r_3", Dict("gln__L_e" => -1.0, "glc__D_e" => 1.0))

add_reactions!(molecular_model, [r_1, r_2, r_3])

model = make_optimization_model(molecular_model, optimizer)
@show model
set_attribute(model, MOI.Silent(), true)

# # FBA
# primal_objective_value, _, solution, time, status = optimize_model(model, print_objective=true)
# @show solution 

# FBA with modified objective
x = model[:x]
@objective(model, Max, sum(x))
primal_objective_value, _, solution, time, status = optimize_model(model)
# check if solution contains loops
#TODO: get stoichiometric matrix from model including bounds and objective
@show size(stoichiometry(molecular_model))
cycles, edge_mapping = ubounded_cycles(S_transform, solution)
@show cycles

unbounded_cycles, unbounded_cycles_original, flux_values = unbounded_cycles_S(cycles, edge_mapping, solution, reaction_mapping)
@show unbounded_cycles_original

# # loopless FBA
# restrict cycle in FBA
add_loopless_constraints(molecular_model, model)
# @show model
a = model[:a]
for cycle in unbounded_cycles_original
    cycle_vars = [a[i] for i in cycle]
    @show cycle_vars
    @constraint(model, sum(cycle_vars) >= 1)
end
# @show model

primal_objective_value, _, solution, time, status = optimize_model(model, "loopless FBA")
cycles, edge_mapping = ubounded_cycles(S_transform, solution)
@show cycles

# # loopless FBA
# add_loopless_constraints(molecular_model, model)
# optimize_model(model, "loopless FBA")

# # loopless FBA max x[25]
# x = model[:x]
# @objective(model, Max, x[25])
# optimize_model(model, "loopless FBA without blocked reaction", print_objective=true)
# reference_flux =  MOI.get.(model, MOI.VariablePrimal(), x)

# # FBA on blocked reaction
# # block reaction
# set_lower_bound(x[1], 0)
# set_upper_bound(x[1], 0)
# optimize_model(model, "loopless FBA with blocked reaction")
# # knock_out_flux =  MOI.get.(model, MOI.VariablePrimal(), x)

# MOMA
# @objective(model, Min, 1/2 * x' * Q * x + L' * x)
# # @show objective_function(model)
# # @show typeof(objective_function(model))
# optimize_model(model, "loopless MOMA")
# moma(model, x, reference_flux)
# unresolved numerical troubles in LP 1022 cannot be dealt with

# MOMA with Boscia 
# moma_boscia(model, x, reference_flux)
# objective value : -5.07598806e6

# FBA with enzyme constraints
# basic_model = deserialize("data/ec_e_coli_core.js")

# m_rids = deserialize("data/metabolic_rids.js")
# t_rids = deserialize("data/transport_rids.js")
# m_gids = deserialize("data/metabolic_gids.js")
# t_gids = deserialize("data/transport_gids.js")

# smodel = make_simplified_enzyme_constrained_model(
#     basic_model,
#     reaction_mass_groups = Dict(
#         "cytosol" => m_rids,
#         "membrane" => t_rids,
#     ),
#     # total mass limit of each group of reactions
#     reaction_mass_group_bounds = Dict(
#         "cytosol" => 75.0,
#         "membrane" => 75.0,
#     ),
# )

# optimizer = SCIP.Optimizer
# model = make_optimization_model(smodel, optimizer)
# set_attribute(model, MOI.Silent(), true)
# optimize_model(model, "enzyme-constrained FBA")

# # ecFBA with loops 
# r_1 = ReactionBidirectional("r_1", Dict("glc__D_e" => -1.0, "gln__L_c" => 1.0))
# r_2 = ReactionBidirectional("r_2", Dict("gln__L_c" => -1.0, "gln__L_e" => 1.0))
# r_3 = ReactionBidirectional("r_3", Dict("gln__L_e" => -1.0, "glc__D_e" => 1.0))

# add_reactions!(basic_model, [r_1, r_2, r_3])
# smodel = make_simplified_enzyme_constrained_model(
#     basic_model,
#     reaction_mass_groups = Dict(
#         "cytosol" => m_rids,
#         "membrane" => t_rids,
#     ),
#     # total mass limit of each group of reactions
#     reaction_mass_group_bounds = Dict(
#         "cytosol" => 75.0,
#         "membrane" => 75.0,
#     ),
# )

# optimizer = SCIP.Optimizer
# model = make_optimization_model(smodel, optimizer)
# set_attribute(model, MOI.Silent(), true)
# optimize_model(model, "enzyme-constrained FBA with loops")

# loopless ecFBA 
# add_loopless_constraints(smodel, model)
# optimize_model(model, "loopless ecFBA")
# DOES NOT WORK YET

# blocked loopless ecFBA

# MOMA loopless ecFBA