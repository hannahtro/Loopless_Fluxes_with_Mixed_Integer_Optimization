using JSON
using SCIP

include("../src/loopless_constraints.jl")

# organisms with thermo infeasible solution:
# iNF517, Ashbya_aceri, Starmerella_bombicola, Tortispora, all yHMPu5000* models

# organism = "Eremothecium_gossypii"
organism = "Starmerella_bombicola_JCM9596"
file_name = "loopless_fba_Gurobi_1800"
dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))

molecular_model = load_model("../molecular_models/ecModels/Classical/emodel_" * organism * "_classical.mat")
model = make_optimization_model(molecular_model, SCIP.Optimizer)
S = stoichiometry(molecular_model)
lb, ub = bounds(molecular_model)
# model = build_fba_model(S, lb, ub, max_reactions=max_reactions)
max_flux_bound = maximum(abs.(vcat(lb, ub)))
m, num_reactions = size(S)
internal_rxn_idxs = [
    ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
    !is_boundary(reaction_stoichiometry(molecular_model, rid))
]

# solution = dict["solution"][1:num_reactions]
# non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-5)], internal_rxn_idxs)
# # non_zero_flux_indices = [idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-4)]
# non_zero_flux_directions = [solution[idx] >= 1e-4 ? 1 : 0 for (idx,val) in enumerate(non_zero_flux_indices)] # TODO: for loop correct ???
# @show thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)

flux = dict["x"]
direction = dict["a"]

direction = round.(direction, digits=5)
@show length(flux), length(direction)
@show thermo_feasible_mu(internal_rxn_idxs, direction, S; scip_tol=0.00001)