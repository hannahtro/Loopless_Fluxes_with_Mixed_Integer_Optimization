using JSON
using SCIP

include("../src/loopless_constraints.jl")
include("../src/cuts_decomposition.jl")

# organisms with thermo infeasible solution:
# iNF517, Ashbya_aceri, Starmerella_bombicola, Tortispora, all yHMPu5000* models

organism = "Eremothecium_sinecaudum"
# organism = "Starmerella_bombicola_JCM9596"
# organism = "Ashbya_aceri"
# file_name = "loopless_fba_Gurobi_1800"
# dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))

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

# # solution = dict["solution"][1:num_reactions]
# # non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-5)], internal_rxn_idxs)
# # # non_zero_flux_indices = [idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-4)]
# # non_zero_flux_directions = [solution[idx] >= 1e-4 ? 1 : 0 for (idx,val) in enumerate(non_zero_flux_indices)] # TODO: for loop correct ???
# # @show thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)

# flux = dict["x"]
# direction = dict["a"]

# direction = round.(direction, digits=5)
# @show length(flux), length(direction)
# @show thermo_feasible_mu(internal_rxn_idxs, direction, S; scip_tol=0.00001)

# load cb solution 
file_name = "combinatorial_benders_fast_Gurobi_36000"
dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))
flux = dict["x"]
direction = dict["a"]

solution_dict = Dict(:G => dict["G"], :μ => dict["μ"])
feasible = thermo_feasible_mu(internal_rxn_idxs, direction, S; scip_tol=0.001, solution_dict=solution_dict)
@show feasible

@show is_feasible(model, flux, direction, S, internal_rxn_idxs, [], lb, ub, tol=0.001, check_cuts=false, check_thermodynamic_feasibility=true, check_indicator=false)

non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-4)], internal_rxn_idxs)
reaction_mapping = Dict()
for (idx, val) in enumerate(internal_rxn_idxs)
    reaction_mapping[val] = idx
end

non_zero_flux_directions = direction[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
@show feasible