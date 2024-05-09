using JSON
using SCIP

include("../src/loopless_constraints.jl")
include("../src/cuts_decomposition.jl")

seed = 10 
organism = "iJR904"
println("--------------------------------------------------------")
println("TEST " * organism)
println("--------------------------------------------------------")

# load model and internal reaction indices
molecular_model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
gecko_model = build_gecko_model(molecular_model, seed)
model = make_optimization_model(gecko_model, SCIP.Optimizer)
@show objective_function(model)
@show objective_sense(model)

S = stoichiometry(gecko_model)
lb, ub = bounds(gecko_model)
# model = build_fba_model(S, lb, ub, max_reactions=max_reactions)
max_flux_bound = maximum(abs.(vcat(lb, ub)))
m, num_reactions = size(S)
internal_rxn_idxs = internal_reactions(gecko_model, enzyme_data=true)

# load ll FBA solution 
println("")
println("ll FBA")
println("--------------------------------------------------------")
file_name = "loopless_fba_gecko_" * string(seed) * "_7200"
dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))
flux = dict["x"]
direction = dict["a"]
@show flux[381]

# check thermodynamic feasibility of all internal reactions 
solution_dict = Dict(:G => dict["G"], :μ => dict["μ"])
file_name="report_ll_fba"
feasible, _, _ = thermo_feasible_mu(internal_rxn_idxs, direction, S, max_flux_bound; scip_tol=0.001, solution_dict=solution_dict, file_name=file_name)
@show feasible
@show is_feasible(model, flux, direction, S, internal_rxn_idxs, [], lb, ub, tol=0.001, check_cuts=false, check_thermodynamic_feasibility=true, check_indicator=false)

reported_violations = filesize(file_name)
@assert reported_violations == 0

# load cb solution 
println("")
println("CB big M")
println("--------------------------------------------------------")
file_name = "combinatorial_benders_gecko_" * string(seed) * "_fast_big_m_7200"
dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))
flux = dict["x"]
direction = dict["a"]
@show flux[381]

# check thermodynamic feasibility of all internal reactions 
# solution_dict = Dict(:G => dict["G"], :μ => dict["μ"])
# file_name="report_cb_bigm"
feasible, sol_G, sol_μ = thermo_feasible_mu(internal_rxn_idxs, direction, S; scip_tol=1e-6)
@show feasible
@show is_feasible(model, flux, direction, S, internal_rxn_idxs, [], lb, ub, tol=0.001, check_cuts=false, check_thermodynamic_feasibility=true, check_indicator=false)

# insert CB solution into ll FBA
add_loopless_constraints(gecko_model, model, max_flux_bound, nullspace_formulation=false, enzyme_data=true)
solution_dict = Dict(:G => sol_G, :μ => sol_μ)
G = model[:G]
μ = model[:μ]
x = model[:x]
a = model[:a]
key_vector = vcat(G, μ, x, a)
value_vector = vcat(solution_dict[:G], solution_dict[:μ], flux, direction)
point = Dict(key_vector .=> value_vector)
report = primal_feasibility_report(model, point, atol=0.001)
file_name="report_cb_bigm"
writedlm(file_name, report)
reported_violations = filesize(file_name)
@assert reported_violations == 0

# # check thermodynamic feasibility of non zero internal reactions 
# non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-3)], internal_rxn_idxs)
# reaction_mapping = Dict()
# for (idx, val) in enumerate(internal_rxn_idxs)
#     reaction_mapping[val] = idx
# end

# G = dict["G"][collect(reaction_mapping[val] for val in non_zero_flux_indices)]
# μ = dict["μ"]
# solution_dict = Dict(:G => G, :μ => μ)
# non_zero_flux_directions = direction[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
# feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001, solution_dict=solution_dict)
# @show feasible
