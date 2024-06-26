using SCIP, JuMP
using COBREXA, Serialization

include("loopless_constraints.jl")
include("optimization_model.jl")

organism = "iAF692"
# organism = "iJR904"
# organism = "iML1515"

optimizer = SCIP.Optimizer
molecular_model = load_model("../molecular_models/" * organism * ".json")
# print_model(molecular_model, organism)

S = stoichiometry(molecular_model)
internal_rxn_idxs = internal_reactions(molecular_model)

model = make_optimization_model(molecular_model, optimizer)
add_loopless_constraints_mu(model, S, internal_rxn_idxs)

# print(model)

objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model)
