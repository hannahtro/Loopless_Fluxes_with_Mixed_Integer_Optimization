using SCIP, GLPK
using COBREXA, Serialization, COBREXA.Everything

model = deserialize("data/ec_e_coli_core.js")

m_rids = deserialize("data/metabolic_rids.js")
t_rids = deserialize("data/transport_rids.js")

m_gids = deserialize("data/metabolic_gids.js")
t_gids = deserialize("data/transport_gids.js")

# simplified enzyme constrained model (SMOMENT)
# adjust the metabolic activity within the cell to respect known enzymatic parameters and enzyme mass constraints
smodel = make_simplified_enzyme_constrained_model(
    model,
    reaction_mass_groups = Dict(
        "cytosol" => m_rids,
        "membrane" => t_rids,
    ),
    # total mass limit of each group of reactions
    reaction_mass_group_bounds = Dict(
        "cytosol" => 75.0,
        "membrane" => 75.0,
    ),
)

sol = flux_balance_analysis(smodel, GLPK.Optimizer)

fluxes = values_dict(:reaction, sol)
fluxes["BIOMASS_Ecoli_core_w_GAM"]
@show flux_summary(fluxes)
enzs = values_dict(:enzyme, sol)
# @show enzs

# Full enzyme constrained model (GECKO formulation)
ecmodel = make_enzyme_constrained_model(
    model,
    gene_product_mass_groups = Dict(
        "cytosol" => m_gids,
        "membrane" => t_gids,
    ),
    gene_product_mass_group_bounds = Dict(
        "cytosol" => 75.0,
        "membrane" => 75.0,
    ),
)

sol = flux_balance_analysis(ecmodel, GLPK.Optimizer)

fluxes = values_dict(:reaction, sol)
fluxes["BIOMASS_Ecoli_core_w_GAM"]
@show flux_summary(fluxes)
enzs = values_dict(:enzyme, sol)
# @show enzs

# LoadError: KeyError: key "PFK#forward#1" not found
# flux_solution = flux_balance_analysis(
#     ecmodel,
#     GLPK.Optimizer;
#     modifications = [add_loopless_constraints()],
# ) #|> values_dict
# @show flux_solution