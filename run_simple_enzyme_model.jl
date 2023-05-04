using Tulip, COBREXA, Serialization, COBREXA.Everything
using Escher, CairoMakie, ColorSchemes

model = deserialize("ec_e_coli_core.js")

m_rids = deserialize("metabolic_rids.js")
t_rids = deserialize("transport_rids.js")

m_gids = deserialize("metabolic_gids.js")
t_gids = deserialize("transport_gids.js")

# simplified enzyme constrained model (SMOMENT)
smodel = make_simplified_enzyme_constrained_model(
    model,
    reaction_mass_groups = Dict(
        "cytosol" => m_rids,
        "membrane" => t_rids,
    ),
    reaction_mass_group_bounds = Dict(
        "cytosol" => 75.0,
        "membrane" => 75.0,
    ),
)

sol = flux_balance_analysis(smodel, Tulip.Optimizer)

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

sol = flux_balance_analysis(ecmodel, Tulip.Optimizer)

fluxes = values_dict(:reaction, sol)
fluxes["BIOMASS_Ecoli_core_w_GAM"]
@show flux_summary(fluxes)
enzs = values_dict(:enzyme, sol)
# @show enzs

f = Figure(resolution = (1200, 800));
ax = Axis(f[1, 1]);
###### PLOT FUNCTION
hidexdecorations!(ax)
hideydecorations!(ax)
escherplot!(ax, joinpath(pkgdir(Escher), "", "e_coli_core.json"))