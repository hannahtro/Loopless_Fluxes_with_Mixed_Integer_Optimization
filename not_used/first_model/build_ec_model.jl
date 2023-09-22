using COBREXA, JSON, CSV, Serialization

include("build_utilities.jl")

# Load base model
model = load_model(ObjectModel, "e_coli_core.json")

# add some transport costs
for rid in [
    "ACALDt", "ACt2r", "PYRt2", "SUCCt3", "ETOHt2r"
]
    model.reactions[rid].gene_associations = [Isozyme(["b2492"]), Isozyme(["b0904"])]
end

# remove spontaneous reaction genes
rm_spontaneous!(model; spont_rxn_genes = ["s0001"])

# classify reactions
metabolic_reactions = String[]
transport_reactions = String[]
other_reactions = String[]
pseudo_reactions = String[]
for rid in reactions(model)
    if is_metabolic_reaction(model, rid) && !is_membrane_reaction(model, rid)
        push!(metabolic_reactions, rid)
    elseif is_transport_reaction(model, rid)
        push!(transport_reactions, rid)
    elseif any(
        f(model, rid) for f in [
            is_atp_maintenance_reaction,
            is_biomass_reaction,
            is_boundary,
            is_exchange_reaction,
            is_pseudo_reaction,
            is_spontaneous_reaction,
        ]
    )
        push!(pseudo_reactions, rid) #  these get no kcat
    else
        push!(other_reactions, rid)
    end
end

# load kcat data
rid_kcat = Dict(
    String(row.react_id) => row.kappmax_KO_ALE_davidi_per_pp_per_s_ensemble_model for
    row in CSV.File("Dataset_S1C_turnonver_n.csv")
)

default_kcat = 65.0
for rid in [metabolic_reactions; transport_reactions; other_reactions]
    add_kcat!(model, rid, rid_kcat, default_kcat)
end

# load molar masses of gene products
unassigned_genes = String[]
gid_mm = JSON.parsefile("protein_masses.json")
for gid in genes(model)
    if haskey(gid_mm, gid)
        model.genes[gid].product_molar_mass = gid_mm[gid]
    else
        push!(unassigned_genes, gid)
    end
end

# load isozyme stoichiometry
add_gene_association_stoichiometry!(model)

# Fix constraints for exchange reactions
for rid in ["EX_glc__D_e"] # unconstrain
    model.reactions[rid].lower_bound = -1000.0
    model.reactions[rid].upper_bound = 1000.0
end

# get genes of metabolic and transporter reactions
transport_gids = String[]
metabolic_gids = String[]
for rid in reactions(model)
    !has_reaction_grr(model, rid) && continue
    if rid in metabolic_reactions
        for grr in reaction_gene_associations(model, rid)
            append!(metabolic_gids, grr)
        end
    elseif rid in [transport_reactions; other_reactions]
        for grr in reaction_gene_associations(model, rid)
            append!(transport_gids, grr)
        end
    end
end
unique!(transport_gids)
unique!(metabolic_gids)

setdiff(genes(model), [transport_gids; metabolic_gids])
intersect(transport_gids, metabolic_gids)

# fix some stuff
isos = model.reactions["ATPS4r"].gene_associations
for i in 1:length(isos)
    all(values(isos[1].gene_product_stoichiometry) .== 1) && deleteat!(model.reactions["ATPS4r"].gene_associations, 1)
end

for rid in reactions(model)
    isos = reaction_isozymes(model, rid)
    isnothing(isos) && continue
    for iso in isos 
        for (k, v) in iso.gene_product_stoichiometry
            if isnothing(model.genes[k].product_molar_mass)
                model.genes[k].product_molar_mass = 40.0 # 
            end
        end
    end 
end

# save model
serialize("ec_e_coli_core.js", model)

# save relevant IDs
serialize("metabolic_gids.js", metabolic_gids)
serialize("transport_gids.js", transport_gids)
serialize("metabolic_rids.js", metabolic_reactions)
serialize("transport_rids.js", [transport_reactions; other_reactions])

for rid in reactions(model)
    if isnothing(model.reactions[rid].gene_associations)
        println(rid)
    end
end