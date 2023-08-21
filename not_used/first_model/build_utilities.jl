"""
Check if reaction has a gene reaction rule assigned to it.
"""
function has_reaction_grr(model, rid)
    grr = reaction_gene_associations(model, rid)
    if isnothing(grr) ||
       isempty(grr) ||
       isnothing(first(grr)) ||
       isempty(first(grr)) ||
       first(first(grr)) == "U"
        return false
    end
    return true
end

"""
Check if reaction is a membrane protein heuristically.
"""
function is_membrane_reaction(model, rid)
    rs = reaction_stoichiometry(model, rid)
    length(unique(last(split(x, "_")) for x in keys(rs))) > 1 && return true
    any(
        in.(["q8h2_c", "q8_c", "mql8_c", "mqn8_c", "2dmmql8_c", "2dmmq8_c"], Ref(keys(rs))),
    ) && return true

    return false
end

"""
If reaction spontaneous, remove all grrs for simplicity  
"""
function rm_spontaneous!(model; spont_rxn_genes = ["s0001"])
    for rid in reactions(model)
        if has_reaction_grr(model, rid)
            gas = model.reactions[rid].gene_associations
            idxs = findall(
                x -> any(in.(spont_rxn_genes, Ref(keys(x.gene_product_stoichiometry)))),
                gas,
            )
            if length(idxs) == length(gas)
                model.reactions[rid].gene_associations = nothing
            else
                deleteat!(model.reactions[rid].gene_associations, idxs)
            end
        end
    end
    delete!.(Ref(model.genes), spont_rxn_genes)
    return nothing
end

"""
Add a turnover number to `model` for reaction `rid` with `default_kcat` if `rid`
is not found in `rid_kcat`. 
"""
function add_kcat!(model, rid, rid_kcat, default_kcat)
    to_hour = 3600.0 * 1e-3 # scale
    gas = model.reactions[rid].gene_associations
    isnothing(gas) && return nothing
    for isozyme in gas
        isozyme.kcat_forward = get(rid_kcat, rid, default_kcat) * to_hour
        isozyme.kcat_backward = get(rid_kcat, rid * "_b", default_kcat) * to_hour
    end
    return nothing
end

"""
Add isozyme stoichiometry.
"""
function add_gene_association_stoichiometry!(model)

    #: protein stoich map, infer from uniprot
    mer_map = Dict{String,Float64}(
        "Homotetramer" => 4,
        "Homodimer" => 2,
        "Homotrimer" => 3,
        "Homohexamer" => 6,
        "Homopentamer" => 5,
        "Homodecamer" => 10,
        "Homooctamer" => 8,
        "Homoheptamer" => 7,
        "Homododecamer" => 12,
        "Homomonomer" => 1,
        "Monomer" => 1,
    )

    #: infer protein stoichiometry from uniprot annotations
    uniprot =
        JSON.parsefile(joinpath("uniprot.json"))
    for rid in reactions(model)
        !has_reaction_grr(model, rid) && continue
        gas = model.reactions[rid].gene_associations
        for isozyme in gas
            if length(isozyme.gene_product_stoichiometry) == 1 # only assign homomers using uniprot data
                gid = first(keys(isozyme.gene_product_stoichiometry))
                isozyme.gene_product_stoichiometry[gid] =
                    max(get(mer_map, last(get(uniprot, gid, ["", ""])), 1.0), 1.0)
            else # assume complexes have uni-stoichiometry for now
                for gid in keys(isozyme.gene_product_stoichiometry)
                    isozyme.gene_product_stoichiometry[gid] = 1.0
                end
            end
        end
    end

    #:fix complex stoichiometry, use ComplexPortal database
    complex =
        JSON.parsefile("complex.json")
    for rid in reactions(model)

        !has_reaction_grr(model, rid) && continue

        gas = model.reactions[rid].gene_associations
        length(first(gas).gene_product_stoichiometry) == 1 && continue # skip monomers, already done

        for iso in gas
            for iso_com in values(complex)
                iso_gids = iso.gene_product_stoichiometry
                int_len = length(intersect(keys(iso_com), keys(iso_gids)))
                if int_len == length(iso_gids) == length(iso_com)
                    for gid in keys(iso_gids)
                        iso_gids[gid] = max(iso_com[gid], 1.0)
                    end
                end
            end
        end

        # this is not perfect, but should work for the majority of isozymes
    end

    return nothing
end
