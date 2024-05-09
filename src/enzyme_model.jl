using COBREXA
using Distributions
using Random 

"""
add enzyme data to JsonModel or StandardModel
"""
function build_gecko_model(model, seed=3, μ=1.2)
    Random.seed!(seed)
    lb, ub = bounds(model)
    M = maximum(abs.(vcat(lb, ub)))

    # remove all non-directional bounds on the model, except for ATPM!!
    # FBA models need bounds on Glucose and Oxygen which is not needed in enzyme models
    @assert "ATPM" in reactions(model)
    for rid in reactions(model)
        rid == "ATPM" && continue
        lb = model.reactions[rid].lb
        ub = model.reactions[rid].ub
        model.reactions[rid].lb = lb < 0 ? -M : 0.0
        model.reactions[rid].ub = ub > 0 ? M : 0.0
    end

    kcat() = 10^rand(Normal(μ,1)) # realistic for kcats with units of 1/s
    mw() = rand(Uniform(10, 100)) #  realistic for proteins with units of kDa

    # first create all the isozymes
    reaction_isozymes = Dict{String,Vector{Isozyme}}()
    for rid in reactions(model)
        isos = Isozyme[]
        grrs = reaction_gene_association(model, rid)
        isnothing(grrs) && continue # only assign isozymes to reactions involving genes    
        for grr in grrs
            kfor = kcat() * 3600 # change unit to 1/h
            kback = kcat() * 3600
            push!(isos, Isozyme(Dict(k => 1 for k in grr), kfor, kback))
        end
        reaction_isozymes[rid] = isos
    end

    gene_product_mass_group = Dict{String,String}()
    gene_product_molar_mass = Dict{String,Float64}()
    gene_product_bounds = Dict{String,Tuple{Float64,Float64}}()
    for gid in genes(model)
        gene_product_mass_group[gid] = "capacity_bound_"*(rand() < 0.5 ? "A" : "B")
        gene_product_molar_mass[gid] = mw()
        gene_product_bounds[gid] = (0.0, 1000.0) # unit in mmol/gDW
    end
    gene_product_mass_group_bound = Dict("capacity_bound_A" => 0.5, "capacity_bound_B" => 0.5) # units g/gDW

    gmodel = make_gecko_model(
        model;
        reaction_isozymes,
        gene_product_bounds,
        gene_product_mass_group,
        gene_product_molar_mass,
        gene_product_mass_group_bound,
    )

    return gmodel
end 
