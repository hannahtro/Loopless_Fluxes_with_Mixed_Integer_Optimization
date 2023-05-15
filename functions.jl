using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe

function print_model(model, name="MODEL")
    println("")
    println(name)
    println("----------------------------------")
    println("number of metabolites : ", length(model.metabolites))
    println("number of reactions : ", length(model.reactions))
    println("number of genes : ", length(model.genes))
    # @show model.annotations
    # @show model.notes
    println("objective function: ", model.objective)
    # @show molecular_model.reactions
    # @show molecular_model.metabolites
    println("")
end

function optimize_model(model, type="FBA"; print_objective=false)
    println("")
    println(type)
    println("----------------------------------")
    if print_objective
        println("objective function : ", objective_function(model))
    end
    optimize!(model)
    println("objective value : ", round(MOI.get(model, MOI.ObjectiveValue()),digits=2))
    println("")
    return model
end

function add_loopless_constraints(molecular_model, model)
    # loopless model
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    # @show length(internal_rxn_idxs)

    N_int = nullspace(Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])) # no sparse nullspace function

    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        @constraint(model, -1000 * (1 - a[cidx]) <= x[ridx])
        @constraint(model, x[ridx] <= 1000 * a[cidx])

        @constraint(
            model,
            -1000 * a[cidx] + (1 - a[cidx]) <= G[cidx]
        )
        @constraint(
            model,
            G[cidx] <=
            -a[cidx] + 1000 * (1 - a[cidx])
        )
    end

    @constraint(model, N_int' * G .== 0)
end

function moma(model, x, reference_flux)
    L = - reference_flux
    Q = I(length(x))
    @objective(model, Min, 1/2 * x' * Q * x + L' * x)
    # @show objective_function(model)
    # @show typeof(objective_function(model))
    optimize_model(model, "loopless MOMA")
end

function moma_boscia(model, x, reference_flux, type="loopless MOMA in Boscia")
    println("")
    println(type)
    println("----------------------------------")
    
    L = - reference_flux
    Q = I(length(x))

    function f(x)
        length_var = size(Q)[1]
        f_x = 1/2 * x[1:length_var]' * Q * x[1:length_var] + L' * x[1:length_var]
    end
    # @show f(ones(length(x)))
    
    function grad!(storage, x)
        length_var = size(Q)[1]
        storage[1:length_var] = Q * x[1:length_var] + L
        storage[length_var+1:length(x)] .= 0
    end
    
    set_objective_sense(model, FEASIBILITY_SENSE)
    moi_model = backend(model)
    lmo = FrankWolfe.MathOptLMO(moi_model)
    x, _, result = Boscia.solve(f, grad!, lmo, verbose=true) 
    println("objective value : ", round(f(x),digits=2))   
    println("")
end