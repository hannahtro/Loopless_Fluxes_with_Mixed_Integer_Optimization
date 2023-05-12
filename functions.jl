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
    println("objective value : ", MOI.get(model, MOI.ObjectiveValue()))
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
