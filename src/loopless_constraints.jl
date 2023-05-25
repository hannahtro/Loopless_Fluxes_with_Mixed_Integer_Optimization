using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe

function add_loopless_constraints(molecular_model, model)
    # loopless model
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    add_loopless_constraints(model, stoichiometry(molecular_model), internal_rxn_idxs)
end

function add_loopless_constraints(model, S, internal_rxn_idxs::Vector{Int64})
    # @show length(internal_rxn_idxs)
    N_int = nullspace(Array(S[:, internal_rxn_idxs])) # no sparse nullspace function

    x = model[:x]
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

    # @show internal_rxn_idxs[1:10]
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

function add_loopless_indicator_constraints(molecular_model, model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    N_int = nullspace(Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])) # no sparse nullspace function

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    a = @variable(model, a[1:length(internal_rxn_idxs)])

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + eps() <= 0})

        @constraint(model, a[cidx] => {-1000 <= G[cidx] <= -1})
        @constraint(model, !a[cidx] => {1 <= G[cidx] <= 1000})
    end

    @constraint(model, N_int' * G .== 0)
end

