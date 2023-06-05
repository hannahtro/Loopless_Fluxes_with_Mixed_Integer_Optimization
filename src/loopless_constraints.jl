using COBREXA, Serialization, COBREXA.Everything
import COBREXA.Everything: add_loopless_constraints
using LinearAlgebra, SparseArrays
using Boscia, FrankWolfe

"""
compute internal reactions of COBREXA model
"""
function add_loopless_constraints(molecular_model, model)
    # loopless model
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    add_loopless_constraints(model, stoichiometry(molecular_model), internal_rxn_idxs)
end

"""
add loopless FBA constraints
"""
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

"""
add loopless FBA constraints using indicator variables instead of big M formulation
"""
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

"""
add constraints to block detected cycles in loopless FBA model
using Boolean expresssions
"""
function block_cycle_constraint(optimization_model, unbounded_cycles, flux_directions, internal_rxn_idxs; vector_formulation=true)
    internal_reactions = Dict(ridx => cidx for (cidx, ridx) in enumerate(internal_rxn_idxs))
    a = optimization_model[:a] 

    # cycle through forward arcs:  a1 ∧ a2 ∧ a3 = 1
    # to block cycle: ¬(a1 ∧ a2 ∧ a3) = ¬a1 ∨ ¬a2 ∨ ¬a3 >= 1
    # ¬a1 = 1 - a1 => 
    # 1-a1 ∨ 1-a2 ∨ 1-a3 >= 1 
    # -a1 ∨ -a2 ∨ -a3 >= 1-3 = -2
    # for backward arc: ¬(¬a1) = a1
    if vector_formulation
        for (idx, cycle) in enumerate(unbounded_cycles)
            # get correct a, because v > a
            cycle_vars = [internal_reactions[i] for i in cycle]
            # reactions that have to be negated
            sum_forward = sum([1 for dir in flux_directions[idx] if dir > 0])
            dir_coef = [dir > 0 ? -1 : 1 for dir in flux_directions[idx]]
            # @show length(cycle_vars), length(dir_coef)
            constraint_coef = Array(sparsevec(cycle_vars, dir_coef, length(a)))
            @constraint(optimization_model, constraint_coef' * a >= 1 - sum_forward)
        end
        # open("../csv/model_vector.lp", "w") do f
        #     print(f, optimization_model)
        # end
    else 
        for (idx, cycle) in enumerate(unbounded_cycles)
            # @show cycle # reactions
            cycle_vars = [a[internal_reactions[i]] for i in cycle]
            bool_blocked_cycle = []
            for (dir_idx, dir) in enumerate(flux_directions[idx])
                if dir > 0
                    push!(bool_blocked_cycle, 1)
                    push!(bool_blocked_cycle, -cycle_vars[dir_idx])
                elseif dir < 0
                    push!(bool_blocked_cycle, cycle_vars[dir_idx])
                end
            end
            # @show bool_blocked_cycle
            @constraint(optimization_model, sum(bool_blocked_cycle) >= 1)
        end
        # open("../csv/model_loop.lp", "w") do f
        #     print(f, optimization_model)
        # end
    end
    # print(optimization_model)
end

function thermo_feasible(unbounded_cycles_original, flux_directions, S)
    #TODO: verify for list of cycles
    for (cycle_idx, cycle) in enumerate(unbounded_cycles_original)
        thermo_feasible_model = Model(SCIP.Optimizer)
        G = @variable(thermo_feasible_model, G[1:length(cycle)]) # approx ΔG for internal reactions

        # add G variables for each reaction in cycle
        for (idx, ridx) in enumerate(cycle)
            # @show idx, ridx
            # @show flux_directions[cycle_idx][idx]
            if flux_directions[cycle_idx][idx] > 0
                @constraint(thermo_feasible_model, -1000 <= G[idx] <= -1)
            elseif flux_directions[cycle_idx][idx] < 0
                @constraint(thermo_feasible_model, 1 <= G[idx] <= 1000)
            end
        end

        # @show Array(S[:, cycle])
        N_int = nullspace(Array(S[:, cycle])) # no sparse nullspace function
        # @show N_int
        @constraint(thermo_feasible_model, N_int' * G .== 0)

        # print(thermo_feasible_model)

        _, _, solution, _, status = optimize_model(thermo_feasible_model)
        # @show solution
        # @show N_int' * solution
        return status == MOI.OPTIMAL
    end
end