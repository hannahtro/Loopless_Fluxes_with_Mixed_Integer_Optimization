using COBREXA, Serialization
using LinearAlgebra, SparseArrays
using Boscia, FrankWolfe
using GLPK, JuMP
using DelimitedFiles

include("optimization_model.jl")

"""
compute internal reactions of COBREXA model
"""
function add_loopless_constraints(molecular_model, model, max_flux_bound=1000; nullspace_formulation=true, reduced=false)
    # loopless model
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(reactions(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    if nullspace_formulation
        add_loopless_constraints(model, stoichiometry(molecular_model), internal_rxn_idxs, max_flux_bound)
    else 
        if !reduced
            add_loopless_constraints_mu(model, stoichiometry(molecular_model), internal_rxn_idxs, max_flux_bound)
        else 
            add_loopless_constraints_mu_reduced(model, stoichiometry(molecular_model), internal_rxn_idxs, max_flux_bound)
        end
    end
end

"""
add loopless FBA constraints
"""
function add_loopless_constraints(model, S, internal_rxn_idxs::Vector{Int64}, max_flux_bound=1000)
    # @show length(internal_rxn_idxs)
    N_int = nullspace(Array(S[:, internal_rxn_idxs])) # no sparse nullspace function

    x = model[:x]
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

    # @show internal_rxn_idxs[1:10]
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        @constraint(model, -max_flux_bound * (1 - a[cidx]) <= x[ridx])
        @constraint(model, x[ridx] <= max_flux_bound * a[cidx])

        @constraint(
            model,
            -max_flux_bound * a[cidx] + (1 - a[cidx]) <= G[cidx]
        )
        @constraint(
            model,
            G[cidx] <=
            -a[cidx] + max_flux_bound * (1 - a[cidx])
        )
    end
    @constraint(model, N_int' * G .== 0)
end

"""
add loopless FBA constraints witout nullspace formulation
"""
function add_loopless_constraints_mu(model, S, internal_rxn_idxs::Vector{Int64}, max_flux_bound=1000)
    # @show length(internal_rxn_idxs)
    S_int = Array(S[:, internal_rxn_idxs])

    x = model[:x]
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    μ = @variable(model, μ[1:size(S)[1]])

    # @show internal_rxn_idxs[1:10]
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        @constraint(model, -max_flux_bound * (1 - a[cidx]) <= x[ridx])
        @constraint(model, x[ridx] <= max_flux_bound * a[cidx])

        @constraint(
            model,
            -max_flux_bound * a[cidx] + (1 - a[cidx]) <= G[cidx]
        )
        @constraint(
            model,
            G[cidx] <=
            -a[cidx] + max_flux_bound * (1 - a[cidx])
        )
    end

    @constraint(model, G' .== μ' * S_int)
end

"""
add relaxed loopless FBA constraints
"""
function add_relaxed_loopless_constraints(model, S, internal_rxn_idxs::Vector{Int64}; nullspace_formulation=false)
    # @show length(internal_rxn_idxs)
    x = model[:x]
    if nullspace_formulation
        N_int = nullspace(Array(S[:, internal_rxn_idxs])) # no sparse nullspace function
        G = @variable(model, G[1:length(internal_rxn_idxs)])

        @constraint(model, N_int' * G .== 0)
    else 
        S_int = Array(S[:, internal_rxn_idxs])
        G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
        μ = @variable(model, μ[1:size(S)[1]])

        @constraint(model, G' .== μ' * S_int)
    end
end

"""
add loopless FBA constraints witout nullspace formulation using less decision variables,
as G is defined by mu and S_int, we do not need G explicitly
"""
function add_loopless_constraints_mu_reduced(model, S, internal_rxn_idxs::Vector{Int64}, max_flux_bound=1000)
    # @show length(internal_rxn_idxs)
    S_int = Array(S[:, internal_rxn_idxs])

    x = model[:x]
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    μ = @variable(model, μ[1:size(S)[1]])

    # @show internal_rxn_idxs[1:10]
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        @constraint(model, -max_flux_bound * (1 - a[cidx]) <= x[ridx])
        @constraint(model, x[ridx] <= max_flux_bound * a[cidx])

        @constraint(
            model,
            -max_flux_bound * a[cidx] + (1 - a[cidx]) <= (μ' * S_int)[cidx]
        )
        @constraint(
            model,
            (μ' * S_int)[cidx] <=
            -a[cidx] + max_flux_bound * (1 - a[cidx])
        )
    end

end

"""
add loopless FBA constraints using indicator variables instead of big M formulation
"""
function add_loopless_indicator_constraints(molecular_model, model, max_flux_bound=1000)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(reactions(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    N_int = nullspace(Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])) # no sparse nullspace function

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    a = @variable(model, a[1:length(internal_rxn_idxs)])

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - 0.000001 >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + 0.000001 <= 0})

        @constraint(model, a[cidx] => {-max_flux_bound <= G[cidx] <= -1})
        @constraint(model, !a[cidx] => {1 <= G[cidx] <= max_flux_bound})
    end

    @constraint(model, N_int' * G .== 0)
end

"""
add loopless FBA constraints using indicator variables instead of big M formulation, without nullspace formulation
"""
function add_loopless_indicator_constraints_mu(molecular_model, model, max_flux_bound=1000)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(reactions(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    S_int = Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    a = @variable(model, a[1:length(internal_rxn_idxs)])
    μ = @variable(model, μ[1:size(stoichiometry(molecular_model))[1]])

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - 0.000001 >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + 0.000001 <= 0})

        @constraint(model, a[cidx] => {-max_flux_bound <= G[cidx] <= -1})
        @constraint(model, !a[cidx] => {1 <= G[cidx] <= max_flux_bound})
    end

    @constraint(model, G' .== μ' * S_int)
end

"""
add loopless FBA constraints using bilinear constraints
"""
function add_loopless_bilinear_constraints(molecular_model, model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(reactions(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    S_int = Array(stoichiometry(molecular_model)[:, internal_rxn_idxs])

    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    μ = @variable(model, μ[1:size(stoichiometry(molecular_model))[1]])
    a = @variable(model, a[1:length(internal_rxn_idxs)])

    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] == 0})
        @constraint(model, !a[cidx] => {G[cidx] * x[ridx] <= -0.000001})
    end

    @constraint(model, G' .== μ' * S_int)
end

"""
add constraints to block detected cycles in loopless FBA model
using Boolean expresssions
"""
function block_cycle_constraint(optimization_model, unbounded_cycles, flux_directions, internal_rxn_idxs, S; vector_formulation=true, shortest_cycles=false, block_limit=10000, nullspace_formulation=false)
    internal_reactions = Dict(ridx => cidx for (cidx, ridx) in enumerate(internal_rxn_idxs))
    a = optimization_model[:a] 

    # filter infeasible cycles
    if nullspace_formulation
        infeasible = [thermo_feasible(cycle, flux_directions[idx], S) ? 0 : 1 for (idx, cycle) in enumerate(unbounded_cycles)]
    else 
        infeasible = [thermo_feasible_mu(cycle, flux_directions[idx], S) ? 0 : 1 for (idx, cycle) in enumerate(unbounded_cycles)]
    end 
    idx_infeasible = findall(x->x==1, infeasible)
    unbounded_cycles = unbounded_cycles[idx_infeasible]
    flux_directions = flux_directions[idx_infeasible]

    # sort by length
    if shortest_cycles 
        cycles_length = [length(c) for c in unbounded_cycles]
        # @show issorted(cycles_length, rev=false)
        cycles_length, unbounded_cycles, flux_directions = getindex.((cycles_length, unbounded_cycles, flux_directions), (sortperm(cycles_length),))    
        # @show issorted(cycles_length, rev=false)
    end 

    # subset of cycles
    unbounded_cycles = unbounded_cycles[1:min(block_limit, length(unbounded_cycles))]

    # cycle through forward arcs:  a1 ∧ a2 ∧ a3 = 1
    # to block cycle: ¬(a1 ∧ a2 ∧ a3) = ¬a1 ∨ ¬a2 ∨ ¬a3 >= 1
    # ¬a1 = 1 - a1 => 
    # 1-a1 ∨ 1-a2 ∨ 1-a3 >= 1 
    # -a1 ∨ -a2 ∨ -a3 >= 1-3 = -2
    # for backward arc: ¬(¬a1) = a1
    num_blocked_cycles = 0
    if vector_formulation
        for (idx, cycle) in enumerate(unbounded_cycles)
            if length(cycle) > 2
                num_blocked_cycles += 1
                # get correct a, because v > a
                cycle_vars = [internal_reactions[i] for i in cycle]
                # reactions that have to be negated
                sum_forward = sum([1 for dir in flux_directions[idx] if dir > 0])
                dir_coef = [dir > 0 ? -1 : 1 for dir in flux_directions[idx]]
                # @show length(cycle_vars), length(dir_coef)
                constraint_coef = Array(sparsevec(cycle_vars, dir_coef, length(a)))
                @constraint(optimization_model, constraint_coef' * a >= 1 - sum_forward)
            end
        end
        # open("../experiments/csv/model_vector.lp", "w") do f
        #     print(f, optimization_model)
        # end
    else 
        for (idx, cycle) in enumerate(unbounded_cycles)
            if length(cycle) > 2
            num_blocked_cycles += 1
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
        end
            # open("../experiments/csv/model_loop.lp", "w") do f
            #     print(f, optimization_model)
            # end
        # print(optimization_model)
    end
    return num_blocked_cycles
end


"""
compute thermodynamic feasibility for a given cycle using the nullspace formulation,
    where the flux direction is captured by binary values
"""
# TODO: cylce should just contain internal reactions
function thermo_feasible(cycle, flux_directions, S, max_flux_bound=1000)
    # warning if flux_directions not integral
    if sum([(i != 1 || 1 !=0) ? 0 : 1 for i in flux_directions]) != 0
        @warn "flux directions should be binary"
    end
    thermo_feasible_model = Model(SCIP.Optimizer)
    G = @variable(thermo_feasible_model, G[1:length(cycle)]) # approx ΔG for internal reactions

    # add G variables for each reaction in cycle
    # @show cycle
    for (idx, _) in enumerate(cycle)
        if isapprox(flux_directions[idx], 0, atol=0.0001)
            @constraint(thermo_feasible_model, -max_flux_bound <= G[idx] <= -1)
        elseif isapprox(flux_directions[idx], 1, atol=0.0001)
            @constraint(thermo_feasible_model, 1 <= G[idx] <= max_flux_bound)
        end
    end

    N_int = nullspace(Array(S[:, cycle])) # no sparse nullspace function
    # @show N_int, Array(S[:, cycle])
    @constraint(thermo_feasible_model, N_int' * G .== 0)

    # print(thermo_feasible_model)

    _, _, solution, _, status = optimize_model(thermo_feasible_model)
    # @show solution, N_int' * solution
    return status == MOI.OPTIMAL
end

"""
compute thermodynamic feasibility for a given cycle without using the nullspace formulation,
    where the flux direction is captured by binary values, 
    the reactions should be internal reactions only
"""
function thermo_feasible_mu(cycle, flux_directions, S, max_flux_bound=1000; scip_tol=1e-6, solution_dict=Dict())
    # warning if flux_directions not integral
    if sum([(i != 1 || 1 !=0) ? 0 : 1 for i in flux_directions]) != 0
        @warn "flux directions should be binary"
    end
    thermo_feasible_model = Model(SCIP.Optimizer)

    MOI.set(thermo_feasible_model, MOI.RawOptimizerAttribute("numerics/feastol"), scip_tol)
    
    S_int = S[:, cycle]
    @variable(thermo_feasible_model, G[1:length(cycle)]) # approx ΔG for internal reactions
    @variable(thermo_feasible_model, μ[1:size(S_int)[1]])

    # add G variables for each reaction in cycle
    for (idx, _) in enumerate(cycle)
        if isapprox(flux_directions[idx], 0, atol=0.0001)
            @constraint(thermo_feasible_model, -max_flux_bound <= G[idx] <= -1 + scip_tol)
        elseif isapprox(flux_directions[idx], 1, atol=0.0001)
            @constraint(thermo_feasible_model, 1 - scip_tol <= G[idx] <= max_flux_bound)
        end
    end

    @constraint(thermo_feasible_model, G' .== μ' * S_int)   

    # print(thermo_feasible_model)

    _, _, solution, _, status = optimize_model(thermo_feasible_model)
    # if status == MOI.OPTIMAL
    #     @show MOI.get.(thermo_feasible_model, MOI.VariablePrimal(), G)
    #     @show MOI.get.(thermo_feasible_model, MOI.VariablePrimal(), μ)
    # end
    
    if !isempty(solution_dict)
        key_vector = vcat(G, μ)
        value_vector = vcat(solution_dict[:G], solution_dict[:μ])
        point = Dict(key_vector .=> value_vector)
        report = primal_feasibility_report(thermo_feasible_model, point, atol=0.001)
        writedlm("report.txt", report)
    end 

    return status == MOI.OPTIMAL
end

"""
returns assignment of G and a for a given solution
"""
function determine_G(S, solution, internal_rxn_idxs, max_flux_bound=1000)
    @assert length(solution) == size(S)[2]
    steady_state =  isapprox.(S * solution[1:size(S)[2]],0, atol=0.0001)
    @assert steady_state == ones(size(S)[1])
    Gibbs_model = Model(SCIP.Optimizer)
    N_int = nullspace(Array(S[:, internal_rxn_idxs])) # no sparse nullspace function
    G = @variable(Gibbs_model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions

    for (idx,ridx) in enumerate(internal_rxn_idxs)
        if solution[ridx] > 0
            @constraint(Gibbs_model, -max_flux_bound <= G[idx] <= -1)
        elseif solution[ridx] < 0
            @constraint(Gibbs_model, 1 <= G[idx] <= max_flux_bound)
        end
    end

    @constraint(Gibbs_model, N_int' * G .== 0)
    _, _, G, _, status = optimize_model(Gibbs_model)

    a = [solution[ridx] > 0 ? 1 : 0 for (idx,ridx) in enumerate(internal_rxn_idxs)] # if s is zero, a can be zero or one

    @show length(G), length(a)
    if length(G)!=length(a)
        @show isnan(G)
        @error "no assignment found for G, invalid flux" 
    end

    return vcat(solution, G, a)
end

"""
returns assignment of G, mu and a for a given solution
"""
function determine_G_mu(S, solution, internal_rxn_idxs, max_flux_bound=1000)
    @assert length(solution) == size(S)[2]
    Gibbs_model = Model(SCIP.Optimizer)
    S_int = Array(S[:, internal_rxn_idxs])
    G = @variable(Gibbs_model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    μ = @variable(Gibbs_model, μ[1:size(S_int)[1]])

    @show length(G), length(solution)

    for (idx,ridx) in enumerate(internal_rxn_idxs)
        if solution[ridx] > 0
            @constraint(Gibbs_model, -max_flux_bound <= G[idx] <= -1)
        elseif solution[ridx] < 0
            @constraint(Gibbs_model, 1 <= G[idx] <= max_flux_bound)
        end
    end

    @constraint(Gibbs_model, G' .== μ' * S_int)   
     _, _, sol, _, status = optimize_model(Gibbs_model)
     @show status

     a = [solution[ridx] > 0 ? 1 : 0 for (idx,ridx) in enumerate(internal_rxn_idxs)] # if s is zero, a can be zero or one

     if length(sol)==1
         @show isnan(sol)
         @error "no assignment found for G, invalid flux" 
     end

    return vcat(solution, sol[1:length(internal_rxn_idxs)], a, sol[length(internal_rxn_idxs)+1:end]) # G, a, μ
end

"""
check loopless violation of a given solution and corresponding binary variables
for internal reactions
"""
function check_loopless_violation(flux, flux_directions, S, internal_rxn_idxs, max_flux_bound=1000)
    model = Model(SCIP.Optimizer)
    S_int = S[:, internal_rxn_idxs]
    N_int = nullspace(Array(S[:, internal_rxn_idxs])) # no sparse nullspace function
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    α = @variable(model, α[1:length(internal_rxn_idxs)])

    internal_fluxes = flux[internal_rxn_idxs]
    # add G variables for each reaction in cycle
    non_zero_reactions = []
    for (idx, direction) in enumerate(flux_directions)
        if !isapprox(internal_fluxes[idx], 0, atol=0.0001)
            push!(non_zero_reactions, idx)
            if isapprox(flux_directions[idx], 1, atol=0.0001)
                @constraint(model, -max_flux_bound <= G[idx] <= 0)
                @constraint(model, α[idx] <= - G[idx])
            elseif isapprox(flux_directions[idx], 0, atol=0.0001)
                @constraint(model, 0 <= G[idx] <= max_flux_bound)
                @constraint(model, α[idx] <= G[idx])
            end
        end
    end

    @constraint(model, α >= 0)       
    @constraint(model, N_int' * G .== 0)

    @objective(model, Max, ones(length(α[non_zero_reactions]))' * α[non_zero_reactions])

    primal_objective_value, _, solution, _, status = optimize_model(model, print_objective=true)
    return primal_objective_value, solution
end
