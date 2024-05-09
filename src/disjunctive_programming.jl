using DisjunctiveProgramming
using HiGHS
using COBREXA
using JSON

include("optimization_model.jl")

function add_disjunctive_loopless_constraints(model, S, internal_rxn_idxs, max_flux_bound; constraints_only_on=NaN)
    m, n = size(S)
    S_int = Array(S[:, internal_rxn_idxs])

    # define variables
    x = model[:x]
    G = @variable(model, G[1:length(internal_rxn_idxs)]) # approx ΔG for internal reactions
    μ = @variable(model, μ[1:size(S)[1]])
    @variable(model, Y[1:2*length(internal_rxn_idxs)], Logical) # Boolean variables

    # variable bounds
    for i in 1:length(internal_rxn_idxs)
        set_lower_bound(G[i], -10000)
        set_upper_bound(G[i], 10000)
    end

    for i in 1:m 
        set_lower_bound(μ[i], -10000)
        set_upper_bound(μ[i], 10000)
    end

    # loopless constraints
    if isnan(constraints_only_on)
        constraint(model, G' .== μ' * S_int)
        for (cidx, ridx) in enumerate(internal_rxn_idxs)
            @constraint(model, x[ridx] >= 0, Disjunct(Y[cidx]))
            @constraint(model, -max_flux_bound <= G[cidx] <= -1, Disjunct(Y[cidx]))
            @constraint(model, x[ridx] <= 0, Disjunct(Y[cidx+length(internal_rxn_idxs)]))
            @constraint(model, 1 <= G[cidx] <= max_flux_bound, Disjunct(Y[cidx+length(internal_rxn_idxs)]))
            @disjunction(model, [Y[cidx], Y[cidx+length(internal_rxn_idxs)]])
        end
    else
        # todo: only on index
        cidx, ridx = [(cidx, ridx) for (cidx, ridx) in enumerate(internal_rxn_idxs) if cidx==constraints_only_on][1]
        @constraint(model, x[ridx] >= 0, Disjunct(Y[cidx]))
        @constraint(model, -max_flux_bound <= G[cidx] <= -1, Disjunct(Y[cidx]))
        @constraint(model, x[ridx] <= 0, Disjunct(Y[cidx+length(internal_rxn_idxs)]))
        @constraint(model, 1 <= G[cidx] <= max_flux_bound, Disjunct(Y[cidx+length(internal_rxn_idxs)]))
        @disjunction(model, [Y[cidx], Y[cidx+length(internal_rxn_idxs)]])
    end
end

"""
optimize model and print process,
returns objective value, solution, time taken and termination status
"""
function optimize_gdp_model(model; time_limit=Inf, print_objective=false, silent=true, mute=true, gdp_method="BigM", big_m_constant=1000)
    if !mute 
        println("")
        println(type)
        println("----------------------------------")
    end
    if print_objective
        println("objective function : ", objective_function(model))
    end
    if !isinf(time_limit)
        set_time_limit_sec(model, time_limit)
    end
    if silent
        set_attribute(model, MOI.Silent(), true)
    else 
        set_attribute(model, MOI.Silent(), false)
    end

    if gdp_method == "BigM"
        optimize!(model, gdp_method = BigM(big_m_constant, true))
    elseif gdp_method == "Hull"
        optimize!(model, gdp_method = Hull())
    elseif gdp_method == "Indicator"
        optimize!(model, gdp_method = Indicator())
    end    
    
    status = termination_status(model)
    time = solve_time(model)

    # @show solution_summary(model)
    # @show status
    if has_values(model)
        primal_objective_value = MOI.get(model, MOI.ObjectiveValue())
        dual_objective_value = MOI.get(model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=2))
            println("")
        end
        solution = [value(var) for var in all_variables(model)]
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
    end

    return primal_objective_value, dual_objective_value, solution, time, status
end

function is_basic_feasible_solution(solution_x, solution_G, internal_rxn_idxs; S, lb, ub)
    for i in 1:length(internal_rxn_idxs)
        m = build_auxiliary_qp(i, solution_x, solution_G, S, lb, ub, internal_rxn_idxs)
        solution_qp, dual_objective_value, solution, time_taken, status = optimize_gdp_model(m; gdp_method="Hull", silent=true)
        _, ridx = [(cidx, ridx) for (cidx, ridx) in enumerate(internal_rxn_idxs) if cidx==i][1]
        @show solution_qp
        @show value(m[:x][ridx])
        @show value(m[:G][i])
        if !isapprox(solution_qp, 0, atol=1e-6)
            return true
        end
    end
    return false
end

function build_auxiliary_qp(i, solution_x, solution_G, S, lb, ub, internal_rxn_idxs; optimizer=SCIP.Optimizer)    
    m, n = size(S)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))

    model = GDPModel(optimizer)
    @variable(model, x[1:n])
    @constraint(model, mb, S * x .== 0) # mass balance
    for i in 1:n
        set_lower_bound(x[i], lb[i])
        set_upper_bound(x[i], ub[i])
    end    
    MOI.set(model, MOI.RelativeGapTolerance(), 1e-6)
    MOI.set(model, MOI.AbsoluteGapTolerance(), 1e-6)

    add_disjunctive_loopless_constraints(model, S, internal_rxn_idxs, max_flux_bound, constraints_only_on=i)

    cidx, ridx = [(cidx, ridx) for (cidx, ridx) in enumerate(internal_rxn_idxs) if cidx==i][1]
    @objective(model, Min, sum((x - solution_x).^2))

    return model
end