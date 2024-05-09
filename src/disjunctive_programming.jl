using DisjunctiveProgramming
using HiGHS
using COBREXA
import COBREXA: make_optimization_model
using JSON

include("optimization_model.jl")
include("../src/loopless_fba.jl")

function dp_data(organism="e_coli_core"; optimizer=HiGHS.Optimizer, time_limit=1800, silent=true, type="dp", json=true, gdp_method="BigM", verbose=true, big_m_constant=1000, tol=1e-6) 
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    if verbose
        print_model(molecular_model, organism)
    end
    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))

    model = make_optimization_model(molecular_model::MetabolicModel, optimizer, true)
    MOI.set(model, MOI.RelativeGapTolerance(), 1e-6)
    MOI.set(model, MOI.AbsoluteGapTolerance(), 1e-6)

    add_disjunctive_loopless_constraints(model, S, internal_rxn_idxs, max_flux_bound)

    if gdp_method == "BigM"
        primal_objective_value, dual_objective_value, solution, time_taken, status = optimize_gdp_model(model; silent=silent, time_limit=time_limit, big_m_constant=big_m_constant)
        type = type * "_BigM_" * string(big_m_constant)
    elseif gdp_method == "Hull"
        primal_objective_value, dual_objective_value, solution, time_taken, status = optimize_gdp_model(model; gdp_method="Hull", silent=silent, time_limit=time_limit)
        type = type * "_Hull"
    elseif gdp_method == "Indicator"
        primal_objective_value, dual_objective_value, solution, time_taken, status = optimize_gdp_model(model; gdp_method="Indicator", silent=silent, time_limit=time_limit)
        type = type * "_Indicator"
    else 
        @error "set gdp method correctly"
    end

    if status == MOI.OPTIMAL
        S = stoichiometry(molecular_model)
        x = value.(model[:x])
        G = value.(model[:G])
        μ = value.(model[:μ])
        
        # check thermodynamic feasibility of solution
        steady_state =  isapprox.(S * x, 0, atol=0.0001)
        @assert steady_state == ones(size(S)[1])
        
        # check thermodynamic feasibility of solution through non zero flux
        flux = value.(model[:x])
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
        non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
        @assert feasible
    else 
        feasible = false
        x = NaN
        G = NaN 
        μ = NaN
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = primal_objective_value
    dict[:dual_bound] = dual_objective_value
    dict[:x] = x
    dict[:G] = G
    dict[:μ] = μ
    dict[:time] = time_taken
    dict[:termination] = status
    dict[:time_limit] = time_limit
    dict[:thermo_feasible] = feasible
    dict[:max_flux_bound] = max_flux_bound
    dict[:objective_function] = objective(molecular_model)
    dict[:sense] = objective_sense(model)

    if optimizer != SCIP.Optimizer
        type = type * "_" * replace(string(optimizer), ".Optimizer"=>"")
    end
    
    file_name = joinpath("json/" * organism * "_" * type * "_" * string(time_limit) * ".json")
    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
end

function make_optimization_model(model::MetabolicModel, optimizer, gdp; sense=MAX_SENSE)
    precache!(model)

    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = GDPModel(optimizer)
    @variable(optimization_model, x[1:n])
    @objective(optimization_model, sense, objective(model)' * x)
    @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    for i in 1:n
        set_lower_bound(x[i], xl[i])
        set_upper_bound(x[i], xu[i])
    end

    return optimization_model
end

function add_disjunctive_loopless_constraints(model, S, internal_rxn_idxs, max_flux_bound)
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
    @constraint(model, G' .== μ' * S_int)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
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