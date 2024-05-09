using COBREXA, Serialization
import COBREXA: add_loopless_constraints, change_optimizer_attribute
using DataFrames, CSV, JSON
using SCIP, JuMP

include("../src/loopless_constraints.jl")
include("optimization_model.jl")

# FBA data 
function cobrexa_fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, json=true)
    # build model
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    internal_rxn_idxs = internal_reactions(molecular_model)

    model = flux_balance_analysis(molecular_model, optimizer)
    status = termination_status(model)
    solved_time = solve_time(model)

    if has_values(model)
        primal_objective_value = MOI.get(model, MOI.ObjectiveValue())
        dual_objective_value = MOI.get(model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=5))
            println("")
        end
        solution = [value(var) for var in all_variables(model)]
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
        non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
        feasible = false
    end

    @show feasible
    dict = Dict{Symbol, Any}()

    dict[:objective_value] = primal_objective_value 
    dict[:dual_bound] = dual_objective_value
    dict[:solution] = solution
    dict[:thermo_feasible] = feasible
    dict[:time] = solved_time
    dict[:termination] = status
    dict[:time_limit] = time_limit

    type = "cobrexa_fba"
    file_name = "json/" * organism * "_" * type * "_" * string(time_limit) * ".json"

    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end

    return primal_objective_value, solution, status
end 

# loopless FBA data
function cobrexa_loopless_fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, json=true)
    # build model
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)
    S = stoichiometry(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)
    model = flux_balance_analysis(
        molecular_model,
        optimizer,
        modifications = [add_loopless_constraints, change_optimizer_attribute(MOI.Silent(), true), change_optimizer_attribute(MOI.TimeLimitSec(), time_limit)]
    )

    status = termination_status(model)
    solved_time = solve_time(model)
    nodes = MOI.get(model, MOI.NodeCount())

    if has_values(model)
        primal_objective_value = MOI.get(model, MOI.ObjectiveValue())
        dual_objective_value = MOI.get(model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=5))
            println("")
        end
        solution = [value(var) for var in all_variables(model)]
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
        non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
        @show feasible
        if status == MOI.OPTIMAL
            @assert feasible
        end 
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
        feasible = false
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = primal_objective_value
    dict[:dual_bound] = dual_objective_value
    dict[:solution] = solution
    dict[:thermo_feasible] = feasible
    dict[:time] = solved_time
    dict[:termination] = status
    dict[:nodes] = nodes
    dict[:time_limit] = time_limit

    type = "cobrexa_loopless_fba"
    file_name = "json/" * organism * "_" * type * "_" * string(time_limit) * ".json"

    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
    
    return primal_objective_value, solution, status
end 

