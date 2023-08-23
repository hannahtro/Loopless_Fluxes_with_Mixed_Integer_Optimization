using COBREXA, Serialization, COBREXA.Everything
using DataFrames, CSV, JSON
using SCIP, JuMP

include("../src/loopless_constraints.jl")
include("optimization_model.jl")

# FBA data 
function cobrexa_fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, json=true)
    # build model
    molecular_model = deserialize("../molecular_models/" * organism * ".js")
    # print_model(molecular_model, organism)

    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    solved_model = flux_balance_analysis(molecular_model, optimizer)
    model = solved_model.result
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
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
    end

    df = DataFrame(
        objective_value=primal_objective_value, 
        dual_bound=dual_objective_value,
        solution=[solution], 
        time=solved_time, 
        termination=status,
        time_limit=time_limit)

    type = "cobrexa_fba"
    if json 
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end

    return primal_objective_value, solution, status
end 

# loopless FBA data
function cobrexa_loopless_fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, json=true)
    # build model
    molecular_model = deserialize("../molecular_models/" * organism * ".js")
    # print_model(molecular_model, organism)

    # TODO: set time limit
    loopless_flux = flux_balance_analysis(
        molecular_model,
        optimizer,
        modifications = [add_loopless_constraints(), modify_optimizer_attribute(MOI.Silent(), true), modify_optimizer_attribute(MOI.TimeLimitSec(), time_limit)]
    )

    model = loopless_flux.result
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
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = primal_objective_value
    dict[:dual_bound] = dual_objective_value
    dict[:solution] = solution
    dict[:time] = solved_time
    dict[:termination] = status
    dict[:nodes] = nodes
    dict[:time_limit] = time_limit

    type = "cobrexa_loopless_fba"
    file_name = joinpath(@__DIR__,"../json/" * organism * "_" * type * "_" * string(time_limit) * ".json")

    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
    
    return primal_objective_value, solution, status
end 

