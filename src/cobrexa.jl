using COBREXA, Serialization, COBREXA.Everything
using DataFrames, CSV
using SCIP, JuMP

# FBA data 
function fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, csv=true)
    # build model
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

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
        solution=solution, 
        time=solved_time, 
        termination=status,
        time_limit=time_limit)

    type = "cobrexa_fba"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if csv 
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=false)    
        end
    end
end 

# loopless FBA data
function loopless_fba_data(organism; optimizer=SCIP.Optimizer, time_limit=1800, mute=true, csv=true)
    # build model
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    @infiltrate
    loopless_flux = flux_balance_analysis(
        molecular_model,
        optimizer,
        modifications = [add_loopless_constraints()]
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

    df = DataFrame(
        objective_value=primal_objective_value, 
        dual_bound=dual_objective_value,
        solution=solution, 
        time=solved_time, 
        termination=status,
        nodes=nodes,
        time_limit=time_limit)

    type = "cobrexa_loopless_fba"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if csv 
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=false)    
        end
    end
end 

