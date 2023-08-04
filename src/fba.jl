using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using DataFrames, CSV

include("optimization_model.jl")

function get_fba_data(organism="iML1515"; time_limit=1800, type = "fba", save_lp=false, csv=true)
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = deserialize("../data/" * organism * ".js")
    print_model(molecular_model)

    model = make_optimization_model(molecular_model, optimizer)
    @show model
    set_attribute(model, MOI.Silent(), true)

    if save_lp
        open("../csv/models/fba_model_" * organism * ".lp", "w") do f
            print(f, model)
        end
    end

    # FBA
    objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)

    df = DataFrame(
        objective_value=objective_fba, 
        dual_bound=dual_bound,
        solution=[vars_fba], 
        time=time_fba, 
        termination=termination_fba,
        time_limit=time_limit)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    if csv 
        if !isfile(file_name)
            CSV.write(file_name, df, append=true, writeheader=true)
        else 
            CSV.write(file_name, df, append=true)    
        end
    end

    # # loopless FBA
    # type = "loopless_fba"
    # add_loopless_constraints(molecular_model, model)
    # @show model
    # set_attribute(model, MOI.Silent(), false)
    # objective_loopless_fba, _, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
    #     optimize_model(model, "loopless FBA", time_limit=time_limit)

    # df = DataFrame(
    #     objective_loopless_fba=objective_loopless_fba, 
    #     vars_loopless_fba=vars_loopless_fba, 
    #     time_loopless_fba=time_loopless_fba, 
    #     termination_loopless_fba=termination_loopless_fba)

    # file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    # if !isfile(file_name)
    #     CSV.write(file_name, df, append=true, writeheader=true)
    # else 
    #     CSV.write(file_name, df, append=true)    
    # end
end