using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra
using Boscia, FrankWolfe
using DataFrames, CSV

include("functions.jl")

function json_to_js(organism="iML1515")
    # could not extract json.gz
    molecular_model = load_model(ObjectModel, "data/" * organism * ".json")
    serialize("data/" * organism * ".js", molecular_model)
end

function get_fba_data(organism="iML1515"; time_limit = 1800)
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = deserialize("data/" * organism * ".js")
    print_model(molecular_model)

    model = make_optimization_model(molecular_model, optimizer)
    @show model
    set_attribute(model, MOI.Silent(), true)

    # FBA
    type = "fba"
    objective_fba, _, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)

    df = DataFrame(
        objective_fba=objective_fba, 
        vars_fba=vars_fba, 
        time_fba=time_fba, 
        termination_fba=termination_fba)

    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end

    # loopless FBA
    type = "loopless_fba"
    add_loopless_constraints(molecular_model, model)
    @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, _, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, "loopless FBA", time_limit=time_limit)

    df = DataFrame(
        objective_loopless_fba=objective_loopless_fba, 
        vars_loopless_fba=vars_loopless_fba, 
        time_loopless_fba=time_loopless_fba, 
        termination_loopless_fba=termination_loopless_fba)

    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end
end

function get_moma_data(organism="iML1515", idx=1; var=NaN, time_limit=Inf, time_limit_fba=1800)
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = deserialize("data/" * organism * ".js")
    print_model(molecular_model)

    model = make_optimization_model(molecular_model, optimizer)
    @show model
    add_loopless_constraints(molecular_model, model)
    @show model

    # block reaction
    if isnan(var)
        var = model[:x]
    end

    set_lower_bound(var[idx], 0)
    set_upper_bound(var[idx], 0)

    # load from csv
    reference_flux = DataFrame(CSV.File(
        joinpath(@__DIR__, "csv/" * organism * "_loopless_fba_" * string(time_limit_fba) * ".csv")
        ))[!,:vars_loopless_fba][1:length(var)]
    primal_objective_value, solution, time, status = moma(model, var, reference_flux, time_limit=time_limit)

    df = DataFrame(
        objective_loopless_moma=primal_objective_value, 
        vars_loopless_moma=solution, 
        time_loopless_moma=time, 
        termination_loopless_moma=status)

    type = "moma"
    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end
end

function get_moma_boscia_data(organism="iML1515", idx=1; var=NaN, time_limit=Inf, time_limit_fba=1800)
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = deserialize("data/" * organism * ".js")
    print_model(molecular_model)

    model = make_optimization_model(molecular_model, optimizer)
    @show model
    add_loopless_constraints(molecular_model, model)
    set_attribute(model, MOI.Silent(), true)

    # block reaction
    if isnan(var)
        var = model[:x]
    end

    set_lower_bound(var[idx], 0)
    set_upper_bound(var[idx], 0)

    # load from csv
    reference_flux = DataFrame(CSV.File(
        joinpath(@__DIR__, "csv/" * organism * "_loopless_fba_" * string(time_limit_fba) * ".csv")
        ))[!,:vars_loopless_fba][1:length(var)]
    primal_objective_value, solution, time, status = moma_boscia(model, var, reference_flux, time_limit=time_limit)

    df = DataFrame(
        objective_loopless_moma=primal_objective_value, 
        vars_loopless_moma=solution, 
        time_loopless_moma=time, 
        termination_loopless_moma=status)

    type = "moma_boscia"
    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end

end