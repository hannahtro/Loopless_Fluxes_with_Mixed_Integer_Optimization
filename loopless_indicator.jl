using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra
using DataFrames, CSV

include("functions.jl")

function loopless_fba_data()
    # loopless FBA
    type = "loopless_fba"
    time_limit = 3600  #64000
    add_loopless_constraints(molecular_model, model)
    @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
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

function loopless_indicator_fba_data(model, molecular_model, organism)
    # loopless FBA
    type = "loopless_indicator_fba"
    time_limit = 600  #64000
    add_loopless_indicator_constraints(molecular_model, model)
    @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
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

organism = "iJR904"

# build model
optimizer = SCIP.Optimizer
molecular_model = deserialize("data/" * organism * ".js")
# print_model(molecular_model, organism)

model = make_optimization_model(molecular_model, optimizer)
# @show model
set_attribute(model, MOI.Silent(), true)

loopless_indicator_fba_data(model, molecular_model, organism)