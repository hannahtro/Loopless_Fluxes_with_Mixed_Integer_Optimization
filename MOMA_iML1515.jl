using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP, Tulip
using LinearAlgebra
using Boscia, FrankWolfe
using DataFrames, CSV

include("functions.jl")

@show VERSION

organism = "iML1515"

# could not extract json.gz
# molecular_model = load_model(ObjectModel, "data/iML1515.json")
# serialize("data/iML1515.js", molecular_model)

# build model
optimizer = SCIP.Optimizer

molecular_model = deserialize("data/" * organism * ".js")
print_model(molecular_model, "Escherichia coli str. K-12 substr. MG1655")

model = make_optimization_model(molecular_model, optimizer)
@show model
set_attribute(model, MOI.Silent(), true)

# FBA
type = "fba"
objective_fba, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)

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

