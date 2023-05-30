using DataFrames, CSV

# compute dual gap with time limit of loopless FBA
function loopless_fba_data(organism; time_limit=1800)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)

    type = "loopless_fba"
    add_loopless_constraints(molecular_model, model)
    @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, _, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit)

    df = DataFrame(
        objective_loopless_fba=objective_loopless_fba, 
        vars_loopless_fba=vars_loopless_fba, 
        time_loopless_fba=time_loopless_fba, 
        termination_loopless_fba=termination_loopless_fba)

    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

# compute dual gap with time limit of loopless FBA with indicators
function loopless_indicator_fba_data(organism; time_limit=1800)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    # loopless FBA

    type = "loopless_indicator_fba"
    add_loopless_indicator_constraints(molecular_model, model)
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

    CSV.write(file_name, df, append=false, writeheader=true)
end

# compute dual gap with time limit of loopless FBA with blocked cycles
# compute dual gap with time limit of loopless FBA with indicators with bocked cycles
 