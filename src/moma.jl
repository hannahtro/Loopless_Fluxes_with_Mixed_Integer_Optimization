using COBREXA, Serialization, COBREXA.Everything
using SCIP, JuMP
using LinearAlgebra
using Boscia, FrankWolfe

"""
perform moma for given reference flux,
    reaction has to be blocked in the model before calling moma()
"""

function moma(model, x, reference_flux; time_limit)
    L = - reference_flux
    Q = I(length(x))
    @objective(model, Min, 1/2 * x' * Q * x + L' * x)
    # @show objective_function(model)
    # @show typeof(objective_function(model))
    primal_objective_value, dual_objective_value, solution, time, status = optimize_model(model, "loopless MOMA"; time_limit=time_limit)

    return primal_objective_value, solution, time, status
end

"""
use boscia to perform moma for given reference flux 
"""
function moma_boscia(model, x, reference_flux, type="loopless MOMA in Boscia"; time_limit=Inf)
    println("")
    println(type)
    println("----------------------------------")
    
    L = - reference_flux
    Q = I(length(x))

    function f(x)
        length_var = size(Q)[1]
        f_x = 1/2 * x[1:length_var]' * Q * x[1:length_var] + L' * x[1:length_var]
    end
    # @show f(ones(length(x)))
    
    function grad!(storage, x)
        length_var = size(Q)[1]
        storage[1:length_var] = Q * x[1:length_var] + L
        storage[length_var+1:length(x)] .= 0
    end
    
    set_time_limit_sec(model, time_limit)
    set_objective_sense(model, FEASIBILITY_SENSE)
    moi_model = backend(model)
    lmo = FrankWolfe.MathOptLMO(moi_model)
    x, _, result = Boscia.solve(f, grad!, lmo, verbose=true, time_limit=time_limit) 
    primal_objective_value = f(x)
    println("objective value : ", round(primal_objective_value, digits=2))   
    println("")

    return primal_objective_value, x, result[:total_time_in_sec], result[:status]
end

function get_moma_data(organism="iML1515", idx=1; var=NaN, time_limit=Inf, time_limit_fba=1800, type = "moma")
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = load_model("data/" * organism * ".json")
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

    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end
end

function get_moma_boscia_data(organism="iML1515", idx=1; var=NaN, time_limit=Inf, time_limit_fba=1800, type = "moma_boscia")
    # build model
    optimizer = SCIP.Optimizer

    molecular_model = load_model("data/" * organism * ".json")
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

    file_name = joinpath(@__DIR__,"csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if !isfile(file_name)
        CSV.write(file_name, df, append=true, writeheader=true)
    else 
        CSV.write(file_name, df, append=true)    
    end
end