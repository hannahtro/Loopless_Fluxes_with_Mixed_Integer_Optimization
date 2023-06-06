using COBREXA, Serialization, COBREXA.Everything
using Serialization
using DataFrames, CSV
using SCIP, JuMP

include("../src/loopless_constraints.jl")
include("../src/optimization_model.jl")
include("../src/cycle_detection.jl")
include("loopless_fba.jl")

function variable_mapping(model::Model, solution)
    # Store a mapping of the variable primal solution
    variable_primal = Dict(x => solution[idx] for (idx,x) in enumerate(all_variables(model)))
    # Now we can loop through our cached solutions and set the starting values.
    for (x, primal_start) in variable_primal
        set_start_value(x, primal_start)
    end
end

function loopless_fba_set_primal(organism; load=true, time_limit=180)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")

    # get values for G and a for solution 
    S = stoichiometry(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)

    type = "solution"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    if load 
        df = CSV.read(file_name, DataFrame)
        solution = df[!,:sol]
    else 
        # TODO: can take long, but is LP
        solution = determine_G(S, solution, internal_rxn_idxs)   
        df_sol = DataFrame(sol = solution)
        CSV.write(file_name, df_sol, append=false, writeheader=true)
    end

    # loopless FBA 
    # add loopless constraints and block cycles
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    add_loopless_constraints(molecular_model, model)
    
    # add primal solution
    variable_mapping(model, solution)

    @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false)

    nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes)

    type = "set_primal"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

organism = "iAF692"

loopless_fba_set_primal(organism, time_limit=1800)
loopless_fba_data(organism, time_limit=1800)
