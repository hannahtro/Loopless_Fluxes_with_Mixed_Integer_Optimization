using COBREXA, Serialization, COBREXA.Everything
using Serialization
using DataFrames, CSV
using SCIP, JuMP

include("loopless_constraints.jl")
include("optimization_model.jl")
include("cycle_detection.jl")
include("loopless_fba.jl")

"""
set solution vector as primal
"""
function variable_mapping(model::Model, solution)
    # Store a mapping of the variable primal solution
    @assert length(all_variables(model)) == length(solution)
    variable_primal = Dict(x => solution[idx] for (idx,x) in enumerate(all_variables(model)))
    # Now we can loop through our cached solutions and set the starting values.
    for (x, primal_start) in variable_primal
        set_start_value(x, primal_start)
    end
end

"""
compute loopless fba after setting primal for a given model
"""
function loopless_fba_set_primal(organism, model, S, internal_rxn_idxs; nullspace_formulation=true, flux=[], time_limit=600)
    if nullspace_formulation 
        solution = determine_G(S, flux, internal_rxn_idxs)   
    else
        solution = determine_G_mu(S, flux, internal_rxn_idxs)   
    end

    if !nullspace_formulation
        @assert length(solution) == length(flux) + 2 * length(internal_rxn_idxs) + size(S)[1]
    else 
        @assert length(solution) == length(flux) + 2 * length(internal_rxn_idxs)
    end 

    if nullspace_formulation 
        add_loopless_constraints(model, S, internal_rxn_idxs)
    else 
        add_loopless_constraints_mu(model, S, internal_rxn_idxs)
    end 

    # add primal solution
    variable_mapping(model, solution)

    # @show model
    type = "set_primal"
    if nullspace_formulation
        type = type * "_nullspace"
    end

    set_attribute(model, MOI.Silent(), true)
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false)

    num_nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=num_nodes)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
    return objective_loopless_fba, time_loopless_fba, num_nodes
end

"""
compute loopless fba after setting primal for a given organism
"""
function loopless_fba_set_primal(organism; flux=[], load=true, nullspace_formulation=false, time_limit=180)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")

    print_model(molecular_model)
    # get values for G and a for solution 
    S = stoichiometry(molecular_model)
    num_m, num_reactions = size(S)
    internal_rxn_idxs = internal_reactions(molecular_model)

    type = "solution"
    if nullspace_formulation
        type = type * "_nullspace"
    end
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    if load 
        df = CSV.read(file_name, DataFrame)
        solution = df[!,:sol]
    else 
        if isempty(flux)
            objective, flux, _, _ = loopless_fba_data(organism, time_limit=1800, nullspace_formulation=nullspace_formulation, csv=false)
            if !nullspace_formulation
                @assert length(flux) == num_reactions + 2 * length(internal_rxn_idxs) + size(S)[1]
            else 
                @assert length(flux) == num_reactions + 2 * length(internal_rxn_idxs)
            end 
            flux = flux[1:num_reactions]
        end
        # TODO: can take long, but is LP
        if nullspace_formulation
            solution = determine_G(S, flux, internal_rxn_idxs)   
        else
            solution = determine_G_mu(S, flux, internal_rxn_idxs)   
        end

        df_sol = DataFrame(sol = solution)
        CSV.write(file_name, df_sol, append=false, writeheader=true)
    end
    # @show length(flux)
    # @show length(solution)

    if !nullspace_formulation
        @assert length(solution) == length(flux) + 2 * length(internal_rxn_idxs) + size(S)[1]
    else 
        @show length(solution), length(internal_rxn_idxs), size(S)[1]
        @assert length(solution) == length(flux) + 2 * length(internal_rxn_idxs)
    end 

    # loopless FBA 
    # add loopless constraints and block cycles
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    if nullspace_formulation 
        add_loopless_constraints(molecular_model, model, nullspace_formulation=true)
    else 
        add_loopless_constraints(molecular_model, model, nullspace_formulation=false)
    end 

    # add primal solution
    variable_mapping(model, solution)

    type = "set_primal"
    if nullspace_formulation
        type = type * "_nullspace"
    end

    # @show model
    set_attribute(model, MOI.Silent(), true)
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false)

    num_nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=num_nodes)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
    return objective_loopless_fba, time_loopless_fba, num_nodes
end


