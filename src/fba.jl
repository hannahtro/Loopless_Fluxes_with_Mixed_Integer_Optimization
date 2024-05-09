using COBREXA, Serialization
using SCIP, JuMP, Gurobi
using LinearAlgebra
using DataFrames, CSV, JSON

include("optimization_model.jl")
include("enzyme_model.jl")
include("loopless_constraints.jl")

function get_fba_data(organism="iML1515"; time_limit=1800, type="fba", save_lp=false, json=true, yeast=false, optimizer=SCIP.Optimizer, enzyme_data=false, seed=1, mean=1.2)
    # build model
    if yeast 
        molecular_model = load_model(StandardModel, "../molecular_models/ecModels/Classical/emodel_" * organism * "_classical.mat")
    else 
        molecular_model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
        print_model(molecular_model, organism)
    end

    if enzyme_data
        molecular_model = build_gecko_model(molecular_model, seed, mean)
    end 

    # filter original reactions correctly
    internal_rxn_idxs = internal_reactions(molecular_model, enzyme_data=enzyme_data)
    @show length(internal_rxn_idxs)

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    # # check for fixed reactions
    # fixed_reactions = [idx for (idx,val) in enumerate(lb) if val==ub[idx]]
    # if !isempty(fixed_reactions)
    #     @warn string(fixed_reactions) * " fixed to " * string(lb[fixed_reactions])
    # end

    # # check for fixed/ exchange (?) reactions
    # occurencs = [(i, count(==(i), S.rowval)) for i in unique(S.rowval)]
    # for (row, sum) in occurencs
    #     if sum == 1
    #         @warn string(row) * " is fixed to zero"
    #     end
    # end

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)

    if save_lp
        write_to_file(model, "../experiments/lp_models/fba_model_" * organism * ".lp")
    end

    # FBA
    objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)
    @show objective_fba
    # test feasibility, filter non-zero internal fluxes, set binaries accordingly
    solution = vars_fba[1:num_reactions]
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible, _, _ = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
    @show feasible

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = objective_fba
    dict[:dual_bound] = dual_bound
    dict[:solution] = vars_fba
    dict[:time] = time_fba
    dict[:termination] = termination_fba
    dict[:time_limit] = time_limit
    dict[:thermo_feasible] = feasible

    if json 
        if optimizer != SCIP.Optimizer
            type = type * "_" * replace(string(optimizer), ".Optimizer"=>"")
        end
        if enzyme_data
            type = type * "_gecko_" * string(seed)
            if mean != 1 
                type = type * "_" * string(mean)
            end
        end

        file_name = "json/" * organism * "_" * type * ".json"

        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
    
    return objective_fba, termination_fba, feasible
end
