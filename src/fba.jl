using COBREXA, Serialization
using SCIP, JuMP, Gurobi
using LinearAlgebra
using DataFrames, CSV, JSON

include("optimization_model.jl")
include("loopless_constraints.jl")

function get_fba_data(organism="iML1515"; time_limit=1800, type="fba", save_lp=false, json=true, yeast=false, optimizer=SCIP.Optimizer)
    # build model
    if yeast 
        molecular_model = load_model("../molecular_models/ecModels/Classical/emodel_" * organism * "_classical.mat")
    else 
        molecular_model = load_model("../molecular_models/" * organism * ".json")
        print_model(molecular_model, organism)
    end
    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)
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
        write_to_file(model, "../lp_models/fba_model_" * organism * ".lp")
    end

    # FBA
    objective_fba, dual_bound, vars_fba, time_fba, termination_fba = optimize_model(model, print_objective=true)

    # test feasibility, filter non-zero internal fluxes, set binaries accordingly
    solution = vars_fba[1:num_reactions]
    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
    non_zero_flux_directions = [solution[idx] >= 1e-5 ? 1 : 0 for idx in non_zero_flux_indices]
    feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S)
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
        file_name = "json/" * organism * "_" * type * ".json"

        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
    # dict  = JSON.parse(open("../json/" * organism * "_" * type * ".json"))  

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

    # file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    # if !isfile(file_name)
    #     CSV.write(file_name, df, append=true, writeheader=true)
    # else 
    #     CSV.write(file_name, df, append=true)    
    # end
end
