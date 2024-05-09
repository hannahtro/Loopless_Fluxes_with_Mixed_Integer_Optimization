using COBREXA, Serialization
using DataFrames, CSV, JSON
using SCIP, JuMP, GLPK, HiGHS

include("loopless_constraints.jl")
include("optimization_model.jl")

"""
compute dual gap with time limit of loopless FBA
"""
function loopless_fba_data(organism; time_limit=1800, silent=true, nullspace_formulation=false, type="loopless_fba", json=true, yeast=false, optimizer=SCIP.Optimizer, max_reactions=[], print_objective=false, scip_tol=1.0e-6)
    # build model
    if yeast 
        molecular_model = load_model("../molecular_models/" * organism * ".xml")
    else 
        molecular_model = load_model("../molecular_models/" * organism * ".json")
        print_model(molecular_model, organism)
    end

    @show max_reactions
    if !isempty(max_reactions)
        molecular_model.objective = Dict{String, Float64}()
        for molecule in max_reactions
            molecular_model.objective[molecule] = 1.0
        end
    end

    model = make_optimization_model(molecular_model, optimizer)
    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    # model = build_fba_model(S, lb, ub, max_reactions=max_reactions)
    internal_rxn_idxs = internal_reactions(molecular_model)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))

    m, num_reactions = size(S)

    reaction_mapping = Dict()
    for (idx, val) in enumerate(internal_rxn_idxs)
        reaction_mapping[val] = idx
    end
    # xl, xu = bounds(molecular_model)
    # model = build_fba_model(S, xl, xu)
    # x = model[:x]
    # @objective(model, MAX_SENSE, objective(molecular_model)' * x)

    # @show model
    # open("../experiments/csv/model_cobrexa_" * organism * ".lp", "w") do f
    #     print(f, model)
    # end

    add_loopless_constraints(molecular_model, model, nullspace_formulation=nullspace_formulation, max_flux_bound)

    # @show model
    # open("../experiments/csv/model_vector_" * organism * ".lp", "w") do f
    #     print(f, model)
    # end

    objective_loopless_fba, dual_bound, _, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, silent=silent, print_objective=print_objective)

    if optimizer == SCIP.Optimizer || optimizer == Gurobi.Optimizer
        nodes = MOI.get(model, MOI.NodeCount())
    else 
        nodes = NaN
    end
    if termination_loopless_fba == MOI.OPTIMAL
        S = stoichiometry(molecular_model)
        x = value.(model[:x])
        a = value.(model[:a])
        G = value.(model[:G])
        if !nullspace_formulation
            μ = value.(model[:μ])
        else 
            μ = NaN
        end 
        
        steady_state =  isapprox.(S * x, 0, atol=0.0001)
        @assert steady_state == ones(size(S)[1])
        # test feasibility, filter non-zero fluxes, set binaries accordingly
	    non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(x) if !isapprox(val, 0, atol=1e-4)], internal_rxn_idxs)
        non_zero_flux_directions = a[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
        @assert feasible
    else 
        feasible = false
        x = NaN
        a = NaN 
        G = NaN 
        μ = NaN
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = objective_loopless_fba
    dict[:dual_bound] = dual_bound
    dict[:x] = x
    dict[:a] = a
    dict[:G] = G
    dict[:μ] = μ
    dict[:time] = time_loopless_fba
    dict[:termination] = termination_loopless_fba
    dict[:nodes] = nodes
    dict[:time_limit] = time_limit
    dict[:nullspace_formulation] = nullspace_formulation
    dict[:thermo_feasible] = feasible
    dict[:max_flux_bound] = max_flux_bound
    dict[:objective_function] = objective(molecular_model)
    dict[:sense] = objective_sense(model)

    if nullspace_formulation
        type = type * "_nullspace"
    end
    if !isempty(max_reactions)
        type = type * "_modified_objective"
    end
    if optimizer != SCIP.Optimizer
        type = type * "_" * replace(string(optimizer), ".Optimizer"=>"")
    end
    
    file_name = joinpath("json/" * organism * "_" * type * "_" * string(time_limit) * ".json")
    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
    return objective_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes
end

"""
compute dual gap with time limit of loopless FBA
"""
function loopless_relaxed_fba_data(organism; time_limit=1800, silent=true, nullspace_formulation=false, type = "loopless_fba_relaxed", csv=true, save_lp=false)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    S = stoichiometry(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)

    add_relaxed_loopless_constraints(model, S, internal_rxn_idxs, nullspace_formulation=nullspace_formulation)

    if nullspace_formulation
        type = type * "_nullspace"
    end

    if save_lp 
        write_to_file(model, "../experiments/csv/models/" * type * "_" * organism * ".lp")
    end 

    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false, silent=silent)

    nodes = MOI.get(model, MOI.NodeCount())
    if termination_loopless_fba == MOI.OPTIMAL
        S = stoichiometry(molecular_model)
        steady_state =  isapprox.(S * vars_loopless_fba[1:size(S)[2]],0, atol=0.0001)
        @assert steady_state == ones(size(S)[1])
    end

    # @show nodes
    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=[vars_loopless_fba], 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes,
        time_limit=time_limit, 
        nullspace_formulation=nullspace_formulation)

    if csv
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
        CSV.write(file_name, df, append=false, writeheader=true)
    end
    return objective_loopless_fba, vars_loopless_fba, time_loopless_fba, nodes
end

"""
compute dual gap with time limit of loopless FBA with indicator for bilinear constraints
"""
function loopless_fba_bilinear_data(organism; time_limit=1800, silent=true, type="loopless_bilinear_fba", csv=true)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # S = stoichiometry(molecular_model)
    # xl, xu = bounds(molecular_model)
    # model = build_fba_model(S, xl, xu)
    # x = model[:x]
    # @objective(model, MAX_SENSE, objective(molecular_model)' * x)

    # @show model
    # open("../experiments/csv/model_cobrexa_" * organism * ".lp", "w") do f
    #     print(f, model)
    # end

    add_loopless_bilinear_constraints(molecular_model, model)

    # @show model
    # open("../experiments/csv/model_vector_" * organism * ".lp", "w") do f
    #     print(f, model)
    # end

    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false, silent=silent)

    nodes = MOI.get(model, MOI.NodeCount())
    if termination_loopless_fba == MOI.OPTIMAL
        S = stoichiometry(molecular_model)
        steady_state =  isapprox.(S * vars_loopless_fba[1:size(S)[2]],0, atol=0.0001)
        @assert steady_state == ones(size(S)[1])
    end

    # @show nodes
    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=[vars_loopless_fba], 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes,
        time_limit=time_limit
    )
    
    file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
    return objective_loopless_fba, vars_loopless_fba, time_loopless_fba, nodes
end

"""
compute dual gap with time limit of loopless FBA with blocked cycles
"""
function loopless_fba_blocked_data(organism; time_limit=180, ceiling=1000, same_objective=true, vector_formulation=true, shortest_cycles=false, block_limit=100, type="loopless_fba_blocked", nullspace_formulation=false, reduced=false, csv=true)
    # load model
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    # compute FBA
    optimizer = SCIP.Optimizer
    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    if !same_objective
        x = model[:x]
        @objective(model, Max, sum(x))
    end
    _, _, solution, _, _ = optimize_model(model, print_objective=false)

    # split hyperarcs
    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform  = split_hyperarcs(S, lb, ub, solution)

    # find cycles, get original reactions
    cycles, edge_mapping, _ = ubounded_cycles(S_transform, solution_transform, ceiling=ceiling)
    # @show cycles
    unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)

    # build model
    # add loopless constraints and block cycles
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    add_loopless_constraints(molecular_model, model, nullspace_formulation=nullspace_formulation, reduced=reduced)
    internal_rxn_idxs = internal_reactions(molecular_model)
    num_blocked_cycles = block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S, vector_formulation=vector_formulation, shortest_cycles=shortest_cycles, block_limit=block_limit, nullspace_formulation=nullspace_formulation)

    # optimize loopless FBA
    if !vector_formulation
        type = type * "_for_loop" 
    end
    if shortest_cycles
        type = type * "_shortest_cycles"
    end
    if nullspace_formulation
        type = type * "_nullspace"
    end
    if reduced
        type = type * "_reduced"
    end
    
    type = type * "_" * string(block_limit)

    # @show model
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
        nodes=nodes,
        num_blocked_cycles=num_blocked_cycles,
        ceiling=ceiling,
        time_limit=time_limit, 
        shortest_cycles=shortest_cycles, 
        nullspace_formulation=nullspace_formulation, 
        block_limit=block_limit)

    if !same_objective
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    end

    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end

"""
compute dual gap with time limit of loopless FBA with indicators
"""
function loopless_indicator_fba_data(organism; time_limit=1800, type="loopless_indicator_fba", nullspace_formulation=true, json=true, silent=true, print_objective=true, yeast=false)
    # build model
    optimizer = SCIP.Optimizer

    if yeast 
        molecular_model = load_model("../molecular_models/" * organism * ".xml")
    else 
        molecular_model = load_model("../molecular_models/" * organism * ".json")
    end
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    internal_rxn_idxs = internal_reactions(molecular_model)

    reaction_mapping = Dict()
    for (idx, val) in enumerate(internal_rxn_idxs)
        reaction_mapping[val] = idx
    end

    # @show model
    # loopless FBA
    if nullspace_formulation    
        add_loopless_indicator_constraints(molecular_model, model)
    else 
        add_loopless_indicator_constraints_mu(molecular_model, model)
    end

    # @show model
    objective_loopless_fba, dual_bound, _, time_loopless_fba, termination_loopless_fba = 
    optimize_model(model, type, time_limit=time_limit, silent=silent, print_objective=print_objective)
    @show objective_loopless_fba, termination_loopless_fba

    if optimizer == SCIP.Optimizer || optimizer == Gurobi.Optimizer
        nodes = MOI.get(model, MOI.NodeCount())
    else 
        nodes = NaN
    end
    if termination_loopless_fba == MOI.OPTIMAL
        S = stoichiometry(molecular_model)
        x = value.(model[:x])
        a = value.(model[:a])
        G = value.(model[:G])
        if !nullspace_formulation
            μ = value.(model[:μ])
        else 
            μ = NaN
        end 
        
        steady_state =  isapprox.(S * x, 0, atol=0.0001)
        @assert steady_state == ones(size(S)[1])
        # test feasibility, filter non-zero fluxes, set binaries accordingly
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(x) if !isapprox(val, 0, atol=1e-4)], internal_rxn_idxs)
        non_zero_flux_directions = a[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
        @assert feasible
    else 
        feasible = false
        x = NaN
        a = NaN 
        G = NaN 
        μ = NaN
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = objective_loopless_fba
    dict[:dual_bound] = dual_bound
    dict[:x] = x
    dict[:a] = a
    dict[:G] = G
    dict[:μ] = μ
    dict[:time] = time_loopless_fba
    dict[:termination] = termination_loopless_fba
    dict[:nodes] = nodes
    dict[:time_limit] = time_limit
    dict[:nullspace_formulation] = nullspace_formulation
    dict[:thermo_feasible] = feasible     
    dict[:objective_function] = objective(molecular_model)
    dict[:sense] = objective_sense(model)

    if nullspace_formulation
        type = type * "_nullspace"
    end
    if optimizer != SCIP.Optimizer
        type = type * "_" * replace(string(optimizer), ".Optimizer"=>"")
    end

    file_name = joinpath("json/" * organism * "_" * type * "_" * string(time_limit) * ".json")
    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
end

"""
compute dual gap with time limit of loopless FBA with indicators with bocked cycles
"""
function loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10, same_objective=true, shortest_cycles=false, block_limit=500, type = "loopless_indicator_fba_blocked", nullspace_formulation=false, csv=true)
    # load model
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    # compute FBA
    optimizer = SCIP.Optimizer
    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    if !same_objective
        x = model[:x]
        @objective(model, Max, sum(x))
    end
    _, _, solution, _, _ = optimize_model(model, print_objective=false)

    # split hyperarcs
    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    objective_function = objective(molecular_model)
    # @show objective_function
    S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform  = split_hyperarcs(S, lb, ub, solution)

    # find cycles, get original reactions
    cycles, edge_mapping, _ = ubounded_cycles(S_transform, solution_transform, ceiling=ceiling)
    # @show cycles
    unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)

    # build model
    # add loopless constraints and block cycles
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    if nullspace_formulation
        add_loopless_indicator_constraints(molecular_model, model)
    else 
        add_loopless_indicator_constraints_mu(molecular_model, model)
    end
    
    internal_rxn_idxs = internal_reactions(molecular_model)
    num_blocked_cycles = block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S, shortest_cycles=shortest_cycles, block_limit=block_limit, nullspace_formulation=nullspace_formulation)

    # optimize loopless FBA
    # @show model
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
        nodes=nodes,
        num_blocked_cycles=num_blocked_cycles,
        ceiling=ceiling,
        time_limit=time_limit, 
        shortest_cycles=shortest_cycles, 
        nullspace_formulation=nullspace_formulation, 
        block_limit=block_limit)

    if nullspace_formulation
        type = type * "_nullspace"
    end
    type = type * "_" * string(block_limit)
    
    if !same_objective
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    end

    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end
