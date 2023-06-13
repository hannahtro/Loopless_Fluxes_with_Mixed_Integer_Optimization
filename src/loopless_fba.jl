using COBREXA, Serialization, COBREXA.Everything
using DataFrames, CSV
using SCIP, JuMP

include("loopless_constraints.jl")
include("optimization_model.jl")
include("cycle_detection.jl")

"""
compute dual gap with time limit of loopless FBA
"""
function loopless_fba_data(organism; time_limit=1800, silent=true, nullspace_formulation=false, type = "loopless_fba", csv=true)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    open("../csv/model_cobrexa_" * organism * ".lp", "w") do f
        print(f, model)
    end

    add_loopless_constraints(molecular_model, model, nullspace_formulation=nullspace_formulation)

    # @show model
    # open("../csv/model_vector_" * organism * ".lp", "w") do f
    #     print(f, model)
    # end

    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=false, silent=silent)

    nodes = MOI.get(model, MOI.NodeCount())
    
    # @show nodes
    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes,
        time_limit=time_limit, 
        nullspace_formulation=nullspace_formulation)

    if !nullspace_formulation
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_nullspace_" * string(time_limit) * ".csv")
    end

    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
    return objective_loopless_fba, vars_loopless_fba, time_loopless_fba, nodes
end

"""
compute dual gap with time limit of loopless FBA with blocked cycles
"""
function loopless_fba_blocked_data(organism; time_limit=180, ceiling=1000, same_objective=true, vector_formulation=true, shortest_cycles=false, block_limit=100, type="loopless_fba_blocked", nullspace_formulation=false)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")
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
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    add_loopless_constraints(molecular_model, model, nullspace_formulation=nullspace_formulation)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    num_blocked_cycles = block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S, vector_formulation=vector_formulation, shortest_cycles=shortest_cycles, block_limit=block_limit)

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
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    end

    CSV.write(file_name, df, append=false, writeheader=true)
end

"""
compute dual gap with time limit of loopless FBA with indicators
"""
function loopless_indicator_fba_data(organism; time_limit=1800, type = "loopless_indicator_fba", nullspace_formulation=true)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    # loopless FBA

    if nullspace_formulation    
        add_loopless_indicator_constraints(molecular_model, model)
    else 
        add_loopless_indicator_constraints_mu(molecular_model, model)
    end

    # @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, dual_objective_value, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, "loopless FBA", time_limit=time_limit)
    
    nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_objective_value,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes,
        time_limit=time_limit, 
        nullspace_formulation=nullspace_formulation)

    if nullspace_formulation
        type = type * "_nullspace"
    end

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

"""
compute dual gap with time limit of loopless FBA with indicators with bocked cycles
"""
function loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10, same_objective=true, shortest_cycles=false, block_limit=500, type = "loopless_indicator_fba_blocked", nullspace_formulation=false)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")
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
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    if nullspace_formulation
        add_loopless_indicator_constraints(molecular_model, model)
    else 
        add_loopless_indicator_constraints_mu(molecular_model, model)
    end
    
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    num_blocked_cycles = block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, S, shortest_cycles=shortest_cycles, block_limit=block_limit)

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
    
    if !same_objective
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    end

    CSV.write(file_name, df, append=false, writeheader=true)
end
