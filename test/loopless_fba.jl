using COBREXA, Serialization, COBREXA.Everything
using Serialization
using DataFrames, CSV
using SCIP, JuMP

include("../src/loopless_constraints.jl")
include("../src/optimization_model.jl")
include("../src/cycle_detection.jl")

# compute dual gap with time limit of loopless FBA
function loopless_fba_data(organism; time_limit=1800)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)

    type = "loopless_fba"
    add_loopless_constraints(molecular_model, model)
    # @show model
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=true)

    nodes = MOI.get(model, MOI.NodeCount())
    
    @show nodes
    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

# compute dual gap with time limit of loopless FBA with blocked cycles
function loopless_fba_blocked_data(organism; time_limit=180, ceiling=10, same_objective=true, vector_formulation=true)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    # compute FBA
    optimizer = SCIP.Optimizer
    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    type = "fba"
    if !same_objective
        x = model[:x]
        @objective(model, Max, sum(x))
    end
    _, _, solution, _, _ = optimize_model(model, print_objective=true)

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
    add_loopless_constraints(molecular_model, model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs, vector_formulation=vector_formulation)

    # optimize loopless FBA
    if vector_formulation
        type = "loopless_fba_blocked"
    else 
        type = "loopless_fba_blocked_for_loop"
    end

    # @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=true)

    nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes)

    if !same_objective
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    else 
        file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    end

    CSV.write(file_name, df, append=false, writeheader=true)
end

# compute dual gap with time limit of loopless FBA with indicators
function loopless_indicator_fba_data(organism; time_limit=1800)
    # build model
    optimizer = SCIP.Optimizer
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    # loopless FBA

    type = "loopless_indicator_fba"
    add_loopless_indicator_constraints(molecular_model, model)
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
        nodes=nodes)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

# compute dual gap with time limit of loopless FBA with indicators with bocked cycles
function loopless_indicator_fba_blocked_data(organism; time_limit=1800, ceiling=10)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    # compute FBA
    optimizer = SCIP.Optimizer
    model = make_optimization_model(molecular_model, optimizer)
    # @show model
    set_attribute(model, MOI.Silent(), true)
    _, _, solution, _, _ = optimize_model(model, print_objective=true)

    # split hyperarcs
    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    objective_function = objective(molecular_model)
    @show objective_function
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
    add_loopless_indicator_constraints(molecular_model, model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    block_cycle_constraint(model, unbounded_cycles_original, flux_directions, internal_rxn_idxs)

    # optimize loopless FBA
    type = "loopless_indicator_fba_blocked"
    # @show model
    set_attribute(model, MOI.Silent(), false)
    objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
        optimize_model(model, type, time_limit=time_limit, print_objective=true)

    nodes = MOI.get(model, MOI.NodeCount())

    df = DataFrame(
        objective_value=objective_loopless_fba, 
        dual_bound=dual_bound,
        solution=vars_loopless_fba, 
        time=time_loopless_fba, 
        termination=termination_loopless_fba,
        nodes=nodes)

    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")

    CSV.write(file_name, df, append=false, writeheader=true)
end

# organism = "iJR904"
# loopless_fba_data(organism, time_limit=3600)

# loopless_indicator_fba_data(organism, time_limit=3600)

# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=500)

# organism = "iAF692"
# loopless_fba_data(organism, time_limit=3600)

# loopless_indicator_fba_data(organism, time_limit=3600)

# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=500)


# organism = "iML151"
# loopless_fba_data(organism, time_limit=3600)

# loopless_indicator_fba_data(organism, time_limit=3600)

# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=50)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=100)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=200)
# loopless_fba_blocked_data(organism, time_limit=3600, ceiling=500)

organism = "iJR904"
loopless_fba_blocked_data(organism, time_limit=10, ceiling=50, vector_formulation=false)
loopless_fba_blocked_data(organism, time_limit=10, ceiling=50, vector_formulation=true)

# loopless_fba_blocked_data(organism, time_limit=600, ceiling=50, same_objective=false, vector_formulation=false)
# loopless_fba_data(organism, time_limit=600)
