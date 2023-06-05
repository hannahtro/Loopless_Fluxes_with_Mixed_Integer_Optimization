using COBREXA, Serialization, COBREXA.Everything
using Serialization
using DataFrames, CSV
using SCIP, JuMP

include("../src/loopless_constraints.jl")
include("../src/optimization_model.jl")
include("../src/cycle_detection.jl")

function set_primal(organism, solution; time_limit=180)
    # load model
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, organism)

    # # FBA
    # optimizer = SCIP.Optimizer
    # model = make_optimization_model(molecular_model, optimizer)
    # # @show model
    # set_attribute(model, MOI.Silent(), true)
    # if !same_objective
    #     x = model[:x]
    #     @objective(model, Max, sum(x))
    # end
    # _, _, solution, _, _ = optimize_model(model, print_objective=false)

    # # check that FBA solution is a feasible loopless FBA solution
    # S = stoichiometry(molecular_model)
    # lb, ub = bounds(molecular_model)
    # S_transform, lb_transform, ub_transform, reaction_mapping, solution_transform  = split_hyperarcs(S, lb, ub, solution)
    # # find cycles, get original reactions
    # cycles, edge_mapping, _ = ubounded_cycles(S_transform, solution_transform, ceiling=ceiling)
    # # @show cycles
    # unbounded_cycles, unbounded_cycles_original, flux_directions = unbounded_cycles_S(cycles, edge_mapping, solution_transform, reaction_mapping)
    # # feasibility check 
    # internal_rxn_idxs = internal_reactions(molecular_model)
    # num_blocked_cycles = thermo_feasible(internal_rxn_idxs, solution, S)
    # @assert num_blocked_cycles == 0 #TODO: just use thermo_feasible function

    # get values for G and a for solution 
    S = stoichiometry(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)
    G = determine_G(S, solution, internal_rxn_idxs)   

    # loopless FBA 
    # add loopless constraints and block cycles
    # molecular_model = deserialize("../data/" * organism * ".js")
    # # print_model(molecular_model, organism)

    # model = make_optimization_model(molecular_model, optimizer)
    # add_loopless_constraints(molecular_model, model)
    
    # add primal solution

    # @show model
    # set_attribute(model, MOI.Silent(), false)
    # objective_loopless_fba, dual_bound, vars_loopless_fba, time_loopless_fba, termination_loopless_fba = 
    #     optimize_model(model, type, time_limit=time_limit, print_objective=false)

    # nodes = MOI.get(model, MOI.NodeCount())

    # df = DataFrame(
    #     objective_value=objective_loopless_fba, 
    #     dual_bound=dual_bound,
    #     solution=vars_loopless_fba, 
    #     time=time_loopless_fba, 
    #     termination=termination_loopless_fba,
    #     nodes=nodes,
    #     num_blocked_cycles=num_blocked_cycles)

    # type = "set_primal"
    # if !same_objective
    #     file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * ".csv")
    # else 
    #     file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * "_" * string(ceiling) * "_same_objective.csv")
    # end

    # CSV.write(file_name, df, append=false, writeheader=true)
end

# read from csv
# TDOO: load FBA solution
file_name = "../csv/iJR904/iJR904_loopless_fba_300.csv"
df = DataFrame(CSV.File(file_name))
solution = df[!,:solution]

organism = "iJR904"
set_primal(organism, solution, time_limit=1800)
