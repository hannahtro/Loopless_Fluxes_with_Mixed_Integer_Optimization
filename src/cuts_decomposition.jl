using Dates 
using HiGHS

include("optimization_model.jl")
include("loopless_constraints.jl")

"""
add an additional cycle to block until the solution is thermodynamically feasible
"""
function no_good_cuts(model, internal_rxn_idxs, S; time_limit=1800)
    x = model[:x]
    m, num_reactions = size(S)

    # add indicator variables 
    a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(model, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(model, !a[cidx] => {x[ridx] + eps() <= 0})
    end
    # @objective(model, Max, 0)
    # println("here")

    start_time = time()
    # println("here")
    dual_bounds = []

    objective_value, dual_bound, solution, _, termination = optimize_model(model)
    # @show solution
    push!(dual_bounds, dual_bound)
    solutions = [round.(solution, digits=5)]
    cuts = []
    # print(model)
    iter = 1
    while !thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S) && time()-start_time < time_limit
        # @show thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
        solution_a = solution[num_reactions+1:end]
        @assert round.(solution_a) == solution_a

        Z = []
        O = []
        for (idx, ridx) in enumerate(internal_rxn_idxs)
            if solution_a[idx] > 0 
                push!(O,idx)
            else 
                push!(Z,idx)
            end
        end 
        # @show Z,O

        cut = @constraint(model, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(internal_rxn_idxs)-1)
        @assert !(cut in cuts)
        push!(cuts,[cut])
        # print(model)

        objective_value, dual_bound, solution, _, termination = optimize_model(model)
        solution = round.(solution, digits=5)
        solution_a = solution[num_reactions+1:end]

        @assert solutions[end][num_reactions+1:end] != solution_a
        @assert sum(solution_a[O]) + sum([1-solution_a[i] for i in Z]) <= length(internal_rxn_idxs)-1
        @assert !(solution in solutions)
        push!(solutions,solution)
        push!(dual_bounds, dual_bound)
        iter += 1
    end
    # @show cuts
    # @show solutions
    end_time = time()
    time_taken = end_time - start_time
    @show time_taken

    return objective_value, dual_bounds, solution, time_taken, termination, iter
end

function no_good_cuts_data(organism; time_limit=1800, csv=true)
    model = deserialize("../data/" * organism * ".js")
    print_model(model, "organism")

    S = stoichiometry(model)
    lb, ub = bounds(model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(model)) if
        !is_boundary(reaction_stoichiometry(model, rid))
    ]

    model = build_fba_model(S, lb, ub)
    objective_value, dual_bounds, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

    thermo_feasible = thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)

    df = DataFrame(
        objective_value=objective_value, 
        dual_bounds=[dual_bounds],
        solution=[solution], 
        time=time, 
        termination=termination,
        time_limit=time_limit, 
        thermo_feasible=thermo_feasible,
        iter=iter)

    type = "no_good_cuts"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end

function combinatorial_benders(master_problem, internal_rxn_idxs, S)
    x = master_problem[:x]
    m, num_reactions = size(S)
    S_int = Array(S[:, internal_rxn_idxs])

    # add indicator variables 
    a = @variable(master_problem, a[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(master_problem, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(master_problem, !a[cidx] => {x[ridx] + eps() <= 0})
    end

    objective_value, dual_bound, solution, _, termination = optimize_model(master_problem)
    solution_a = solution[num_reactions+1:end]

    # build sub problem to solution 
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)

    G = @variable(sub_problem, G[1:length(internal_rxn_idxs)])
    μ = @variable(sub_problem, μ[1:size(S)[1]])
    for (idx,val) in enumerate(solution_a)
        if val == 0
            @constraint(sub_problem, G[idx] <= -1) #TODO: check loopless fba formulation
        elseif val == 1
            @constraint(sub_problem, G[idx] >= 1)      
        end
    end
    @constraint(sub_problem, G' .== μ' * S_int)

    objective_value, dual_bound, solution, _, termination = optimize_model(sub_problem, silent=false)

    # if infeasible, block ray from dual program in master sub_problem
    @show termination
    @show primal_status(sub_problem)
    @show dual_status(sub_problem)
    @show has_duals(sub_problem)

    @show MOI.get(sub_problem, MOI.ObjectiveBound())
    @show result_count(sub_problem)
    @show dual_objective_value(sub_problem)
    # repeat until solution is feasible
end