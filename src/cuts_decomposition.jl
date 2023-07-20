using Dates 
using HiGHS
using Dualization 

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

    start_time = time()
    dual_bounds = []

    objective_value, dual_bound, solution, _, termination = optimize_model(model)

    push!(dual_bounds, dual_bound)
    solutions = [round.(solution, digits=5)]
    cuts = []

    iter = 1
    while !thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S) && time()-start_time < time_limit
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

        cut = @constraint(model, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(internal_rxn_idxs)-1)
        @assert !(cut in cuts)
        push!(cuts,[cut])

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

    end_time = time()
    time_taken = end_time - start_time
    # @show time_taken

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

"""
build master problem of combinatorial Benders decomposition with FBA constraints and indicator variables,
maps indicator variables to flux direction
"""
function build_master_problem(master_problem, internal_rxn_idxs)
    set_attribute(master_problem, MOI.Silent(), true)
    x = master_problem[:x]

    # add indicator variables 
    a = @variable(master_problem, a[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(master_problem, a[cidx] => {x[ridx] - eps() >= 0})
        @constraint(master_problem, !a[cidx] => {x[ridx] + eps() <= 0})
    end
end

"""
build sub problem of combinatorial Benders decomposition including the thermodynamic constraints on the indicator variables 
for a given solution to the master problem and the minimal infeasible subset C
"""
function build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)
    # @show solution_a
    # @show C
    set_attribute(sub_problem, MOI.Silent(), true)
    set_objective_sense(sub_problem, MAX_SENSE)
    S_int = Array(S[:, internal_rxn_idxs])

    # G = @variable(sub_problem, G[1:length(internal_rxn_idxs)])
    μ = @variable(sub_problem, μ[1:size(S_int)[1]])
    constraint_list = []
    for (idx,val) in enumerate(solution_a) 
        if val == 0 && (idx in C)
            c = @constraint(sub_problem, (S_int' * μ)[idx] <= -1) #TODO: check loopless fba formulation
            push!(constraint_list,c)
        elseif val == 1 && (idx in C)
            c = @constraint(sub_problem, (S_int' * μ)[idx] >= 1)      
            push!(constraint_list,c)
        else
            @assert (idx in C) == false
        end
    end
    # c_matrix = @constraint(sub_problem, G .== S_int' * μ)
    return constraint_list #, c_matrix
end

"""
returns a minimal infeasible subset of reactions for a give soltion and stoichiometric matrix
"""
# TODO: compute several MISs at once
function compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs; fast=true, time_limit=1800, silent=true)
    if !fast
        # not a MIS
        # C = [idx for (idx,val) in enumerate(solution_a) if val==1]
        C = [idx for (idx,val) in enumerate(solution_a)]
        # C = [1,2,3]
    else 
        A = deepcopy(S_int)
        for (idx,a) in enumerate(solution_a)
            if a == 0
                A'[idx,:] = - A'[idx,:] # update rows
            end
        end
        b = [1 for i in 1:length(internal_rxn_idxs)]

        model = Model(HiGHS.Optimizer)
        μ = @variable(model, μ[1:size(S_int)[1]])
        @constraint(model, A' * μ .>= b)
        @objective(model, Max, 0)

        optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
        mis_model = dualize(model, optimizer, dual_names = DualNames("dual","dual_constraints"))
        set_attribute(mis_model, MOI.Silent(), true)

        @constraint(mis_model, b' * all_variables(mis_model) == 1)
        @objective(mis_model, Min, sum(all_variables(mis_model)))

        # print(mis_model)
        if silent
            set_attribute(mis_model, MOI.Silent(), true)
        else 
            set_attribute(mis_model, MOI.Silent(), false)
        end
        set_time_limit_sec(mis_model, time_limit)
        optimize!(mis_model)

        if termination_status(mis_model) != MOI.OPTIMAL
            println("MIS problem not feasible")
            C = []
        else
            solution_mis = [value(var) for var in all_variables(mis_model)]
            C = [idx for (idx,val) in enumerate(solution_mis) if !(isapprox(val,0))]
        end

        # # λ'Aμ ≥ λ'b should be violated
        # # λ solution to sub problem, A constructed in fast MIS search, 
        # # μ flux values of master problem, b constructed in fast MIS search
        # @assert !(solution_mis' * A * solution_master[internal_rxn_idxs] >= solution_mis' * b)
    end
    return C
end

"""
adds combinatorial Benders' cut to the master problem, by forcing a different assignment of the indicator variables
of the reactions in the minimal infeasible subset C
"""
function add_combinatorial_benders_cut(master_problem, solution_a, C)
    a = master_problem[:a]
    Z = []
    O = []
    for idx in C
        if solution_a[idx] > 0 
            push!(O,idx)
        else 
            push!(Z,idx)
        end
    end 
    # @show Z,O
    if isempty(Z)
        @constraint(master_problem, sum(a[O]) <= length(C)-1)
    elseif isempty(O)
        @constraint(master_problem, sum([1-a[i] for i in Z]) <= length(C)-1)
    else 
        @constraint(master_problem, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(C)-1)
    end
end 

"""
solve problem by splitting it into a master problem with indicator variables and a linear sub problem based 
on a solution to the master problem and minimal infeasible subsets. The sub problem 
"""
function combinatorial_benders(master_problem, internal_rxn_idxs, S; max_iter=Inf, fast=true, time_limit=1800, silent=true, save_model=false, organism="")
    @show fast
    _, num_reactions = size(S)

    start_time = time()
    dual_bounds = []
    objective_values = []

    # solve master problem
    build_master_problem(master_problem, internal_rxn_idxs)   
    if save_model
        open("../csv/model_" * organism * ".lp", "w") do f
            print(f, master_problem)
        end
    end 
    @show objective_function(master_problem)
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem)
    solution_master = round.(solution_master, digits=6)
    solutions = [solution_master]
    solution_a = solution_master[num_reactions+1:end]
    push!(dual_bounds, dual_bound_master)
    push!(objective_values, objective_value_master)

    # compute corresponding MIS
    S_int = Array(S[:, internal_rxn_idxs])
    # @show size(S_int)
    C = compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs, fast=fast, time_limit=time_limit)

    # build sub problem to master solution 
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)

    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
    # @show C

    # add Benders' cut if subproblem is infeasible
    iter = 1
    while termination_sub == MOI.INFEASIBLE && iter <= max_iter && time()-start_time < time_limit
        @show iter
        @assert primal_status(sub_problem) == MOI.NO_SOLUTION
        @assert dual_status(sub_problem) == MOI.INFEASIBILITY_CERTIFICATE
        @assert has_duals(sub_problem)

        # add CB cut to MP and solve MP
        add_combinatorial_benders_cut(master_problem, solution_a, C)
        if save_model
            open("../csv/model_" * organism * "_" * string(iter) * ".lp", "w") do f
                print(f, master_problem)
            end
        end 
        println("_______________")
        println("master problem")
        objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
        solution_master = round.(solution_master, digits=5)
        @assert !(solution_master in solutions)
        @assert termination_master == MOI.OPTIMAL
        push!(solutions, solution_master)
        push!(dual_bounds, dual_bound_master)
        push!(objective_values, objective_value_master)
        solution_a = solution_master[num_reactions+1:end]
        # @show solution_a

        # compute corresponding MIS
        println("_______________")
        println("compute MIS")
        C = compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs, fast=fast, time_limit=time_limit, silent=silent)
        if isempty(C)
            feasible = thermo_feasible_mu(internal_rxn_idxs, solution_master[internal_rxn_idxs], S)
            @assert feasible
            sub_problem = Model(optimizer)
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, internal_rxn_idxs)
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
            @assert termination_sub == MOI.OPTIMAL
        else 
            # build sub problem to master solution 
            sub_problem = Model(optimizer)
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)
            println("_______________")
            println("sub problem")
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
            iter += 1
        end
    end
    end_time = time()
    time_taken = end_time - start_time
    solution = vcat(solution_master, solution_sub)
    # @show termination_sub

    if termination_sub == MOI.OPTIMAL
        feasible = thermo_feasible_mu(internal_rxn_idxs, solution[internal_rxn_idxs], S)
        @assert feasible
    end 

    return objective_value_master, objective_values, dual_bounds, solution, time_taken, termination_sub, iter
end

function combinatorial_benders_data(organism; time_limit=1800, csv=true, max_iter=Inf, fast=true, silent=true, save_model=false)
    @show fast
    molecular_model = deserialize("../data/" * organism * ".js")
    print_model(molecular_model, "organism")

    S = stoichiometry(molecular_model)
    # lb, ub = bounds(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    # model = build_fba_model(S, lb, ub, optimizer=SCIP.Optimizer)
    master_problem = make_optimization_model(molecular_model, SCIP.Optimizer)
   
    objective_value, objective_values, dual_bounds, solution, time, termination, iter = combinatorial_benders(master_problem, internal_rxn_idxs, S, max_iter=max_iter, fast=fast, silent=silent, save_model=save_model, organism=organism)

    @show termination
    @show objective_value
    df = DataFrame(
        objective_value=objective_value, 
        dual_bounds=[dual_bounds],
        objective_values=[objective_values],
        solution=[solution], 
        time=time, 
        termination=termination,
        time_limit=time_limit, 
        iter=iter)

    type = "combinatorial_benders"
    if fast
        type = type * "_fast"
    end
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end