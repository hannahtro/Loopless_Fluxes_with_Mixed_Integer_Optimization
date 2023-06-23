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
    # print(master_problem)
end

# TODO: pass MIS
function build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)
    set_attribute(sub_problem, MOI.Silent(), true)
    set_objective_sense(sub_problem, MAX_SENSE)
    S_int = Array(S[:, internal_rxn_idxs])

    G = @variable(sub_problem, G[1:length(internal_rxn_idxs)])
    μ = @variable(sub_problem, μ[1:size(S)[1]])
    constraint_list = []
    for (idx,val) in enumerate(solution_a)
        if val == 0
            c = @constraint(sub_problem, G[idx] <= -1) #TODO: check loopless fba formulation
            push!(constraint_list,c)
        elseif val == 1
            c = @constraint(sub_problem, G[idx] >= 1)      
            push!(constraint_list,c)
        end
    end
    c_matrix = @constraint(sub_problem, G' .== μ' * S_int)
    # push!(constraint_list,c)
    # print(sub_problem)
    return constraint_list, c_matrix
end

# TODO: implement fast MIS Search
# TODO: block zeros and ones??
# TODO: compute several MISs at once
function compute_MIS(solution_a, S_int; fast=true)
    # build mis_search_primal
    # build mis_search_dual
    # add objective with weights
    # C = 1:length(solution_a)
    if !fast
        C = [idx for (idx,val) in enumerate(solution_a) if val==1]
    else 
        @show solution_a
        @show S_int
        @show S_int'

        A = deepcopy(S_int')
        for (idx,a) in enumerate(solution_a)
            if a == 0
                A[idx,:] = - A[idx,:] # update rows
            end
        end
        @show S_int'
        @show A
        # @show A[3,:] #TODO verify idex

        b = [1 for i in 1:size(S_int)[1]]

        @show b
        C = []

        mis_model = Model(SCIP.Optimizer)
        @variable(mis_model, λ[1:length(b)])
        @constraint(mis_model, λ .>= 0)
        @constraint(mis_model, A' * λ .== 0)
        @constraint(mis_model, b'*λ==1)
        @objective(mis_model, Min, sum(λ))

        print(mis_model)
        optimize!(mis_model)

        @show termination_status(mis_model)
        @show MOI.get(mis_model, MOI.ObjectiveValue())
        @show [value(var) for var in all_variables(mis_model)]


        # @show MOI.get(sub_problem, MOI.ObjectiveBound())
        # # @show result_count(sub_problem)
        # @show dual_objective_value(sub_problem)

        # dual_values = []
        # for c in constraint_list
        #     push!(dual_values, dual(c))
        # end
        # append!(dual_values, dual.(c_matrix))
        # @show dual_values

        # dual_sub_problem = dualize(sub_problem, HiGHS.Optimizer; dual_names = DualNames("dual", ""))
        # print(dual_sub_problem)
        # optimize!(dual_sub_problem)
        # @show MOI.get(dual_sub_problem, MOI.ObjectiveValue())

        # obj = objective_function(dual_sub_problem)
        # @constraint(dual_sub_problem, obj==1)
        # @objective(dual_sub_problem, Min, 0)

        # print(dual_sub_problem)
        # optimize!(dual_sub_problem)
        # @show MOI.get(dual_sub_problem, MOI.ObjectiveValue())
        # repeat until solution is feasible
        # @show C
    end
    return C
end

function add_combinatorial_benders_cut(master_problem, solution_a, C; fast=true)
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

function combinatorial_benders(master_problem, internal_rxn_idxs, S; max_iter=Inf, fast=true)
    @show fast
    m, num_reactions = size(S)

    # solve master problem
    build_master_problem(master_problem, internal_rxn_idxs)   
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem)
    solutions = [solution_master]
    solution_a = solution_master[num_reactions+1:end]
    
    # compute corresponding MIS
    S_int = Array(S[:, internal_rxn_idxs])
    C = compute_MIS(solution_a, fast=fast, S_int)

    # build sub problem to master solution 
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    constraint_list, c_matrix = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)

    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=false)

    # if infeasible, add Benders' cut
    iter = 1
    while termination_sub == MOI.INFEASIBLE && iter < max_iter
        @assert primal_status(sub_problem) == MOI.NO_SOLUTION
        @assert dual_status(sub_problem) == MOI.INFEASIBILITY_CERTIFICATE
        @assert has_duals(sub_problem)

        add_combinatorial_benders_cut(master_problem, solution_a, C) 
        objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem)
        @assert !(solution_master in solutions)
        push!(solutions, solution_master)
        solution_a = solution_master[num_reactions+1:end]
        # @show solution_a

        sub_problem = Model(optimizer)
        constraint_list, c_matrix = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)
        objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=false)

        iter += 1
    end
end