using Dates 
using HiGHS
using Dualization 
using JSON

include("utils.jl")
include("optimization_model.jl")
include("loopless_constraints.jl")

# """
# add an additional cycle to block until the solution is thermodynamically feasible
# """
# function no_good_cuts(model, internal_rxn_idxs, S; time_limit=1800)
#     x = model[:x]
#     m, num_reactions = size(S)

#     # TODO: replace by build master problem
#     # add indicator variables 
#     build_master_problem(model, internal_rxn_idxs)
#     # a = @variable(model, a[1:length(internal_rxn_idxs)], Bin)
#     # for (cidx, ridx) in enumerate(internal_rxn_idxs)
#     #     # add indicator 
#     #     @constraint(model, a[cidx] => {x[ridx] >= -0.000001})
#     #     @constraint(model, !a[cidx] => {x[ridx] <= 0.000001})
#     # end
#     # @objective(model, Max, 0)

#     start_time = time()
#     dual_bounds = []

#     objective_value, dual_bound, solution, _, termination = optimize_model(model)

#     if termination != MOI.OPTIMAL
#         end_time = time()
#         time_taken = end_time - start_time
#         return objective_value, dual_bounds, solution, time_taken, termination, 0
#     end

#     solution_a = solution[num_reactions+1:end]
#     push!(dual_bounds, dual_bound)
#     solutions = [round.(solution, digits=5)]
#     cuts = []

#     iter = 1
#     while !thermo_feasible_mu(internal_rxn_idxs, solution_a, S) && time()-start_time < time_limit
#         @show iter
#         @assert isapprox(round.(solution_a), solution_a, atol=1e-4)

#         C = [i for i in 1:length(internal_rxn_idxs)]
#         add_combinatorial_benders_cut(model, solution_a, C, cuts)
#         # Z = []
#         # O = []
#         # for (idx, ridx) in enumerate(internal_rxn_idxs)
#         #     if solution_a[idx] > 1e-4 
#         #         push!(O, idx)
#         #     else 
#         #         push!(Z, idx)
#         #     end
#         # end 

#         # cut = @constraint(model, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(internal_rxn_idxs) - 1)
#         # @assert !(cut in cuts)
#         # push!(cuts,[cut])

#         objective_value, dual_bound, solution, _, termination = optimize_model(model, time_limit=time_limit)
        
#         if termination != MOI.OPTIMAL
#             if termination == MOI.TIME_LIMIT
#                 @warn "master problem cannot be solved"
#             end
#             end_time = time()
#             time_taken = end_time - start_time
#             feasible = thermo_feasible_mu(internal_rxn_idxs, solution_a, S)
#             return objective_value, dual_bounds, solution, time_taken, termination, iter, feasible
#         end

#         solution = round.(solution, digits=5)
#         solution_a = solution[num_reactions+1:end]

#         @assert solutions[end][num_reactions+1:end] != solution_a # ensures that solutions differ
#         # @assert sum(solution_a[O]) + sum([1-solution_a[i] for i in Z]) <= length(internal_rxn_idxs) - 1
#         @assert !(solution in solutions)
#         push!(solutions,solution)
#         push!(dual_bounds, dual_bound)
#         iter += 1
#     end

#     end_time = time()
#     time_taken = end_time - start_time
#     # @show time_taken
#     feasible = thermo_feasible_mu(internal_rxn_idxs, solution_a, S)
#     if time_taken < time_limit
#         @assert feasible
#     end
#     return objective_value, dual_bounds, solution, time_taken, termination, iter, feasible
# end

# function no_good_cuts_data(organism; time_limit=1800, csv=true)
#     model = deserialize("../molecular_models/" * organism * ".js")
#     print_model(model, "organism")

#     S = stoichiometry(model)
#     m, num_reactions = size(S)

#     lb, ub = bounds(model)
#     internal_rxn_idxs = [
#         ridx for (ridx, rid) in enumerate(variables(model)) if
#         !is_boundary(reaction_stoichiometry(model, rid))
#     ]

#     model = build_fba_model(S, lb, ub)
#     objective_value, dual_bounds, solution, time, termination, iter, thermo_feasible = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)
#     @show thermo_feasible

#     if termination == MOI.OPTIMAL
#         thermo_feasible = thermo_feasible_mu(internal_rxn_idxs, solution[num_reactions+1:end], S)
#     else 
#         thermo_feasible = false
#     end

#     df = DataFrame(
#         objective_value=objective_value, 
#         dual_bounds=[dual_bounds],
#         solution=[solution], 
#         time=time, 
#         termination=termination,
#         time_limit=time_limit, 
#         thermo_feasible=thermo_feasible,
#         iter=iter)

#     type = "no_good_cuts"
#     file_name = joinpath(@__DIR__,"../experiments/csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
#     if csv
#         CSV.write(file_name, df, append=false, writeheader=true)
#     end
# end

"""
build master problem of combinatorial Benders decomposition with FBA constraints and indicator variables,
maps indicator variables to flux direction
"""
function build_master_problem(master_problem, internal_rxn_idxs, max_flux_bound=1000; big_m=false, indicator=true)
    @show indicator, big_m
    set_attribute(master_problem, MOI.Silent(), true)
    x = master_problem[:x]

    # add indicator variables 
    a = @variable(master_problem, a[1:length(internal_rxn_idxs)], Bin)

    if big_m
        println("BIG M ADDED")
        for (cidx, ridx) in enumerate(internal_rxn_idxs)
            @constraint(master_problem, -max_flux_bound * (1 - a[cidx]) <= x[ridx])
            @constraint(master_problem, x[ridx] <= max_flux_bound * a[cidx])
        end
    end

    if indicator
        println("INDICATOR ADDED")
        # open("../csv/master_problem.lp", "w") do f
        #     print(f, master_problem)
        # end
        for (cidx, ridx) in enumerate(internal_rxn_idxs)
            # add indicator 
            @constraint(master_problem, a[cidx] => {x[ridx] >= -0.0000001})
            @constraint(master_problem, !a[cidx] => {x[ridx] <= 0.0000001})
        end
        # open("../csv/master_problem_with_binaries.lp", "w") do f
        #     print(f, master_problem)
        # end
        # write_to_file(master_problem, "../csv/models/cb_master_iAF692.mps")
    end
    
    println("model with " * string(length(x)) * " flux variables and " * string(num_variables(master_problem)) * " variables in total")
end

"""
build master problem of combinatorial Benders decomposition with FBA constraints and indicator variables
using active_on_one only
"""
function build_master_problem_complementary(master_problem, internal_rxn_idxs)
    set_attribute(master_problem, MOI.Silent(), true)
    x = master_problem[:x]

    # add indicator variables 
    a = @variable(master_problem, a[1:length(internal_rxn_idxs)], Bin)
    b = @variable(master_problem, b[1:length(internal_rxn_idxs)], Bin)
    for (cidx, ridx) in enumerate(internal_rxn_idxs)
        # add indicator 
        @constraint(master_problem, a[cidx] => {x[ridx] >= -0.000001})
        @constraint(master_problem, b[cidx] => {x[ridx] <= 0.000001})
        # complementary indicator variable
        @constraint(master_problem, b[cidx] == 1-a[cidx])
    end
    append!(a, b)
    # print(master_problem)
    return a
end

"""
build sub problem of combinatorial Benders decomposition including the thermodynamic constraints on the indicator variables 
for a given solution to the master problem and the minimal infeasible subset C
"""
function build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C_list; tol=0.000001)
    # @show solution_a
    # @show C_list

    set_attribute(sub_problem, MOI.Silent(), true)
    set_objective_sense(sub_problem, MAX_SENSE)
    S_int = Array(S[:, internal_rxn_idxs])

    G = @variable(sub_problem, G[1:length(internal_rxn_idxs)])
    μ = @variable(sub_problem, μ[1:size(S_int)[1]])
    constraints_list = []
    
    for C in C_list
        # @show C
        # solution_a = round.(solution_a, digits=6)
        constraint_list = []
        for (idx,val) in enumerate(solution_a) 
            if isapprox(val, 0, atol=tol) && (idx in C)
                c = @constraint(sub_problem, 1 <= G[idx])    
                push!(constraint_list,c)
            elseif isapprox(val, 1, atol=tol) && (idx in C)
                c = @constraint(sub_problem, G[idx] <= -1)    
                push!(constraint_list,c)
            else
                # @assert (idx in C) == false
                if (idx in C) == true
                    @warn "variable " * string(idx) * " does not have binary value: " * string(solution_a[idx])
                end
            end
        end
        push!(constraints_list, constraint_list)
    end
    # @show length(constraint_list) 
    c_matrix = @constraint(sub_problem, G .== S_int' * μ)
    # print(sub_problem)
    return constraints_list #, c_matrix
end

"""
returns a minimal infeasible subset of reactions for a given solution and stoichiometric matrix
"""
# TODO: compute several MISs at once
function compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs; fast=true, time_limit=1800, silent=true, multiple_mis=0, mis_solver=HiGHS.Optimizer)
    if !fast
        # not a MIS
        # C = [idx for (idx,val) in enumerate(solution_a) if val==1]
        C = [idx for (idx,val) in enumerate(solution_a)]
        # C = [1,2,3]
        C_list = [C]
        termination_mis = MOI.OPTIMAL
    else 
        A = deepcopy(S_int)
        for (idx,a) in enumerate(solution_a)
            if isapprox(a, 0, atol=0.000001)
                A'[idx,:] = - A'[idx,:] # update rows
            end
        end
        b = [1 for i in 1:length(internal_rxn_idxs)]

        model = Model(mis_solver)
        μ = @variable(model, μ[1:size(S_int)[1]])
        @constraint(model, A' * μ .>= b)
        @objective(model, Max, 0)

        optimizer = optimizer_with_attributes(mis_solver, "presolve" => "off")
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
        termination_mis = termination_status(mis_model)

        if termination_mis != MOI.OPTIMAL
            println("MIS problem not feasible")
            @show termination_mis
            @assert termination_mis != MOI.OTHER_ERROR
            C = []
        else
            solution_mis = [value(var) for var in all_variables(mis_model)]
            C = [idx for (idx,val) in enumerate(solution_mis) if !(isapprox(val,0))]
        end

        C_list = []
        if !isempty(C)
            push!(C_list, C)
        end  

        # @show multiple_mis
        if multiple_mis > 0 && !isempty(C)
            for i in 1:multiple_mis
                if i > multiple_mis
                    return C_list
                end
                coefs = ones(length(all_variables(mis_model)))
                coefs[i] = 0
                @objective(mis_model, Min, coefs' * all_variables(mis_model))
                optimize!(mis_model)
                termination_mis = termination_status(mis_model)

                if termination_mis != MOI.OPTIMAL
                    println("MIS problem not feasible")
                    # @show termination_status(mis_model)
                    C = []
                else
                    solution_mis = [value(var) for var in all_variables(mis_model)]
                    C = [idx for (idx,val) in enumerate(solution_mis) if !(isapprox(val,0))]
                    push!(C_list, C)
                end
            end
        end # TODO: return unique(C_list)
        # print(mis_model)
        # # λ'Aμ ≥ λ'b should be violated
        # # λ solution to sub problem, A constructed in fast MIS search, 
        # # μ flux values of master problem, b constructed in fast MIS search
        # @assert !(solution_mis' * A * solution_master[internal_rxn_idxs] >= solution_mis' * b)
    end
    # @show C_list
    # @show unique(C_list)
    return unique(C_list), termination_mis
end

"""
adds combinatorial Benders' cut to the master problem, by forcing a different assignment of the indicator variables
of the reactions in the minimal infeasible subset C.

C contains indices of internal reactions in 1:length(internal_rxn_idxs)
"""
function add_combinatorial_benders_cut(master_problem, solution_a, C_list, cuts)
    for C in C_list
        a = master_problem[:a]
        Z = []
        O = []
        for idx in C
            if solution_a[idx] > 0.000001 
                push!(O, idx)
            else 
                push!(Z, idx)
            end
        end 
        # @show Z,O
        if isempty(Z)
            c = @constraint(master_problem, sum(a[O]) <= length(C)-1)
        elseif isempty(O)
            c = @constraint(master_problem, sum([1-a[i] for i in Z]) <= length(C)-1)
        else 
            c = @constraint(master_problem, sum(a[O]) + sum([1-a[i] for i in Z]) <= length(C)-1)
        end
        push!(cuts, (O, Z, length(C)))
    end
end 

"""
adds combinatorial Benders' cut to the master problem, by forcing a different assignment of the indicator variables
of the reactions in the minimal infeasible subset C using MOI instead of JuMP
"""
function add_combinatorial_benders_cut_moi(ch, solution_a, C, a)
    # @infiltrate
    # @show a
    # m, num_reactions = size(ch.S)
    # @show solution
    # solution_a = solution[1:length(ch.internal_rxn_idxs)]
    # @show solution_a
    # solution_flux = solution[length(ch.internal_rxn_idxs)+1:length(ch.internal_rxn_idxs)+num_reactions]
    # @show solution_flux
    # @show ch.S * solution_flux == zeros(m)
    master_problem = ch.o 
    solution_a = round.(solution_a, digits=5)
    @assert !(solution_a in ch.solutions)

    # @assert !(solution_a in ch.solutions)
    push!(ch.solutions, solution_a)
    # SCIP.SCIPwriteTransProblem(
    #     ch.o,
    #     "trans_problem.lp",
    #     C_NULL,
    #     SCIP.TRUE
    # )
    # @show ch.solutions
    # @show a
    # @show solution_a
    # print(master_problem)
    no_constraints_before = MOI.get(master_problem, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}())
    Z = []
    O = []
    for idx in C
        if solution_a[idx] > 0.00001 
            push!(O,idx)
        else 
            push!(Z,idx)
        end
    end 
    # @show Z,O
    # @show a[Z], a[O]
    if isempty(Z)
        F = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(length(O)), a[O]), 0.0)
        S = MOI.LessThan(Float64(length(C)-1))
    elseif isempty(O)
        F = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(-ones(length(Z)), a[Z]), 0.0)
        S = MOI.LessThan(Float64(length(C)-1-length(Z)))
    else 
        # ATTENTION: constraint added on complementary variable v not a
        F = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(vcat(ones(length(O)), -ones(length(Z))), vcat(a[O], a[Z])), 0.0)
        S = MOI.LessThan(Float64(length(C)-1-length(Z)))
    end

    coeffs = [i.coefficient for i in F.terms]
    # @show coeffs
    # @show coeffs' * vcat(solution_a[Z], solution_a[O])
    # @show (length(C)-1-length(Z))
    # @infiltrate
    c = MOI.add_constraint(master_problem, F, S)
    push!(ch.cuts, c)
    no_constraints_after = MOI.get(master_problem, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}())
    @assert no_constraints_before < no_constraints_after

    # SCIP.SCIPwriteTransProblem(
    #     ch.o,
    #     "trans_problem_cut_added.lp",
    #     C_NULL,
    #     SCIP.TRUE
    # )
end 

"""
solve problem by splitting it into a master problem with indicator variables and a linear sub problem based 
on a solution to the master problem and minimal infeasible subsets. The sub problem 
"""
function combinatorial_benders(master_problem, internal_rxn_idxs, S, lb, ub; max_iter=Inf, fast=true, time_limit=1800, silent=true, multiple_mis=0, big_m=false, save_model=false, subproblem_solver=HiGHS.Optimizer, mis_solver=HiGHS.Optimizer, indicator=false)
    @assert indicator || big_m 

    _, num_reactions = size(S)
    start_time_cb = time()
    dual_bounds = []
    objective_values = []
    cuts = []
    times_master_problem = []
    times_sub_problem = []
    times_mis_problem = []

    @assert multiple_mis >= 0
    multiple_mis = Int(round(0.01 * multiple_mis * num_reactions))
    @show multiple_mis

    # dictionary to map internal reaction ids to index for thermodynamic feasibility variable indices
    reaction_mapping = Dict()
    for (idx, val) in enumerate(internal_rxn_idxs)
        reaction_mapping[val] = idx
    end

    # optimal_solution = parse_array_as_string(first(CSV.read("../experiments/csv/" * organism * "_combinatorial_benders_fast_600.csv", DataFrame),1)[!,:solution])
    # optimal_solution = parse_array_as_string(first(CSV.read("../experiments/csv/" * organism * "_combinatorial_benders_fast_600.csv", DataFrame),1)[!,:solution])
    
    # solve master problem
    # objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    if big_m 
        max_flux_bound = maximum(abs.(vcat(lb, ub)))
        build_master_problem(master_problem, internal_rxn_idxs, max_flux_bound, big_m=big_m, indicator=indicator)   
    else
        build_master_problem(master_problem, internal_rxn_idxs, big_m=big_m, indicator=indicator)   
    end

    # write_to_file(master_problem, "../csv/models/cb_master_iAF692.mof.json")
    if save_model
        # write_to_file(master_problem, "../experiments/csv/models/cb_master_iAF692.mof.json")
        open("../experiments/csv/model_" * organism * ".lp", "w") do f
            print(f, master_problem)
        end
    end 

    start_time = time()
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    end_time = time()
    push!(times_master_problem, end_time - start_time)
    println("master problem solved")
    # solution_master = round.(solution_master, digits=6)
    solutions = [round.(solution_master, digits=5)]
    if length(solution_master) == 1
        if isnan(solution_master) # no solution found
            @warn "master problem cannot be solved prior to adding any cuts"
            end_time = time()
            time_taken = end_time - start_time
            return NaN, NaN, NaN, NaN, time_taken, termination_master, 0, cuts
        end
    end
    solution_a = value.(master_problem[:a])
    push!(dual_bounds, dual_bound_master)
    push!(objective_values, objective_value_master)

    # compute corresponding MIS
    S_int = Array(S[:, internal_rxn_idxs])
    # @show size(S_int)
    start_time = time()
    C_list, mis_model_termination = compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs, fast=fast, time_limit=time_limit, multiple_mis=multiple_mis, mis_solver=mis_solver)
    end_time = time()
    push!(times_mis_problem, end_time - start_time)

    # build sub problem to master solution 
    optimizer = optimizer_with_attributes(subproblem_solver, "presolve" => "off")
    sub_problem = Model(subproblem_solver)
    constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C_list)
    
    start_time = time()
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
    end_time = time()
    push!(times_sub_problem, end_time - start_time)
    println("sub problem solved")

    # @show solution_a
    # @show C
    # @show termination_sub
    # @show solution_sub

    # add Benders' cut if subproblem is infeasible
    iter = 1
    while termination_sub == MOI.INFEASIBLE && iter < max_iter && time()-start_time_cb < time_limit
        @show iter
        @assert primal_status(sub_problem) == MOI.NO_SOLUTION
        # @assert dual_status(sub_problem) == MOI.INFEASIBILITY_CERTIFICATE
        # @assert has_duals(sub_problem)

        # add CB cut to MP and solve MP
        add_combinatorial_benders_cut(master_problem, solution_a, C_list, cuts)
        if save_model
            open("../experiments/csv/model_" * organism * "_" * string(iter) * ".lp", "w") do f
                print(f, master_problem)
            end
        end 

        # test if optimal solution is still feasible
        # @assert is_feasible(master_problem.moi_backend.optimizer.model, optimal_solution[1:num_reactions], optimal_solution[num_reactions+1:num_reactions+length(internal_rxn_idxs)], S, internal_rxn_idxs, cuts, lb, ub, tol=0.000001, check_thermodynamic_feasibility=false)

        # println("_______________")
        # println("master problem")
        start_time = time()
        objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
        end_time = time()
        push!(times_master_problem, end_time - start_time)
        @assert !(round.(solution_master, digits=6) in solutions)
        if termination_master != MOI.OPTIMAL
            @warn "master problem cannot be solved"
            @show termination_master
            end_time = time()
            time_taken = end_time - start_time
            return missing, objective_values, dual_bounds, missing, missing, missing, missing, missing, time_taken, termination_sub, iter, cuts, times_master_problem, times_sub_problem, times_mis_problem
        end
        push!(solutions, round.(solution_master, digits=5))
        push!(dual_bounds, dual_bound_master)
        push!(objective_values, objective_value_master)
        solution_a = value.(master_problem[:a])
        # @show solution_master

        # compute corresponding MIS
        # println("_______________")
        # println("compute MIS")
        # @show solution_a

        # check termination status of MIS computation
        start_time = time()
        C_list, mis_model_termination = compute_MIS(solution_a, S_int, solution_master, internal_rxn_idxs, fast=fast, time_limit=time_limit, silent=silent, multiple_mis=multiple_mis, mis_solver=mis_solver)
        end_time = time()
        push!(times_mis_problem, end_time - start_time)
        @show length(C_list)
        if isempty(C_list)
            if mis_model_termination != MOI.INFEASIBLE # should be TIME_LIMIT
                @show mis_model_termination
                @assert mis_model_termination == MOI.TIME_LIMIT
            else
                flux = value.(master_problem[:x])
                non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
                non_zero_flux_directions = solution_a[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
                feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
                @assert feasible
                sub_problem = Model(subproblem_solver)
                if !isinf(time_limit)
                    set_time_limit_sec(sub_problem, time_limit)
                end              
                constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, [internal_rxn_idxs])
                start_time = time()
                objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
                end_time = time()
                push!(times_sub_problem, end_time - start_time)
                @assert termination_sub == MOI.OPTIMAL
            end
        else 
            # build sub problem to master solution 
            sub_problem = Model(subproblem_solver)
            if !isinf(time_limit)
                set_time_limit_sec(sub_problem, time_limit)
            end 
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C_list)
            # println("_______________")
            # println("sub problem")
	        # write_to_file(sub_problem, "lp_models/cb_subproblem_highs.mps")
            start_time = time()
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
            end_time = time()
            push!(times_sub_problem, end_time - start_time)
            @show termination_sub
            @assert termination_sub == MOI.INFEASIBLE || termination_sub == MOI.OPTIMAL
            # @show solution_sub
            iter += 1
        end
    end

    @show iter
    end_time_cb = time()
    time_taken = end_time_cb - start_time_cb
    solution = vcat(solution_master, solution_sub)
    # @show termination_sub
    if has_values(master_problem)
        a = value.(master_problem[:a])
        x = value.(master_problem[:x])
    else 
        a = NaN 
        x = NaN
    end

    if has_values(sub_problem)
        G = value.(sub_problem[:G])
        μ = value.(sub_problem[:μ])
    else 
        G = NaN 
        μ = NaN
    end

    if termination_sub == MOI.OPTIMAL
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(x) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
        non_zero_flux_directions = a[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
        @assert feasible
    end 

    return objective_value_master, objective_values, dual_bounds, solution, x, a, G, μ, time_taken, termination_sub, iter, cuts, times_master_problem, times_sub_problem, times_mis_problem
end

function combinatorial_benders_data(organism; time_limit=1800, json=true, max_iter=Inf, fast=true, silent=true, optimizer=SCIP.Optimizer, subproblem_solver=HiGHS.Optimizer, store_optimal_solution=false, scip_tol=1.0e-6, yeast=false, multiple_mis=0, big_m=false, indicator=true)
    @show fast
    @assert multiple_mis >= 0
    if big_m + indicator == 2
        @warn "indicator and big M formulation in master problem"
    end 

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

    # dictionary to map internal reaction ids to decision variable indices of thermodynamic feasibility constraints
    reaction_mapping = Dict()
    for (idx, val) in enumerate(internal_rxn_idxs)
        reaction_mapping[val] = idx
    end

    # model = build_fba_model(S, lb, ub, optimizer=SCIP.Optimizer)
    master_problem = make_optimization_model(molecular_model, optimizer)
    if !isinf(time_limit)
        set_time_limit_sec(master_problem, time_limit)
    end    
    if optimizer == SCIP.Optimizer
        MOI.set(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"), scip_tol)
        @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    end
    # MOI.set(master_problem, MOI.RawOptimizerAttribute("presolving/maxrounds"), 0)

    objective_value, objective_values, dual_bounds, solution, x, a, G, μ, time, termination, iter, cuts, times_master_problem, times_sub_problem, times_mis_problem = combinatorial_benders(master_problem, internal_rxn_idxs, S, lb, ub; max_iter=max_iter, fast=fast, silent=silent, time_limit=time_limit, multiple_mis=multiple_mis, big_m=big_m, subproblem_solver=subproblem_solver, indicator=indicator)
    # optimal_solution = get_scip_solutions(master_problem.moi_backend.optimizer.model, number=1)
    
    if store_optimal_solution
        open(organism * "_optimal_solution.txt", "a") do io
            println(io, solution)
        end
    end

    @show termination # sub problem => thermodynamic feasibility
    @show objective_value
    # test feasibiliy
    if termination == MOI.OPTIMAL
        thermo_feasible = is_feasible(master_problem.moi_backend.optimizer.model, x, a, S, internal_rxn_idxs, cuts, lb, ub, tol=0.00001, check_indicator=false)
	    @show thermo_feasible
        # @assert thermo_feasible        
        # @assert is_feasible(master_problem.moi_backend.optimizer.model, round.(solution_flux, digits=6), solution_direction, S, internal_rxn_idxs, cuts, lb, ub, tol=0.00001)
    else 
        thermo_feasible = false 
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = objective_value
    dict[:dual_bounds] = dual_bounds
    dict[:objective_values] = objective_values
    dict[:solution] = solution
    dict[:x] = x
    dict[:a] = a
    dict[:G] = G
    dict[:μ] = μ
    dict[:time] = time
    dict[:termination] = termination
    dict[:time_limit] = time_limit
    dict[:thermo_feasible] = thermo_feasible
    dict[:iter] = iter
    dict[:scip_tol] = scip_tol
    dict[:cuts] = length(cuts)
    dict[:times_master_problem] = times_master_problem
    dict[:times_sub_problem] = times_sub_problem
    dict[:times_mis_problem] = times_mis_problem

    type = "combinatorial_benders"
    if fast
        type = type * "_fast"
    end
    if indicator && big_m
        type = type * "_indicator_and"
    end
    if big_m
        type = type * "_big_m"
    end
    if multiple_mis != 0
        type = type * "_" * string(multiple_mis) * "_mis"
    end
    if optimizer != SCIP.Optimizer
        type = type * "_" * replace(string(optimizer), ".Optimizer"=>"")
    end
    file_name = "json/" * organism * "_" * type * "_" * string(time_limit) * ".json"
    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
end

"""
check whether solution is feasible, within a tolerance
"""
function is_feasible(o, solution_flux, solution_direction, S, internal_rxn_idxs, cuts, lb, ub; check_steady_state=true, check_bounds=true, check_thermodynamic_feasibility=true, check_cuts=true, check_indicator=true, tol=0.000001)
    m, num_reactions = size(S)
    # check steady state assumption 
    if check_steady_state
        steady_state = collect(S) * solution_flux
        # @show steady_state
        for (idx, s) in enumerate(steady_state) 
            if !isapprox(s, 0, atol=tol) 
                @show idx, s
                println("solution does not respect steady state assumption")
                return false
            end
        end
    end
    # check bounds 
    if check_bounds
        if !solution_within_bounds(solution_flux, lb, ub, tol=tol)
            println("solution does not respect reaction bounds")
            return false
        end
    end
    # check thermo feasiblity 
    if check_thermodynamic_feasibility
        reaction_mapping = Dict()
        for (idx, val) in enumerate(internal_rxn_idxs)
            reaction_mapping[val] = idx
        end
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution_flux) if !isapprox(val, 0, atol=1e-6)], internal_rxn_idxs)
        non_zero_flux_directions = solution_direction[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
        feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
        if !feasible 
            println("solution is not thermodynamically feasible")
            return false 
        end
    end
    # check Benders' cuts 
    if check_cuts
        for (idx, (O, Z, length_C)) in enumerate(cuts)
            # @show O, Z, length_C
            if isempty(Z)
                # @show sum(solution_direction[O]), length_C-1
                feasible = sum(solution_direction[O]) <= length_C-1 + tol
            elseif isempty(O)
                # @show sum([1-solution_direction[i] for i in Z]), length_C-1
                feasible = sum([1-solution_direction[i] for i in Z]) <= length_C-1 + tol
            else 
                # @show sum(solution_direction[O]) + sum([1-solution_direction[i] for i in Z]), length_C-1
                # @show solution_direction[O]
                # @show solution_direction[Z]
                feasible = sum(solution_direction[O]) + sum([1-solution_direction[i] for i in Z]) <= length_C-1 + tol
            end
            if !feasible             
                println("solution does not respect Benders' cut")
                return false
            end
        end
    end
    # check indicator variables
    if check_indicator
        solution_flux_internal_rxns = solution_flux[internal_rxn_idxs]
        # @show solution_flux_internal_rxns, solution_direction
        for (idx, a) in enumerate(solution_direction)
            if isapprox(a, 1, atol=tol)
                if solution_flux_internal_rxns[idx] < 0 + tol
                    @show solution_flux_internal_rxns[idx]
                    println("solution does not respect indicator constraint")
                    return false 
                end 
            elseif isapprox(a, 0, atol=tol)
                if solution_flux_internal_rxns[idx] > 0 - tol
                    @show solution_flux_internal_rxns[idx]
                    println("solution does not respect indicator constraint")
                    return false 
                end 
            end 
        end
    end
    return true
end
