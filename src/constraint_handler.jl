using SCIP
using CSV, DataFrames

include("cuts_decomposition.jl")

mutable struct ThermoFeasibleConstaintHandler{} <: SCIP.AbstractConstraintHandler
    o
    ncalls::Int
    internal_rxn_idxs::Vector{Int64} 
    S::Matrix{Float64}
    vars    
    binvars
    feasible_solutions
    solutions
    cuts
    multiple_mis
end

# check if primal solution candidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    println("CHECK")

    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    solution_flux = SCIP.sol_values(ch.o, ch.vars)

    S_int = ch.S[:, ch.internal_rxn_idxs]

    # assert that the variables are binary 
    non_binary_values = [i for i in solution_direction if !(isapprox(i, 0, atol=1e-4) || isapprox(i, 1, atol=1e-4))]
    if !isempty(non_binary_values)
        @warn "binary variables do not have binary value assigned"
        @show maximum(non_binary_values), minimum(non_binary_values)
        return SCIP.SCIP_INFEASIBLE
    end 

    # build sub problem to master solution 
    C_list, termination_mis = compute_MIS(solution_direction, S_int, [], ch.internal_rxn_idxs, fast=true)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, C_list)
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    @show termination_sub

    thermo_feasible = thermo_feasible_mu(ch.internal_rxn_idxs, solution_direction, ch.S)
    @show thermo_feasible
    if thermo_feasible
        @assert termination_sub == MOI.OPTIMAL
        push!(ch.solutions, vcat(solution_direction, solution_flux))
        return SCIP.SCIP_FEASIBLE
    else 
        return SCIP.SCIP_INFEASIBLE
    end 
end

function SCIP.enforce_lp_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible)
    println("LP SOL")

    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    S_int = ch.S[:, ch.internal_rxn_idxs]

    # assert that the variables are binary 
    non_binary_values = [i for i in solution_direction if !(isapprox(i, 0, atol=1e-4) || isapprox(i, 1, atol=1e-4))]
    if !isempty(non_binary_values)
        @warn "binary variables do not have binary value assigned"
        @show maximum(non_binary_values), minimum(non_binary_values)
        return SCIP.SCIP_INFEASIBLE
    end 

    # build sub problem to master solution 
    C_list, termination_mis = compute_MIS(solution_direction, S_int, [], ch.internal_rxn_idxs, fast=true)
    @show termination_mis
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, C_list)
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    @show termination_sub

    # if thermo_feasible_mu(ch.internal_rxn_idxs, solution_direction, ch.S)
    #     @assert termination_sub == MOI.OPTIMAL
    # end 

    # if isempty(C_list)
    #     optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    #     sub_problem = Model(optimizer)
    #     build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, ch.internal_rxn_idxs)
    #     objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    #     @assert termination_sub == MOI.OPTIMAL
    # end

    if termination_sub == MOI.OPTIMAL 
        return SCIP.SCIP_FEASIBLE
    else 
        # add CB cut to master problem
        SCIP_status = add_cb_cut(ch, solution_direction, C_list)
        return SCIP_status
    end
end

function SCIP.enforce_pseudo_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible, objinfeasible)
    println("PSEUDO SOL")
    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    S_int = ch.S[:, ch.internal_rxn_idxs]

    # assert that the variables are binary 
    if !isempty([1 for i in solution_direction if !(isapprox(i, 0, atol=1e-4) || isapprox(i, 1, atol=1e-4))])
        @warn "binary variables do not have binary value assigned"
        return SCIP.SCIP_INFEASIBLE
    end 

    # build sub problem to master solution 
    C_list, termination_mis = compute_MIS(solution_direction, S_int, [], ch.internal_rxn_idxs, fast=true, multiple_mis=ch.multiple_mis)
    @show termination_mis
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, C_list)
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    @show termination_sub

    # if thermo_feasible_mu(ch.internal_rxn_idxs, solution_direction, ch.S)
    #     @assert termination_sub == MOI.OPTIMAL
    # end 

    # if isempty(C_list)
    #     optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    #     sub_problem = Model(optimizer)
    #     build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, ch.internal_rxn_idxs)
    #     objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    #     @assert termination_sub == MOI.OPTIMAL
    # end

    if termination_sub == MOI.OPTIMAL 
        return SCIP.SCIP_FEASIBLE
    else 
        # add CB cut to master problem
        SCIP_status = add_cb_cut(ch, solution_direction, C_list)
        return SCIP_status
    end
end

# add combinatorial Benders cut if solution infeasible
function add_cb_cut(ch::ThermoFeasibleConstaintHandler, solution_direction, C_list)
    println("BENDERS CUT")

    no_constraints_before = length(ch.cuts)

    add_combinatorial_benders_cut_moi(ch, solution_direction, C_list, ch.binvars[1:length(ch.internal_rxn_idxs)])
    ch.ncalls += 1

    no_constraints_after = length(ch.cuts)

    @show no_constraints_before, no_constraints_after
    @show ch.cuts
    # check efficacy 
    for c_idx in no_constraints_before+1:no_constraints_after-no_constraints_before
        a = [t.coefficient for t in MOI.get(ch.o, MOI.ConstraintFunction(), ch.cuts[c_idx]).terms]
        x = [solution_direction[i] for i in C_list[c_idx]]
        @assert MOI.get(ch.o, MOI.ConstraintFunction(), ch.cuts[c_idx]).constant == 0
        b = MOI.get(ch.o, MOI.ConstraintSet(), ch.cuts[c_idx]).upper
        @show a, x, b
        @show sum(a.*x), b
        @assert (sum(a.*x) -b >= 1e-3)
    end

    return SCIP.SCIP_CONSADDED
end

# lock binary variables of master problem
# thermodynamic feasibility deals with the direction only (represented by binary variables)
# variables should not be changed in any direction: SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg)
function SCIP.lock(ch::ThermoFeasibleConstaintHandler, constraint, locktype, nlockspos, nlocksneg)
    println("LOCK")

    for a in ch.binvars
        ai::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.o, a)
        ai == C_NULL && continue
        SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(ch.o, ai, SCIP.SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end

function constraint_handler_data(organism; time_limit=1800, json=true, silent=true, mute=true, big_m=true, multiple_mis=0)
    molecular_model = load_model("../molecular_models/" * organism * ".json")
    print_model(molecular_model, organism)

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model)
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    
    # select percentage of reactions
    @assert multiple_mis >= 0
    multiple_mis_no = Int(round(0.01 * multiple_mis * num_reactions))
    @show multiple_mis_no

    reaction_mapping = Dict()
    for (idx, val) in enumerate(internal_rxn_idxs)
        reaction_mapping[val] = idx
    end

    # extract objective
    objective_func = objective_function(model)
    objective_func_vars = [i.index for i in objective_func.terms.keys]

    # optimize_model(model, mute=false)

    scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, time_limit=time_limit, objective_func_vars=objective_func_vars, objective_func_coeffs=objective_func.terms.vals, silent=silent, big_m=big_m)
    # print(scip_model)

    MOI.set(scip_model, MOI.Silent(), false)

    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [], multiple_mis_no)
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch", enforce_priority=-7000000, check_priority=-7000000)
    MOI.optimize!(scip_model)
    
    # # filter constraint coefficients and vars
    # constraints = [MOI.get(ch.o, MOI.ConstraintFunction(), i) for i in ch.cuts]
    # constraints_tuples = []
    # for c in constraints 
    #     @show c
    #     constraint_tuples = []
    #     for t in c.terms
    #         push!(constraint_tuples, (t.coefficient, t.variable))
    #     end 
    #     push!(constraints_tuples, constraint_tuples)
    # end 

    status = MOI.get(scip_model, MOI.TerminationStatus())
    time = MOI.get(scip_model, MOI.SolveTimeSec())
    result_status = MOI.get(scip_model, MOI.PrimalStatus())
    @show status, result_status

    if result_status != MOI.NO_SOLUTION
        primal_objective_value = MOI.get(scip_model, MOI.ObjectiveValue())
        # TODO: get correct dual bound
        dual_objective_value = MOI.get(scip_model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=2))
            println("")
        end
        # solution = [MOI.get(ch.o, MOI.VariablePrimal(1), MOI.VariableIndex(i)) for i in 1:length(internal_rxn_idxs) + num_reactions]
        solution_direction = SCIP.sol_values(ch.o, ch.binvars)
        solution_flux = SCIP.sol_values(ch.o, ch.vars)
        solution = vcat(solution_flux, solution_direction)

        # test feasibility, filter non-zero fluxes, set binaries accordingly
        non_zero_flux_indices = intersect([idx for (idx, val) in enumerate(solution_flux) if !isapprox(val, 0, atol=1e-5)], internal_rxn_idxs)
        non_zero_flux_directions = solution_direction[collect(reaction_mapping[val] for val in non_zero_flux_indices)] 
        thermo_feasible = thermo_feasible_mu(non_zero_flux_indices, non_zero_flux_directions, S; scip_tol=0.001)
        
        if result_status == MOI.OPTIMAL
            @assert thermo_feasible
            
            # verify that flux directions are binary
            @assert isempty([1 for i in solution_direction if !(isapprox(i, 0, atol=1e-4) || isapprox(i, 1, atol=1e-4))])
            
            # test steady state constraint 
            steady_state =  isapprox.(S * solution_flux, 0, atol=0.0001)
            @assert steady_state == ones(size(S)[1])
        end
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
        solution_direction = NaN
        solution_flux = NaN
        thermo_feasible = missing
    end

    dict = Dict{Symbol, Any}()
    dict[:objective_value] = primal_objective_value 
    dict[:dual_bound] = dual_objective_value
    dict[:solution] = [solution] 
    dict[:x] = solution_flux
    dict[:a] = solution_direction
    dict[:time] = time 
    dict[:termination] = status
    dict[:time_limit] = time_limit 
    dict[:ncalls] = ch.ncalls
    dict[:thermo_feasible] = thermo_feasible

    type = "constraint_handler"
    if !big_m
        tyoe = type * "_indicator"
    end

    if multiple_mis != 0
        type = type * "_" * string(multiple_mis) * "_mis"
    end 

    file_name = joinpath(@__DIR__,"../experiments/json/" * organism * "_" * type * "_" * string(time_limit) * ".json")
    if json 
        open(file_name, "w") do f
            JSON.print(f, dict) 
        end
    end
end

"""
analyze model with constraint handler after solve
"""
function check_solutions(ch, lb, ub)
    # check for optimal feasible solution
    feasible_solutions = [sol for (idx,sol) in enumerate(ch.feasible_solutions) if [solution_within_bounds(sol, lb, ub) for sol in ch.feasible_solutions][idx]]
    optimal_solution_objective_value, optimal_solution_idx  = findmax([eval_objective(ch, sol) for sol in feasible_solutions])
    optimal_solution = feasible_solutions[optimal_solution_idx]
    feasible = thermo_feasible(ch.internal_rxn_idxs, optimal_solution, ch.S)
    @assert feasible
    return optimal_solution_objective_value, optimal_solution
end 

"""
check whether solution is feasible in SCIP directly
"""
function is_feasible_scip(model, solution)
    # build solution
    vars = SCIP.SCIPgetOrigVars(model)
    # @show vars
    # nvars = MOI.get(model, MOI.NumberOfVariables())
    # @show nvars
    nvars = SCIP.SCIPgetNOrigVars(model)
    # @show nvars
    var_vec = unsafe_wrap(Array, vars, nvars)
    sol = SCIP.create_empty_scipsol(model.inner.scip[], C_NULL)

    # ATTENTION: SCIPgetVars != MOI.NumberOfVariables
    for (idx, var) in enumerate(var_vec)
        SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(
            model,
            sol,
            var,
            solution[idx], 
        )
    end

    # retcode = SCIP.SCIPwriteTransProblem(model, "scip.lp", C_NULL, SCIP.FALSE)
    # retcode = SCIP.SCIPwriteOrigProblem(model, "scip.lp", C_NULL, SCIP.FALSE)

    # @infiltrate

    # WARNING: calls SCIP.check again and leads to infinity loop
    # checks solution for feasibility in original problem without adding it to the solution store
    feasible = Ref{SCIP.SCIP_Bool}(SCIP.TRUE)
    SCIP.SCIPcheckSolOrig(
        model, 
        sol, 
        feasible,
        SCIP.TRUE, 
        SCIP.TRUE,
    )

    # # checks solution for feasibility without adding it to the solution store
    feasible = Ref{SCIP.SCIP_Bool}(SCIP.TRUE)
    SCIP.SCIPcheckSol(
        model, 
        sol, 
        SCIP.TRUE, 
        SCIP.TRUE,
        SCIP.TRUE,
        SCIP.TRUE,
        SCIP.TRUE,
        feasible
    )
    @show feasible[]
end

"""
compute the objective value for a solution
"""
function eval_objective(ch, solution)
    f = MOI.get(ch.o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    coeffs = [i.coefficient for i in f.terms]
    var_idxs = [i.variable.value for i in f.terms]
    # @show coeffs, var_idxs
    return coeffs' * solution[var_idxs .- length(ch.binvars)]
end

"""
add solution to SCIP optimizer
"""
# BUG: variable order does not match the order of SCIP's stored solutions
function add_solution(model, solution)
    println("ADD SOLUTION")
    @show SCIP.SCIPgetStage(model)
    if SCIP.SCIPgetStage(model) != SCIP.LibSCIP.SCIP_STAGE_SOLVED
        # @infiltrate
        get_scip_solutions(model)

        solution = round.(solution, digits=5)
        @show solution

        # build solution
        vars = SCIP.SCIPgetOrigVars(model)
        nvars = MOI.get(model, MOI.NumberOfVariables()) 
        # nvars = SCIP.SCIPgetNOrigVars(model)
        var_vec = unsafe_wrap(Array, vars, nvars)
        sol = SCIP.create_empty_scipsol(model.inner.scip[], C_NULL)

        # ATTENTION: SCIPgetVars != MOI.NumberOfVariables
        for (idx, var) in enumerate(var_vec)
            SCIP.@SCIP_CALL SCIP.SCIPsetSolVal(
                model,
                sol,
                var,
                solution[idx], 
            )
        end

        # submit solution
        stored = Ref{SCIP.SCIP_Bool}(SCIP.FALSE)
        SCIP.@SCIP_CALL SCIP.SCIPtrySolFree(
            model,
            Ref(sol),
            SCIP.TRUE,
            SCIP.FALSE,
            SCIP.FALSE,
            SCIP.FALSE,
            SCIP.FALSE,
            stored,
        )
        @show stored[] # 0x00000001
        get_scip_solutions(model)
    end
end

"""
show all scip solutions
"""
function get_scip_solutions(o::SCIP.Optimizer; number=Inf)
    vars = [MOI.VariableIndex(i) for i in 1:MOI.get(o, MOI.NumberOfVariables())]
    sols_vec =
        unsafe_wrap(Vector{Ptr{Cvoid}}, SCIP.LibSCIP.SCIPgetSols(o), SCIP.LibSCIP.SCIPgetNSols(o))
    solutions = []
    if isempty(sols_vec)
        println("no solutions stored")
    else 
        solutions = [SCIP.sol_values(o, vars, sol) for sol in sols_vec]
        solutions = [round.(sol, digits=5) for sol in solutions]
    end
    
    if !isinf(number)
        return solutions[1]
    end
    return solutions
end
