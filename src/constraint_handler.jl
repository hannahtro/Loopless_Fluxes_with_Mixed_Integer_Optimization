using SCIP
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
end

# check if primal solution candidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    println("CHECK")
    # check sub problem for binary variables
    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    # build sub problem to master solution 
    C = compute_MIS(solution_direction, (ch.S[:, ch.internal_rxn_idxs]), [], ch.internal_rxn_idxs)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, C)
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    if termination_sub == MOI.OPTIMAL 
        return SCIP.SCIP_FEASIBLE
    else 
        return SCIP.SCIP_INFEASIBLE
    end
end

function SCIP.enforce_lp_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible)
    println("LP SOL")
    # check sub problem for binary variables and add CB cut to master problem
    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    # build sub problem to master solution 
    C = compute_MIS(solution_direction, (ch.S[:, ch.internal_rxn_idxs]), [], ch.internal_rxn_idxs)
    optimizer = optimizer_with_attributes(HiGHS.Optimizer, "presolve" => "off")
    sub_problem = Model(optimizer)
    build_sub_problem(sub_problem, ch.internal_rxn_idxs, ch.S, solution_direction, C)
    objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem)
    if termination_sub == MOI.OPTIMAL 
        return SCIP.SCIP_FEASIBLE
    else 
        add_cb_cut(ch, solution_direction)
        return SCIP.SCIP_CONSADDED
    end

end

# function SCIP.enforce_pseudo_sol(
#         ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss,
#         solinfeasible, objinfeasible,
#     )
#     println("PSEUDO SOL")
#     solution = SCIP.sol_values(ch.o, ch.vars)
#     @show solution
#     feasible = thermo_feasible_mu(ch.internal_rxn_idxs, round.(solution[ch.internal_rxn_idxs],digits=5), ch.S)
#     # @show solution[ch.internal_rxn_idxs]
#     # @assert length(constraints) == 0
#     if !feasible      
#         add_cb_cut(ch)
#         return SCIP.SCIP_CONSADDED
#     end 
#     return SCIP.SCIP_FEASIBLE
# end

# add combinatorial Benders cut if solution infeasible
function add_cb_cut(ch::ThermoFeasibleConstaintHandler, solution_direction)
    println("BENDERS CUT")
    C = compute_MIS(solution_direction, (ch.S[:, ch.internal_rxn_idxs]), [], ch.internal_rxn_idxs)
    add_combinatorial_benders_cut_moi(ch, solution_direction, C, ch.binvars[1:length(ch.internal_rxn_idxs)])
    ch.ncalls += 1
end

# lock binary variables of master problem
function SCIP.lock(ch::ThermoFeasibleConstaintHandler, constraint, locktype, nlockspos, nlocksneg)
    println("LOCK")
    for a in ch.binvars
        ai::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.o, a)
        ai == C_NULL && continue
        SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(ch.o, ai, SCIP.SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end

function constraint_handler_data(organism; time_limit=1800, csv=true, silent=true, mute=true)
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, "organism")

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    
    # extract objective
    objective_func = objective_function(model)
    objective_func_vars = [i.index for i in objective_func.terms.keys]

    scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, time_limit=time_limit, objective_func_vars=objective_func_vars, objective_func_coeffs=objective_func.terms.vals, silent=silent)
    # print(scip_model)
    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [])
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch")
    MOI.optimize!(scip_model)

    status = MOI.get(scip_model, MOI.TerminationStatus())
    time = MOI.get(scip_model, MOI.SolveTimeSec())
    result_status = MOI.get(scip_model, MOI.PrimalStatus())
    @show status, result_status
    if result_status != MOI.NO_SOLUTION
        primal_objective_value = MOI.get(scip_model, MOI.ObjectiveValue())
        solution = [MOI.get(ch.o, MOI.VariablePrimal(1), MOI.VariableIndex(i)) for i in 1:length(internal_rxn_idxs) + num_reactions]
        
        # TODO: get correct dual bound
        dual_objective_value = MOI.get(scip_model, MOI.ObjectiveBound())
        if !mute
            println("objective value : ", round(primal_objective_value, digits=2))
            println("")
        end
        # solution = MOI.get(scip_model, MOI.VariablePrimal(), vcat(ch.vars, ch.binvars))
    else 
        primal_objective_value = NaN
        dual_objective_value = NaN
        solution = NaN
    end

    df = DataFrame(
        objective_value=primal_objective_value, 
        dual_bound=dual_objective_value,
        solution=[solution], 
        time=time, 
        termination=status,
        time_limit=time_limit, 
        ncalls=ch.ncalls)

    type = "constraint_handler"
    file_name = joinpath(@__DIR__,"../csv/" * organism * "_" * type * "_" * string(time_limit) * ".csv")
    if csv
        CSV.write(file_name, df, append=false, writeheader=true)
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
check whether a solution respects the upper and lower bounds of a reaction
"""
function solution_within_bounds(solution, lb, ub; tol=0.00001)
    for (idx,sol) in enumerate(solution)
        if sol < lb[idx] - tol || sol > ub[idx] + tol
            return false
        end
    end
    return true
end

"""
check whether solution is feasible, within a tolerance
"""
# TODO: add tolerance to different checks
function is_feasible(ch, solution_flux, solution_direction; check_steady_state=true, check_bounds=true, check_thermodynamic_feasibility=true, check_cuts=true, check_indicator=true, tol=0.00001)
    # check steady state assumption 
    if check_steady_state
        steady_state = ch.S * solution_flux
        # @show steady_state
        for s in steady_state 
            if !isapprox(s, 0, atol=tol) 
                println("solution does not respect steady state assumption")
                return false
            end
        end
    end
    # check bounds 
    if check_bounds
        if !solution_within_bounds(solution_flux, ch.lb, ch.ub, tol=tol)
            println("solution does not respect reaction bounds")
            return false
        end
    end
    # check thermo feasiblity 
    if check_thermodynamic_feasibility
        feasible = thermo_feasible_mu(ch.internal_rxn_idxs, round.(solution_flux[ch.internal_rxn_idxs],digits=5), ch.S)
        if !feasible 
            println("solution is not thermodynamically feasible")
            return false 
        end
    end
    # check Benders' cuts 
    if check_cuts
        for cut in ch.cuts
            upper = MOI.get(ch.o, MOI.ConstraintSet(), cut).upper
            func = MOI.get(ch.o, MOI.ConstraintFunction(), cut)
            var_indices = [term.variable.value for term in func.terms]
            coefficients = [term.coefficient for term in func.terms]
            if !(coefficients' * solution_direction[var_indices] + func.constant <= upper)
                println("solution does not respect Benders' cut")
                return false
            end
        end
    end
    # check indicator variables
    if check_indicator
        solution_flux_internal_rxns = solution_flux[ch.internal_rxn_idxs]
        @show solution_flux_internal_rxns, solution_direction
        for (idx, a) in enumerate(solution_direction)
            if isapprox(a, 1, atol=tol)
                if solution_flux_internal_rxns[idx] < 0 - tol
                    println("solution does not respect indicator constraint")
                    return false 
                end 
            elseif isapprox(a, 0, atol=tol)
                if solution_flux_internal_rxns[idx] > 0 + tol
                    println("solution does not respect indicator constraint")
                    return false 
                end 
            end 
        end
    end
    return true
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
function get_scip_solutions(o::SCIP.Optimizer)
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
    @show solutions
    return solutions
end
