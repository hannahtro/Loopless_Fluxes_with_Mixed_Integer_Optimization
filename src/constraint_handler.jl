using SCIP
include("cuts_decomposition.jl")

mutable struct ThermoFeasibleConstaintHandler{} <: SCIP.AbstractConstraintHandler
    o
    ncalls::Int
    internal_rxn_idxs::Vector{Int64} 
    S::Matrix{Float64}
    vars    
    binvars
    solutions
    feasible_solutions
end

# check if primal solution candidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    println("CHECK")
    solution = SCIP.sol_values(ch.o, ch.vars)
    # solution_master_flux = round.(SCIP.sol_values(ch.o, ch.vars),digits=5)

    @show solution[ch.internal_rxn_idxs]
    # @show SCIP.sol_values(ch.o, ch.binvars)
    feasible = thermo_feasible_mu(ch.internal_rxn_idxs, round.(solution[ch.internal_rxn_idxs],digits=5), ch.S)
    @show feasible
    if !feasible      
        return SCIP.SCIP_INFEASIBLE
    end
    # check if solution is feasible
    @infiltrate
    model_temp = Model(SCIP.Optimizer)
    MOI.copy_to(model_temp, ch.o)
    cont_vars = Dict(VariableRef(model_temp, ch.vars[i])  => solution[i] for i in 1:length(ch.vars))
    bin_vars = Dict(VariableRef(model_temp, ch.binvars[i]) => SCIP.sol_values(ch.o, ch.binvars)[i] for i in 1:length(ch.binvars))
    point = merge(cont_vars, bin_vars)
    @show primal_feasibility_report(model_temp, point, skip_missing=true)
    push!(ch.feasible_solutions, solution)
    @show ch.S * solution
    return SCIP.SCIP_FEASIBLE
end

function SCIP.enforce_lp_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible)
    println("LP SOL")
    solution = SCIP.sol_values(ch.o, ch.vars)
    @show solution
    feasible = thermo_feasible_mu(ch.internal_rxn_idxs, round.(solution[ch.internal_rxn_idxs],digits=5), ch.S)
    # @show solution[ch.internal_rxn_idxs]
    # @assert length(constraints) == 0
    @show feasible
    # @show MOI.get(ch.o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    if !feasible      
        add_cb_cut(ch)
        return SCIP.SCIP_CONSADDED
    end 
    push!(ch.feasible_solutions,solution)
    @show ch.S * solution
    return SCIP.SCIP_FEASIBLE
end

function SCIP.enforce_pseudo_sol(
        ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss,
        solinfeasible, objinfeasible,
    )
    println("PSEUDO SOL")
    solution = SCIP.sol_values(ch.o, ch.vars)
    @show solution
    feasible = thermo_feasible_mu(ch.internal_rxn_idxs, round.(solution[ch.internal_rxn_idxs],digits=5), ch.S)
    # @show solution[ch.internal_rxn_idxs]
    # @assert length(constraints) == 0
    if !feasible      
        add_cb_cut(ch)
        return SCIP.SCIP_CONSADDED
    end 
    return SCIP.SCIP_FEASIBLE
end

# add combinatorial Benders cut if solution infeasible
function add_cb_cut(ch::ThermoFeasibleConstaintHandler)
    println("BENDERS CUT")
    internal_rxn_idxs = ch.internal_rxn_idxs
    S = ch.S 
    S_int = Array(S[:, internal_rxn_idxs])

    solution_master_flux = round.(SCIP.sol_values(ch.o, ch.vars),digits=5)
    solution_master_direction = SCIP.sol_values(ch.o, ch.binvars)[1:length(internal_rxn_idxs)]
    solution_master = SCIP.sol_values(ch.o, [MOI.VariableIndex(i) for i in 1:MOI.get(ch.o, MOI.NumberOfVariables())])
    # @show solution
    @show solution_master_flux
    @show solution_master_direction

    feasible = thermo_feasible_mu(internal_rxn_idxs, solution_master_flux[internal_rxn_idxs], S)

    @assert !feasible
    # compute corresponding MIS
    # println("_______________")
    # println("compute MIS")
    m, num_reactions = size(S)
    # solution_a = solution_master[num_reactions+1:end]
    C = compute_MIS(solution_master_direction, S_int, [], internal_rxn_idxs, fast=true, time_limit=600, silent=true)
    # @show ch.S * solution_master_flux
    @show C
    if isempty(C)
        @infiltrate
    else 
        # build sub problem to master solution 
        sub_problem = Model(SCIP.Optimizer)
        constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_master_direction, C)
        # println("_______________")
        # println("sub problem")
        # print(sub_problem)
        objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=true, time_limit=600)
        add_combinatorial_benders_cut_moi(ch, solution_master, C, ch.binvars[1:length(internal_rxn_idxs)])
        ch.ncalls += 1
    end
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
    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [])
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch")
    MOI.optimize!(scip_model)

    status = MOI.get(scip_model, MOI.TerminationStatus())
    time = MOI.get(scip_model, MOI.SolveTimeSec())
    result_status = MOI.get(scip_model, MOI.PrimalStatus())
    @show status, result_status
    if result_status != MOI.NO_SOLUTION
        primal_objective_value, solution = check_solutions(ch, lb, ub)
        # primal_objective_value = MOI.get(scip_model, MOI.ObjectiveValue())
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
function solution_within_bounds(solution, lb, ub)
    for (idx,sol) in enumerate(solution)
        if sol < lb[idx] || sol > ub[idx]
            return false
        end
    end
    return true
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