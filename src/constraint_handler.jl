using SCIP
include("cuts_decomposition.jl")

mutable struct ThermoFeasibleConstaintHandler{} <: SCIP.AbstractConstraintHandler
    o
    ncalls::Int
    internal_rxn_idxs::Vector{Int64} 
    S::Matrix{Int64}
    vars    
    binvars
end

# check if primal solution candidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    println("CHECK")
    @show ch.binvars
    @show ch.vars
    @show sol
    solution_flux = SCIP.sol_values(ch.o, ch.vars)
    solution_direction = SCIP.sol_values(ch.o, ch.binvars)
    append!(solution_flux, solution_direction)
    @show solution_flux
    feasible = thermo_feasible_mu(ch.internal_rxn_idxs, solution_flux[ch.internal_rxn_idxs], ch.S)
    @show thermo_feasible_mu(ch.internal_rxn_idxs, solution_flux[ch.internal_rxn_idxs], ch.S)
    if !feasible      
        return SCIP.SCIP_INFEASIBLE
    end
    return SCIP.SCIP_FEASIBLE
    
end

function SCIP.enforce_lp_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible)
    println("LP SOL")
    # @assert length(constraints) == 0
    add_cb_cut(ch)
    return SCIP.SCIP_CONSADDED
end

function SCIP.enforce_pseudo_sol(
        ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss,
        solinfeasible, objinfeasible,
    )
    println("PSEUDO SOL")
    @assert length(constraints) == 0
    add_cb_cut(ch)
    return SCIP.SCIP_CONSADDED
end

# add combinatorial Benders cut if solution infeasible
function add_cb_cut(ch::ThermoFeasibleConstaintHandler)
    println("BENDERS CUT")
    solution = SCIP.sol_values(ch.o, ch.vars) # TODO: solution of a and x
    internal_rxn_idxs = ch.internal_rxn_idxs
    S = ch.S 

    feasible = thermo_feasible_mu(internal_rxn_idxs, solution[internal_rxn_idxs], S)

    if !feasible
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
            add_combinatorial_benders_cut(master_problem, solution_a, C)
            ch.ncalls += 1
        else 
            # build sub problem to master solution 
            sub_problem = Model(optimizer)
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_a, C)
            println("_______________")
            println("sub problem")
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=silent, time_limit=time_limit)
            return SCIP.SCIP_FEASIBLE
        end
        return SCIP.SCIP_CONSADDED
    end
    return SCIP.SCIP_FEASIBLE
end

# lock binary variables of master problem
function SCIP.lock(ch::ThermoFeasibleConstaintHandler, constraint, locktype, nlockspos, nlocksneg)
    println("LOCK")
    for a in ch.binvars
        ai::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.o, a)
        # @show ai
        ai == C_NULL && continue
        SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(ch.o, ai, SCIP.SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg)
    end
end

