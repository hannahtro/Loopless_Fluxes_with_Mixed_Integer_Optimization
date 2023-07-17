using SCIP
include("cuts_decomposition.jl")

mutable struct ThermoFeasibleConstaintHandler{} <: SCIP.AbstractConstraintHandler
    o
    ncalls::Int
    internal_rxn_idxs::Vector{Int64} 
    S::Matrix{Int64}
end

# check if primal solution cnadidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    print(ch.o)
    # @show sol
    # @show sol[ch.internal_rxn_idxs]
    # feasible = thermo_feasible_mu(ch.internal_rxn_idxs, sol[ch.internal_rxn_idxs], S)
    feasible = false
    if !feasible
        return SCIP.SCIP_INFEASIBLE
    end
    return SCIP.SCIP_FEASIBLE
end

function SCIP.enforce_lp_sol(ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss, solinfeasible)
    # @assert length(constraints) == 0
    add_cb_cut(ch)
    return SCIP_CONSADDED
end

function SCIP.enforce_pseudo_sol(
        ch::ThermoFeasibleConstaintHandler, constraints, nusefulconss,
        solinfeasible, objinfeasible,
    )
    @assert length(constraints) == 0
    add_cb_cut(ch)
    return SCIP_CONSADDED
end

# add combinatorial Benders cut if solution infeasible
function add_cb_cut(ch::ThermoFeasibleConstaintHandler)
    solution = SCIP.sol_values(ch.o, ch.vars)
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

# TODO:lock binary variables of master problem
function SCIP.lock(ch::ThermoFeasibleConstaintHandler, constraint, locktype, nlockspos, nlocksneg)
    # z::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.o, ch.epivar)
    # if z != C_NULL
    #     SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(ch.o, z, SCIP.SCIP_LOCKTYPE_MODEL, nlockspos, nlocksneg)
    # end
    # for x in ch.vars
    #     xi::Ptr{SCIP.SCIP_VAR} = SCIP.var(ch.o, x)
    #     xi == C_NULL && continue
    #     SCIP.@SCIP_CALL SCIP.SCIPaddVarLocksType(ch.o, xi, SCIP.SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg)
    # end
end

