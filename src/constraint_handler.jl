using SCIP
include("cuts_decomposition.jl")

mutable struct ThermoFeasibleConstaintHandler{} <: SCIP.AbstractConstraintHandler
    o
    ncalls::Int
    internal_rxn_idxs::Vector{Int64} 
    S::Matrix{Int64}
    vars    
    binvars
    solutions
end

# check if primal solution candidate is thermodynamically feasible
function SCIP.check(ch::ThermoFeasibleConstaintHandler, constraints::Vector{Ptr{SCIP.SCIP_CONS}}, sol::Ptr{SCIP.SCIP_SOL}, checkintegrality::Bool, checklprows::Bool, printreason::Bool, completely::Bool; tol=1e-6)
    println("CHECK")
    solution = SCIP.sol_values(ch.o, ch.vars)
    push!(ch.solutions, solution)
    @show solution
    @show SCIP.sol_values(ch.o, ch.binvars)
    m, num_reactions = size(ch.S)
    complementary_vals = SCIP.sol_values(ch.o, [MOI.VariableIndex(i) for i in length(ch.internal_rxn_idxs)+num_reactions+1:2*length(ch.internal_rxn_idxs)+num_reactions])
    @show complementary_vals
    @assert SCIP.sol_values(ch.o, ch.binvars) + complementary_vals == ones(length(ch.internal_rxn_idxs))

    @show ch.S * solution
    @show ch.solutions, solution
    # @assert !(solution in ch.solutions)
    # @show solution[ch.internal_rxn_idxs]
    # @show thermo_feasible_mu(ch.internal_rxn_idxs, solution[ch.internal_rxn_idxs], ch.S)
    feasible = thermo_feasible_mu(ch.internal_rxn_idxs, solution[ch.internal_rxn_idxs], ch.S)
    @show feasible
    if !feasible      
        add_cb_cut(ch)
        # return SCIP.SCIP_CONSADDED
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
    internal_rxn_idxs = ch.internal_rxn_idxs
    S = ch.S 
    S_int = Array(S[:, internal_rxn_idxs])

    solution_master_flux = SCIP.sol_values(ch.o, ch.vars)
    solution_master_direction = solution = SCIP.sol_values(ch.o, ch.binvars)[1:length(internal_rxn_idxs)]

    @show solution_master_flux
    @show solution_master_direction

    feasible = thermo_feasible_mu(internal_rxn_idxs, solution_master_flux[internal_rxn_idxs], S)

    if !feasible
        # compute corresponding MIS
        println("_______________")
        println("compute MIS")
        m, num_reactions = size(S)
        # solution_a = solution_master[num_reactions+1:end]
        C = compute_MIS(solution_master_direction, S_int, [], internal_rxn_idxs, fast=true, time_limit=600, silent=true)
        @show C
        if isempty(C)
            feasible = thermo_feasible_mu(internal_rxn_idxs, solution_master_flux[internal_rxn_idxs], S)
            if !feasible
                @show solution_master_flux
                @show solution_master_direction
                return SCIP.SCIP_INFEASIBLE
            end
            sub_problem = Model(optimizer)
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_master_direction, internal_rxn_idxs)
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=true, time_limit=600)
            @assert termination_sub == MOI.OPTIMAL
            return SCIP.SCIP_FEASIBLE
        else 
            # build sub problem to master solution 
            sub_problem = Model(SCIP.Optimizer)
            constraint_list = build_sub_problem(sub_problem, internal_rxn_idxs, S, solution_master_direction, C)
            println("_______________")
            println("sub problem")
            # print(sub_problem)
            objective_value_sub, dual_bound_sub, solution_sub, _, termination_sub = optimize_model(sub_problem, silent=true, time_limit=600)
            add_combinatorial_benders_cut_moi(ch.o, solution_master_direction, C, ch.binvars[1:length(internal_rxn_idxs)])
            ch.ncalls += 1
            # prit
            return SCIP.SCIP_CONSADDED
        end
    end
    return SCIP.SCIP_INFEASIBLE
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

