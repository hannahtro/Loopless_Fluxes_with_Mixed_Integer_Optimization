using Test 
using DataFrames
using CSV
using Infiltrator

include("../src/cuts_decomposition.jl")
include("../src/constraint_handler.jl")

@testset "simple model" begin
    println("--------------------------------------------------------")
    println("TEST SIMPLE MODEL")
    println("--------------------------------------------------------")
    S = [[1,0,0,0] [-1,1,0,0] [0,-1,1,0] [-1,0,1,0] [0,0,-1,0] [-1,0,0,1] [0,0,1,-1]]
    lb = [0,-10,-10,-10,0,0,0]
    ub = [20,30,30,30,20,10,10]
    m, num_reactions = size(S)
    @show m, num_reactions

    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    internal_rxn_idxs = [2,3,4,6,7]
    println("--------------------------------------------------------")

    println("no good cuts")
    # no good cuts
    objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S)
    solution_flux = solution[1:num_reactions]
    @show solution_flux
    solution_direction = solution[num_reactions+1:end]
    @show solution_direction
    nonzero_flux_idxs = [idx for (idx,i) in enumerate(solution_flux[internal_rxn_idxs]) if !isapprox(0, i, atol=0.001)]
    @show nonzero_flux_idxs
    @test thermo_feasible_mu(internal_rxn_idxs[nonzero_flux_idxs], solution_direction[nonzero_flux_idxs], S)
    println("--------------------------------------------------------")

    println("combinatorial Benders")
    # combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    objective_value_cb, objective_values_cb, dual_bounds_cb, solution_cb, time_cb, termination_cb, iter_cb = combinatorial_benders(model, internal_rxn_idxs, S, fast=false)
    @show objective_value_cb, solution_cb
    @test termination_cb == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_cb[internal_rxn_idxs], S)
    @test feasible
    @show objective_value, objective_value_cb
    @test round.(objective_values_cb, digits=5) == round.(sort(objective_values_cb, rev=true), digits=5)
    # @test isapprox(objective_value, objective_value_cb, atol=0.0001)

    println("--------------------------------------------------------")

    println("fast combinatorial Benders")
    # fast combinatorial Benders'
    model = build_fba_model(S, lb, ub, set_objective=true)
    # print(model)
    objective_value_fast, objective_values_fast, dual_bounds_fast, solution_fast, time_fast, termination_fast, iter_fast = combinatorial_benders(model, internal_rxn_idxs, S, fast=true)
    @test termination_fast == MOI.OPTIMAL 
    feasible = thermo_feasible(internal_rxn_idxs, solution_fast[internal_rxn_idxs], S)
    @test feasible
    @test round.(objective_values_fast, digits=5) == round.(sort(objective_values_fast, rev=true), digits=5)

    @test iter_cb >= iter_fast
    # @test time >= time_fast
    @test isapprox(objective_value_cb, objective_value_fast)
    @test solution_cb[1:num_reactions] == solution_fast[1:num_reactions]
    @show objective_value_fast, solution_fast
    println("--------------------------------------------------------")

    println("constraint handler")
    # test constraint handler    
    scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, silent=true)
    # print(scip_model)

    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [])
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch", enforce_priority=-99999999, check_priority=-99999999)

    MOI.optimize!(scip_model)
    @test MOI.get(scip_model, MOI.TerminationStatus()) == MOI.OPTIMAL
    objective_value_ch = MOI.get(scip_model, MOI.ObjectiveValue())
    solution_ch = [MOI.get(ch.o, MOI.VariablePrimal(1), MOI.VariableIndex(i)) for i in 1:length(internal_rxn_idxs) + num_reactions]
    @show objective_value_ch
    @show solution_ch 
    feasible = thermo_feasible(internal_rxn_idxs, solution_ch[internal_rxn_idxs], S)
    @test feasible
    @test isapprox(objective_value_ch,objective_value_fast)
    
    # print SCIP solution
    # SCIP.SCIPprintSol(ch.o, SCIP.SCIPgetBestSol(ch.o), C_NULL, SCIP.TRUE)
end

# TODO: no good cuts approach does not terminate in 200 iterations: verify that solution is eventually found
@testset "iAF692" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST iAF692")
    println("--------------------------------------------------------")
    organism = "iAF692"
    molecular_model = deserialize("../data/" * organism * ".js")
    # print_model(molecular_model, "organism")

    S = stoichiometry(molecular_model)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = [
        ridx for (ridx, rid) in enumerate(variables(molecular_model)) if
        !is_boundary(reaction_stoichiometry(molecular_model, rid))
    ]

    # model = build_fba_model(S, lb, ub, optimizer=SCIP.Optimizer)
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    time_limit = 2
    objective_value, dual_bound, solution, time, termination, iter = no_good_cuts(model, internal_rxn_idxs, S, time_limit=time_limit)

    try 
        thermo_feasible_mu(internal_rxn_idxs,solution[internal_rxn_idxs], S)
    catch 
    else 
        @test time >= time_limit
    end

    println("combinatorial Benders")
    # combinatorial Benders'
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    combinatorial_benders(model, internal_rxn_idxs, S, max_iter=5, fast=false)
    println("--------------------------------------------------------")

    println("fast combinatorial Benders")
    # fast combinatorial Benders'
    model = make_optimization_model(molecular_model, SCIP.Optimizer)
    combinatorial_benders(model, internal_rxn_idxs, S, max_iter=5, fast=true)
    println("--------------------------------------------------------")
    
    println("constraint handler")
    # test constraint handler
    # extract objective
    objective_func = objective_function(model)
    objective_func_vars = [i.index for i in objective_func.terms.keys]
           
    scip_model, bin_vars, flux_vars = build_fba_indicator_model_moi(S, lb, ub, internal_rxn_idxs, set_objective=true, time_limit=10, objective_func_vars=objective_func_vars, objective_func_coeffs=objective_func.terms.vals)
    # print(scip_model)
    ch = ThermoFeasibleConstaintHandler(scip_model, 0, internal_rxn_idxs, S, flux_vars, bin_vars, [], [], [])
    SCIP.include_conshdlr(scip_model, ch; needs_constraints=false, name="thermodynamically_feasible_ch", enforce_priority=-7000000, check_priority=-7000000)
    MOI.optimize!(scip_model)
    # @test MOI.get(scip_model, MOI.TerminationStatus()) == MOI.TIME_LIMIT
end

constraint_handler_data("iAF692", time_limit=600)
# no_good_cuts_data("iAF692", time_limit=3600)

println("--------------------------------------------------------")
# combinatorial_benders_data("iAF692", time_limit=1800, csv=true, fast=false, silent=false)
# combinatorial_benders_data("iAF692", time_limit=1800, csv=true, fast=true, silent=false)
# combinatorial_benders_data("iJR904", time_limit=1800, csv=true, fast=false, silent=false)
# #combinatorial_benders_data("iJR904", time_limit=1800, csv=true, fast=true, silent=false)
# combinatorial_benders_data("iML1515", time_limit=1800, csv=true, fast=false, silent=false)
# #combinatorial_benders_data("iML1515", time_limit=1800, csv=true, fast=true, silent=false)

# constraint_handler_data("iAF692", csv=true, silent=false)
# combinatorial_benders_data("iAF692", time_limit=1800, csv=true, fast=true, silent=true)


organisms = ["iAF692", "iJR904", "iML1515", "e_coli_core", "iNF517", "iSB619", "iNJ661", "iCN900"]

for organism in organisms
    type = "no_good_cut"
    try 
        no_good_cuts_data(organism, time_limit=600)
    catch e 
        println(e)
        file = organism * "_" * type
        open(file * ".txt","a") do io
            println(io, e)
        end
    end

    type = "cb"
    try 
        combinatorial_benders_data(organism, time_limit=600, fast=false)
    catch e 
        println(e)
        file = organism * "_" * type
        open(file * ".txt","a") do io
            println(io, e)
        end
    end

    type = "cb_fast"
    try 
        combinatorial_benders_data(organism, time_limit=600, fast=true)
    catch e 
        println(e)
        file = organism * "_" * type
        open(file * ".txt","a") do io
            println(io, e)
        end
    end
    
    type = "ch"
    try 
        constraint_handler_data(organism, time_limit=600)
    catch e 
        println(e)
        file = organism * "_" * type
        open(file * ".txt","a") do io
            println(io, e)
        end
    end
end

