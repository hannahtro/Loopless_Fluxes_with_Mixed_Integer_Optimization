using Test, Infiltrator
using SCIP, JuMP 

include("../src/enzyme_model.jl")
include("../src/fba.jl")
include("../src/loopless_fba.jl")
include("../src/cuts_decomposition.jl")
include("../src/cobrexa.jl")

println("============================================================")
println("GECKO")
println("============================================================")

seed = 10
@testset "e_coli_core" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST ecoli core")
    println("--------------------------------------------------------")

    # build gecko model
    organism = "e_coli_core"
    model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
    gecko_model = build_gecko_model(model, seed)

    # solve fba with COBREXA
    println("COBREXA FBA")
    println("--------------------------------------------------------")
    opt_model = flux_balance_analysis(gecko_model, SCIP.Optimizer)
    sol = flux_dict(gecko_model, opt_model)
    # @show flux_summary(sol)
    objective_value_cobrexa = objective_value(opt_model)
    @show objective_value_cobrexa
    @test termination_status(opt_model) == MOI.OPTIMAL
    println("--------------------------------------------------------")

    # own fba COBREXA setup
    primal_objective_value, solution, status = cobrexa_fba_data(organism, json=false, enzyme_data=true, mute=true, seed=seed)
    @show primal_objective_value, status
    @test isapprox(objective_value_cobrexa, primal_objective_value, atol=0.001)

    # solve fba 
    println("FBA")
    println("--------------------------------------------------------")
    objective_value_fba, termination, feasible = get_fba_data(organism; json=false, enzyme_data=true, seed=seed)
    @show objective_value_fba, termination, feasible
    @test isapprox(objective_value_cobrexa, objective_value_fba, atol=0.001)
    println("--------------------------------------------------------")

    # solve ll fba
    println("ll-FBA")
    println("--------------------------------------------------------")
    objective_value_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes, termination, feasible = loopless_fba_data(organism; json=false, enzyme_data=true, seed=seed)
    @show objective_value_loopless_fba, termination, feasible
    @test feasible
    @show objective_value_loopless_fba <= objective_value_fba
    @test objective_value_loopless_fba <= objective_value_fba + 0.00001
    println("--------------------------------------------------------")

    # CB 
    println("CB")
    println("--------------------------------------------------------")
    objective_value_cb, termination, feasible = combinatorial_benders_data(organism,json=false, enzyme_data=true, seed=seed)
    @show objective_value_cb, termination, feasible
    @test feasible
    # @test isapprox(objective_value_loopless_fba, objective_value_cb, atol=0.001)
end 

# @testset "iNF517" begin
#     println("")
#     println("--------------------------------------------------------")
#     println("TEST iNF517")
#     println("--------------------------------------------------------")

#     # build gecko model
#     organism = "iNF517"
#     model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
#     gecko_model = build_gecko_model(model, seed)

#     # # solve fba with COBREXA
#     # println("COBREXA FBA")
#     # println("--------------------------------------------------------")
#     # opt_model = flux_balance_analysis(gecko_model, SCIP.Optimizer)
#     # sol = flux_dict(gecko_model, opt_model)
#     # flux_summary(sol)
#     # @test termination_status(opt_model) == MOI.OPTIMAL
#     # objective_value_fba_cobrexa = objective_value(opt_model)
#     # println("--------------------------------------------------------")

#     # solve fba 
#     println("FBA")
#     println("--------------------------------------------------------")
#     objective_value_fba, termination, feasible = get_fba_data(organism; json=false, enzyme_data=true, save_lp=false, seed=seed)
#     @show objective_value_fba, termination, feasible
#     # @show objective_value_fba, objective_value_fba_cobrexa
#     # @test isapprox(objective_value_fba, objective_value_fba_cobrexa, atol=0.001)
#     println("--------------------------------------------------------")

#     # solve ll fba
#     println("ll-FBA")
#     println("--------------------------------------------------------")
#     objective_value_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes, termination, feasible = loopless_fba_data(organism; json=false, enzyme_data=true, nullspace_formulation=false, silent=true, save_model=false, seed=seed)
#     @show objective_value_loopless_fba, termination, feasible, nodes
#     @test feasible
#     @test objective_value_loopless_fba <= objective_value_fba + 0.00001
#     println("--------------------------------------------------------")

#     # CB 
#     println("CB")
#     println("--------------------------------------------------------")
#     objective_value_cb, termination, feasible = combinatorial_benders_data(organism, json=false, enzyme_data=true, big_m=true, silent=true, save_model=true, seed=seed)
#     @show objective_value_cb, termination, feasible
#     @test feasible
#     @show objective_value_loopless_fba, objective_value_cb
#     @test isapprox(objective_value_loopless_fba, objective_value_cb, atol=0.001)
# end 

# @testset "iML1515" begin
#     println("")
#     println("--------------------------------------------------------")
#     println("TEST iML1515")
#     println("--------------------------------------------------------")

#     # build gecko model
#     organism = "iML1515"
#     model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
#     gecko_model = build_gecko_model(model, seed=10)

#     # solve fba with COBREXA
#     println("COBREXA FBA")
#     println("--------------------------------------------------------")
#     opt_model = flux_balance_analysis(gecko_model, SCIP.Optimizer)
#     sol = flux_dict(gecko_model, opt_model)
#     flux_summary(sol)
#     @test termination_status(opt_model) == MOI.OPTIMAL
#     objective_value_fba_cobrexa = objective_value(opt_model)
#     println("--------------------------------------------------------")

#     # solve fba 
#     println("FBA")
#     println("--------------------------------------------------------")
#     objective_value_fba, termination, feasible = get_fba_data(organism; json=false, enzyme_data=true)
#     @show objective_value_fba, termination, feasible
#     @test isapprox(objective_value_fba, objective_value_fba_cobrexa, atol=0.001)
#     println("--------------------------------------------------------")

#     # # solve ll fba
#     # println("ll-FBA")
#     # println("--------------------------------------------------------")
#     # objective_value_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes, termination, feasible = loopless_fba_data(organism; json=false, enzyme_data=true, nullspace_formulation=false, silent=false)
#     # @show objective_value_loopless_fba, termination, feasible, nodes
#     # @test feasible
#     # @test objective_value_loopless_fba <= objective_value + 0.00001
#     # println("--------------------------------------------------------")

#     # CB 
#     println("CB")
#     println("--------------------------------------------------------")
#     objective_value_cb, termination, feasible, x = combinatorial_benders_data(organism, json=false, enzyme_data=true, big_m=true)
#     @show objective_value_cb, termination, feasible
#     @test feasible
#     @test objective_value_cb <= objective_value_fba + 0.00001
#     # @infiltrate
#     # @test isapprox(objective_value_loopless_fba, objective_value_cb, atol=0.001)
# end 

# @testset "iND750" begin
#     println("")
#     println("--------------------------------------------------------")
#     println("TEST iND750")
#     println("--------------------------------------------------------")

#     # build gecko model
#     organism = "iND750"
#     model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
#     gecko_model = build_gecko_model(model, seed=10)

#     # solve fba with COBREXA
#     println("COBREXA FBA")
#     println("--------------------------------------------------------")
#     opt_model = flux_balance_analysis(gecko_model, SCIP.Optimizer)
#     sol = flux_dict(gecko_model, opt_model)
#     flux_summary(sol)
#     @test termination_status(opt_model) == MOI.OPTIMAL
#     objective_value_fba_cobrexa = objective_value(opt_model)
#     println("--------------------------------------------------------")

#     # solve fba 
#     println("FBA")
#     println("--------------------------------------------------------")
#     objective_value_fba, termination, feasible = get_fba_data(organism; json=false, enzyme_data=true)
#     @show objective_value_fba, termination, feasible
#     @test isapprox(objective_value_fba, objective_value_fba_cobrexa, atol=0.001)
#     println("--------------------------------------------------------")

#     # # solve ll fba
#     # println("ll-FBA")
#     # println("--------------------------------------------------------")
#     # objective_value_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes, termination, feasible = loopless_fba_data(organism; json=false, enzyme_data=true, nullspace_formulation=false, silent=false)
#     # @show objective_value_loopless_fba, termination, feasible, nodes
#     # @test feasible
#     # @test objective_value_loopless_fba <= objective_value + 0.00001
#     # println("--------------------------------------------------------")

#     # CB 
#     println("CB")
#     println("--------------------------------------------------------")
#     objective_value_cb, termination, feasible = combinatorial_benders_data(organism, json=false, enzyme_data=true)
#     @show objective_value_cb, termination, feasible
#     @test feasible
#     # @test isapprox(objective_value_loopless_fba, objective_value_cb, atol=0.001)
#     @test objective_value_cb <= objective_value_fba + 0.00001
# end 

@testset "iJR904" begin
    println("")
    println("--------------------------------------------------------")
    println("TEST iJR904")
    println("--------------------------------------------------------")

    # build gecko model
    organism = "iJR904"
    model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
    gecko_model = build_gecko_model(model, seed)

    # solve fba with COBREXA
    println("COBREXA FBA")
    println("--------------------------------------------------------")
    opt_model = flux_balance_analysis(gecko_model, SCIP.Optimizer)
    sol = flux_dict(gecko_model, opt_model)
    flux_summary(sol)
    @test termination_status(opt_model) == MOI.OPTIMAL
    objective_value_fba_cobrexa = objective_value(opt_model)
    println("--------------------------------------------------------")

    # solve fba 
    println("FBA")
    println("--------------------------------------------------------")
    objective_value_fba, termination, feasible = get_fba_data(organism; json=false, enzyme_data=true, seed=seed)
    @show objective_value_fba, termination, feasible
    @test isapprox(objective_value_fba, objective_value_fba_cobrexa, atol=0.001)
    println("--------------------------------------------------------")

    # # solve ll fba
    # println("ll-FBA")
    # println("--------------------------------------------------------")
    # objective_value_loopless_fba, vcat(x,a,G), time_loopless_fba, nodes, termination, feasible = loopless_fba_data(organism; json=false, enzyme_data=true, nullspace_formulation=false, silent=true, seed=seed)
    # @show objective_value_loopless_fba, termination, feasible, nodes
    # @test feasible
    # @test objective_value_loopless_fba <= objective_value + 0.00001
    # println("--------------------------------------------------------")

    # CB 
    # println("CB")
    # println("--------------------------------------------------------")
    # objective_value_cb, termination, feasible = combinatorial_benders_data(organism, json=false, enzyme_data=true, big_m=true, seed=seed, silent=true)
    # @show objective_value_cb, termination, feasible
    # @test feasible
    # # @test isapprox(objective_value_loopless_fba, objective_value_cb, atol=0.001)
    # @test objective_value_cb <= objective_value_fba + 0.00001
end 

