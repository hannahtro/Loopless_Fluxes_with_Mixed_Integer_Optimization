using COBREXA
using JuMP 
using SCIP, Gurobi
using CSV, DataFrames

include("../src/enzyme_model.jl")
include("../src/optimization_model.jl")
include("../src/cuts_decomposition.jl")

# ME to reproduce AssertionError("objective_value_master <= objective_value_fba + 0.0001")
function reproduce_error(organism, seed)
    println("")
    @show organism, seed 

    # load data and setup
    println("SCIP OPTIMIZER") 
    optimizer = SCIP.Optimizer
    time_limit = 600
    silent = true

    molecular_model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
    molecular_model = build_gecko_model(molecular_model, seed)

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model, enzyme_data=true)

    # build dict
    dict = Dict{Symbol, Any}()
    dict[:organism] = organism
    dict[:seed] = seed
    solution_dict = Dict{Symbol, Any}()

    # solve FBA 
    println("FBA")
    master_problem = make_optimization_model(molecular_model, optimizer)
    @show objective_sense(master_problem)
    MOI.set(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"), 1e-8)
    @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_fba, _, _, _, _ = optimize_model(master_problem, silent=silent)
    @show objective_value_fba
    dict[:objective_value_fba] = objective_value_fba
    println("")
    solution_dict[:objective_value_scip] = objective_value_fba
    solution_dict[:x_scip] = [value(var) for var in all_variables(master_problem)]

    # solve master problem with indicator variables
    println("MP with indicator constraints")
    master_problem = make_optimization_model(molecular_model, optimizer)
    build_master_problem(master_problem, internal_rxn_idxs)   
    @show objective_sense(master_problem)
    MOI.set(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"), 1e-8)
    @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    @show objective_value_master, objective_value_fba
    if !isapprox(objective_value_master, objective_value_fba, atol=0.0001)
        @warn "objective values do not match"
    end 
    dict[:objective_value_mp_indicator] = objective_value_master
    println("") 

    # solve master problem with big M 
    println("MP with big M")
    master_problem = make_optimization_model(molecular_model, optimizer)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))
    build_master_problem(master_problem, internal_rxn_idxs, max_flux_bound, big_m=true)   
    @show objective_sense(master_problem)
    MOI.set(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"), 1e-8)
    @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    @show objective_value_master, objective_value_fba
    if !isapprox(objective_value_master, objective_value_fba, atol=0.0001)
        @warn "objective values do not match"
    end 
    dict[:objective_value_mp_bigm] = objective_value_master
    println("") 

    # GUROBI
    println("Gurobi optimizer")
    optimizer = Gurobi.Optimizer
    # solve FBA 
    master_problem = make_optimization_model(molecular_model, optimizer)
    @show objective_sense(master_problem)
    set_optimizer_attribute(master_problem, "FeasibilityTol", 1e-6)
    set_optimizer_attribute(master_problem, "OptimalityTol", 1e-6)
    # @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_fba, _, _, _, _ = optimize_model(master_problem)
    @show objective_value_fba 
    println("")
    solution_dict[:objective_value_gurobi] = objective_value_fba
    solution_dict[:x_gurobi] = [value(var) for var in all_variables(master_problem)]

    # solve master problem with indicator variables
    println("MP with indicator constraints")
    master_problem = make_optimization_model(molecular_model, optimizer)
    build_master_problem(master_problem, internal_rxn_idxs)   
    @show objective_sense(master_problem)
    set_optimizer_attribute(master_problem, "FeasibilityTol", 1e-6)
    set_optimizer_attribute(master_problem, "OptimalityTol", 1e-6)
    # @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    @show objective_value_master, objective_value_fba
    if !isapprox(objective_value_master, objective_value_fba, atol=0.0001)
        @warn "objective values do not match"
    end 
    println("")

    # solve master problem with big M 
    println("MP with big M")    
    master_problem = make_optimization_model(molecular_model, optimizer)
    max_flux_bound = maximum(abs.(vcat(lb, ub)))
    build_master_problem(master_problem, internal_rxn_idxs, max_flux_bound, big_m=true)   
    @show objective_sense(master_problem)
    set_optimizer_attribute(master_problem, "FeasibilityTol", 1e-6)
    set_optimizer_attribute(master_problem, "OptimalityTol", 1e-6)
    # @show MOI.get(master_problem, MOI.RawOptimizerAttribute("numerics/feastol"))
    objective_value_master, dual_bound_master, solution_master, _, termination_master = optimize_model(master_problem, time_limit=time_limit, silent=silent)
    @show objective_value_master, objective_value_fba
    if !isapprox(objective_value_master, objective_value_fba, atol=0.0001)
        @warn "objective values do not match"
    end     

    # open(organism * "_scip_vs_gurobi.json", "w") do f
    #     JSON.print(f, solution_dict) 
    # end

    return dict
end

function filter_csv(file_name; bigm=false, indicator=false, nonzero=false)
    data, header = readdlm(file_name, ',', header=true);
    df = DataFrame(data, vec(header))

    df[!,:objective_value_fba] = round.(df[:,:objective_value_fba], digits=6)
    df[!,:objective_value_mp_indicator] = round.(df[:,:objective_value_mp_indicator], digits=6)
    df[!,:objective_value_mp_bigm] = round.(df[:,:objective_value_mp_bigm], digits=6)

    if nonzero 
        df = filter(:objective_value_fba => !=(0), df)
    end 

    if indicator
        function approx_indicator_filter(objective_value_fba, objective_value_mp_indicator)::Bool
            isapprox(objective_value_mp_indicator, objective_value_fba, atol=0.0001)
        end
        @show filter([:objective_value_fba, :objective_value_mp_indicator] => !approx_indicator_filter, df)

        df = filter([:objective_value_fba, :objective_value_mp_indicator] => !approx_indicator_filter, df)
        sort!(df, [:organism])

        file_name = "csv/assert_error_enzyme_models_indicator.csv"
        CSV.write(file_name, df, append=false, writeheader=true)
    end 

    if bigm
        function approx_bigm_filter(objective_value_fba, objective_value_mp_bigm)::Bool
            isapprox(objective_value_mp_bigm, objective_value_fba, atol=0.0001)
        end
        @show filter([:objective_value_fba, :objective_value_mp_bigm] => !approx_bigm_filter, df)

        df = filter([:objective_value_fba, :objective_value_mp_bigm] => approx_bigm_filter, df)
        sort!(df, [:organism])

        file_name = "csv/assert_error_enzyme_models_bigm.csv"
        CSV.write(file_name, df, append=false, writeheader=true)
    end 
    return df
end 

function check_feasibility(organism, seed)
    println("SCIP OPTIMIZER") 
    optimizer = SCIP.Optimizer
    time_limit = 600
    silent = false

    molecular_model = load_model(StandardModel, "../molecular_models/" * organism * ".json")
    molecular_model = build_gecko_model(molecular_model, seed)

    S = stoichiometry(molecular_model)
    m, num_reactions = size(S)
    lb, ub = bounds(molecular_model)
    internal_rxn_idxs = internal_reactions(molecular_model, enzyme_data=true)

    master_problem = make_optimization_model(molecular_model, optimizer)
    @show objective_sense(master_problem)

    dict = JSON.parse(open(organism * "_scip_vs_gurobi.json"))
    x = master_problem[:x]
    point = Dict(x .=> dict["x_gurobi"])
    report = primal_feasibility_report(master_problem, point, atol=0.0000001)
    file_name = organism * "_report_scip_gurobi.txt"
    writedlm(file_name, report)
end 

### get objective values for different models and seeds
# organisms = [
#     "iAF692", 
#     "iJR904", 
#     "iML1515", 
#     "e_coli_core",
#     "iNF517",
#     "iSB619",
#     # "iNJ661",
#     "iCN900",
#     "iAF1260",
#     # "iEK1008",
#     "iJO1366",
#     "iMM904",
#     "iSDY_1059",
#     "iSFV_1184",
#     "iSF_1195",
#     "iS_1188",
#     "iSbBS512_1146"
# ]

# seeds = [1, 2, 3]

# df = DataFrame(
#     organism = String[], 
#     seed = Int64[],
#     objective_value_fba = Float64[],
#     objective_value_mp_indicator = Float64[],
#     objective_value_mp_bigm = Float64[]
# )

# for seed in seeds
#     for organism in organisms
#         dict = reproduce_error(organism, seed)
#         push!(df, dict)
#     end 
# end


### build filtered df
# file_name = "assert_error_enzyme_models.csv"
# CSV.write("csv/" * file_name, df, append=false, writeheader=true)

# file_name = "csv/assert_error_enzyme_models.csv"
# df = filter_csv(file_name, nonzero=false, indicator=false, bigm=true)


### store solutions
organism = "iAF692"
seed = 2
reproduce_error(organism, seed)

organism = "iAF1260"
seed = 1
reproduce_error(organism, seed)




# check_feasibility("iAF1260", 1)
# check_feasibility("iAF692", 2)
