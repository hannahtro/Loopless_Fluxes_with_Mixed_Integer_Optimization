using DataFrames, JSON, CSV

"""
compare ll-FBA, no good cuts and fast combinatorial Benders
"""
function solver_data(organisms; no_good_cuts=false, fba=false, cobrexa=false, cb=false, cb_big_m=false, yeast=false, nullspace=false, mis_numbers=[], time_limit=1800, solver="SCIP")
    if solver == "" || solver == "SCIP"
        file = string(time_limit) * ".json"
    else 
        file = solver * "_" * string(time_limit) * ".json"
    end 

    df = DataFrame(
            organism = String[], 
            time_limit= Int64[],
            termination_ll_fba = String[], 
            objective_value_ll_fba = Float64[], 
            time_ll_fba = Float64[],
            feasibility_ll_fba = Bool[]
    )

    if no_good_cuts 
        df[!, "termination_no_good_cuts"] = String[] 
        df[!, "objective_value_no_good_cuts"] = Float64[] 
        df[!, "time_no_good_cuts"] = Float64[] 
        df[!, "feasibility_no_good_cuts"] = Bool[]
        df[!, "cuts_no_good_cuts"] = Int64[]
        df[!, "iter_no_good_cuts"] = Int64[]
    end 

    if cb 
        df[!, "termination_cb"] = String[] 
        df[!, "objective_value_cb"] = Float64[] 
        df[!, "time_cb"] = Float64[] 
        df[!, "feasibility_cb"] = Bool[]
        # df[!, "cuts_cb"] = Int64[]
        df[!, "iter_cb"] = Int64[]
    end

    if cb_big_m 
        df[!, "termination_cb_big_m"] = String[] 
        df[!, "objective_value_cb_big_m"] = Float64[] 
        df[!, "time_cb_big_m"] = Float64[] 
        df[!, "feasibility_cb_big_m"] = Bool[]
        df[!, "cuts_cb_big_m"] = Int64[]
        df[!, "cuts_cb_big_m"] = Int64[]
    end

    if nullspace 
        df[!, "termination_ll_fba_nullspace"] = String[] 
        df[!, "objective_value_ll_fba_nullspace"] = Float64[] 
        df[!, "time_ll_fba_nullspace"] = Float64[]
        df[!, "feasibility_ll_fba_nullspace"] = Bool[]
    end

    if fba 
        df[!, "termination_fba"] = String[] 
        df[!, "objective_value_fba"] = Float64[] 
        df[!, "time_fba"] = Float64[]
        df[!, "feasibility_fba"] = Bool[]
    end

    if cobrexa
        df[!, "termination_ll_fba_cobrexa"] = String[] 
        df[!, "objective_value_ll_fba_cobrexa"] = Float64[] 
        df[!, "time_ll_fba_cobrexa"] = Float64[] 
        df[!, "feasibility_cobrexa"] = Bool[]
    end 

    for mis in mis_numbers
        df[!, "termination_cb_mis_" * string(mis)] = String[] 
        df[!, "objective_value_cb_mis_" * string(mis)] = Float64[] 
        df[!, "time_cb_mis_" * string(mis)] = Float64[] 
        df[!, "feasibility_cb_mis_" * string(mis)] = Bool[]
        df[!, "cuts_cb_mis_" * string(mis)] = Int64[]
        df[!, "iter_cb_mis_" * string(mis)] = Int64[]
        df[!, "times_master_problem_mis_" * string(mis)] = Float64[]
        df[!, "times_sub_problem_mis_" * string(mis)] = Float64[]
        df[!, "times_mis_problem_mis_" * string(mis)] = Float64[]
    end 
    
    # @show df

    for organism in organisms
        @show organism
        dict_organism = Dict{Symbol, Any}()
        dict_organism[:organism] = organism

        # read ll-FBA without nullspace formulation data
        dict = JSON.parse(open("json/" * organism * "_loopless_fba_" * file))
        dict_organism[:termination_ll_fba] = dict["termination"]
        dict_organism[:objective_value_ll_fba] = dict["objective_value"]
        dict_organism[:time_ll_fba] = dict["time"]
        dict_organism[:time_limit] = dict["time_limit"]
        dict_organism[:feasibility_ll_fba] = dict["thermo_feasible"]

        if nullspace
            # read ll-FBA data
            dict = JSON.parse(open("json/" * organism * "_loopless_fba_nullspace_" * file))
            dict_organism[:termination_ll_fba_nullspace] = dict["termination"]
            dict_organism[:objective_value_ll_fba_nullspace] = dict["objective_value"]
            dict_organism[:time_ll_fba_nullspace] = dict["time"]
            dict_organism[:feasibility_ll_fba_nullspace] = dict["thermo_feasible"]
        end 

        if no_good_cuts
            # read no good cuts data 
            if organism == "iAF692"
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_" * string(time_limit) * "_1e-5.json"))
            else
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_" * file))
            end

            # set termination status to TIME_LIMIT, INFEASIBLE, OPTIMAL
            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "OTHER_ERROR" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >=  dict["time_limit"]
                    dict_organism[:termination_no_good_cuts] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_no_good_cuts] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_no_good_cuts] = "OPTIMAL"
            end
            dict_organism[:objective_value_no_good_cuts] = dict["objective_value"]
            dict_organism[:time_no_good_cuts] = dict["time"]
            dict_organism[:feasibility_no_good_cuts] = dict["thermo_feasible"]
            dict_organism[:cuts_no_good_cuts] = dict["cuts"]
            dict_organism[:iter_no_good_cuts] = dict["iter"]
        end 

        if cb 
            # read fast CB data
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * file))
            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >= dict["time_limit"]
                    dict_organism[:termination_cb] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_cb] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_cb] = "OPTIMAL"
            else 
                dict_organism[:termination_cb] = dict["termination"]
            end
            dict_organism[:objective_value_cb] = dict["objective_value"]
            dict_organism[:time_cb] = dict["time"]
            dict_organism[:feasibility_cb] = dict["thermo_feasible"]
            # dict_organism[:cuts_cb] = dict["cuts"]
            dict_organism[:iter_cb] = dict["iter"]
        end 

        if cb_big_m
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * file))
            dict_organism[:termination_cb_big_m] = dict["termination"]
            dict_organism[:objective_value_cb_big_m] = dict["objective_value"]
            dict_organism[:time_cb_big_m] = dict["time"]
            dict_organism[:feasibility_cb_big_m] = dict["thermo_feasible"]
            dict_organism[:cuts_cb_big_m] = dict["cuts"]
            dict_organism[:iter_cb_big_m] = dict["iter"]
        end 

        if fba 
            dict = JSON.parse(open("json/" * organism * "_fba_" * solver * ".json"))
            dict_organism[:termination_fba] = dict["termination"]
            dict_organism[:objective_value_fba] = dict["objective_value"]
            dict_organism[:time_fba] = dict["time"]
            dict_organism[:feasibility_fba] = dict["thermo_feasible"]
        end 

        for mis in mis_numbers
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * string(mis) * "_mis_" * file))
            dict_organism[Symbol("termination_cb_mis_" * string(mis))] = dict["termination"]
            dict_organism[Symbol("objective_value_cb_mis_" * string(mis))] = dict["objective_value"]
            dict_organism[Symbol("time_cb_mis_" * string(mis))] = dict["time"]
            dict_organism[Symbol("feasibility_cb_mis_" * string(mis))] = dict["thermo_feasible"]
            dict_organism[Symbol("cuts_cb_mis_" * string(mis))] = dict["cuts"]
            dict_organism[Symbol("iter_cb_mis_" * string(mis))] = dict["iter"]
            dict_organism[Symbol("times_master_problem_mis_" * string(mis))] = mean(dict["times_master_problem"])
            dict_organism[Symbol("times_sub_problem_mis_" * string(mis))] = mean(dict["times_sub_problem"])
            dict_organism[Symbol("times_mis_problem_mis_" * string(mis))] = mean(dict["times_mis_problem"])
        end 

        for (key, value) in dict_organism
            if isnothing(value)
                dict_organism[key] = NaN
            end
        end
        # @show dict_organism
        push!(df, dict_organism)
    end
    
    # round objective value, cut time above time limit and convert to int
    for row in eachrow(df)
        if row.time_ll_fba > 1800
            row.time_ll_fba = 1800
        end

        if cobrexa
            if row.time_ll_fba_cobrexa > 1800
                row.time_ll_fba_cobrexa = 1800
            end
        end 

        if nullspace
            if row.time_ll_fba_nullspace > 1800
                row.time_ll_fba_nullspace = 1800
            end
        end
  
        if no_good_cuts
            if row.time_no_good_cuts > 1800
                row.time_no_good_cuts = 1800
            end
        end 
        if cb 
            if row.time_cb > 1800
                row.time_cb = 1800
            end
        end
    end


    df[!, :objective_value_ll_fba] = round.(df[!, :objective_value_ll_fba], digits=3)
    df[!, :time_ll_fba] = round.(df[!, :time_ll_fba], digits=0)
    df[!, :time_ll_fba] = Int.(df[!, :time_ll_fba])

    if no_good_cuts
        df[!, :objective_value_no_good_cuts] = round.(df[!, :objective_value_no_good_cuts], digits=3)
        df[!, :time_no_good_cuts] = round.(df[!, :time_no_good_cuts], digits=0)
        df[!, :time_no_good_cuts] = Int.(df[!, :time_no_good_cuts])
    end

    if cb
        df[!, :objective_value_cb] = round.(df[!, :objective_value_cb], digits=3)
        df[!, :time_cb] = round.(df[!, :time_cb], digits=0)
        df[!, :time_cb] = Int.(df[!, :time_cb])
    end

    if fba 
    end 

    # @show df[!, [:objective_value_ll_fba, :objective_value_no_good_cuts, :objective_value_cb]]
    # @show df[!, [:termination_ll_fba, :termination_no_good_cuts, :termination_cb]]
    # @show df[!, [:time_ll_fba, :time_no_good_cuts, :time_cb]]
    
    if yeast
        file_name = "results_yeast_" * solver * ".csv"
    else
        file_name = "results_bigg_" * solver * ".csv"
    end 
    CSV.write("csv/" * file_name, df, append=false, writeheader=true)
end

# organisms = ["iAF692", "e_coli_core", "iJR904", "iML1515", "iNF517", "iNJ661", "iCN900"] # "iSB619" not feasible
organisms = [
    "iAF692", # recompute for 1e-5
    "iJR904", 
    "iML1515", 
    "e_coli_core",
    "iNF517",
    # "iSB619", # AssertionError("feasible")
    "iNJ661",
    "iCN900",
    "iAF1260",
    "iEK1008",
    "iJO1366",
    "iMM904",
    "iSDY_1059",
    "iSFV_1184", # recompute on cluster
    "iSF_1195",
    "iS_1188",
    "iSbBS512_1146" # recompute on cluster
]
mis_numbers = [5, 10, 20, 30]
solver_data(organisms, time_limit=1800, yeast=false, cb=true, fba=false, cb_big_m=false, mis_numbers=mis_numbers)

# organisms = [
#     "Hanseniaspora_uvarum",
#     "yHMPu5000035696_Hanseniaspora_singularis",
#     "yHMPu5000034963_Hanseniaspora_clermontiae",
#     "yHMPu5000035695_Hanseniaspora_pseudoguilliermondii",
#     "yHMPu5000035684_Kloeckera_hatyaiensis",
#     "Eremothecium_sinecaudum",
#     "yHMPu5000035659_Saturnispora_dispora",
#     "Tortispora_caseinolytica",
#     "Starmerella_bombicola_JCM9596",
#     "Eremothecium_gossypii",
#     "Ashbya_aceri"]
# solver_data(organisms, time_limit=36000, fba=true, yeast=true, cb=true, cb_big_m=true, solver="Gurobi")