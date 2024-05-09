using DataFrames, JSON, CSV

"""
compare ll-FBA, no good cuts and fast combinatorial Benders
"""
function solver_data(organisms; no_good_cuts=false, fba=false, cobrexa=false, cb=false, cb_big_m=false, ch=false, yeast=false, nullspace=false, ch_mis=false, time_limit=1800, solver="SCIP")
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
    end 

    if cb 
        df[!, "termination_cb"] = String[] 
        df[!, "objective_value_cb"] = Float64[] 
        df[!, "time_cb"] = Float64[] 
        df[!, "feasibility_cb"] = Bool[]
    end

    if cb_big_m 
        df[!, "termination_cb_big_m"] = String[] 
        df[!, "objective_value_cb_big_m"] = Float64[] 
        df[!, "time_cb_big_m"] = Float64[] 
        df[!, "feasibility_cb_big_m"] = Bool[]
    end

    if ch 
        df[!, "termination_ch"] = Union{String, Missing}[]
        df[!, "objective_value_ch"] = Union{Float64, Missing}[]
        df[!, "time_ch"] = Union{Float64, Missing}[] 
        df[!, "feasibility_ch"] = Union{Bool, Missing}[]
    end

    if ch_mis
        df[!, "termination_ch_mis_5"] = Union{String, Missing}[]
        df[!, "objective_value_ch_mis_5"] = Union{Float64, Missing}[]
        df[!, "time_ch_mis_5"] = Union{Float64, Missing}[] 
        df[!, "feasibility_ch_mis_5"] = Union{Bool, Missing}[]
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
        end 

        if cb_big_m
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * file))
            dict_organism[:termination_cb_big_m] = dict["termination"]
            dict_organism[:objective_value_cb_big_m] = dict["objective_value"]
            dict_organism[:time_cb_big_m] = dict["time"]
            dict_organism[:feasibility_cb_big_m] = dict["thermo_feasible"]
        end 

        if ch
            try 
                dict = JSON.parse(open("json/" * organism * "_constraint_handler_" * file))
                if dict["thermo_feasible"] == NaN || isnothing(dict["thermo_feasible"])
                    dict["thermo_feasible"] = missing
                end

            catch e         
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => missing,
                    "time" => missing, 
                    "thermo_feasible" => missing,
                )
            end  
            dict_organism[:termination_ch] = dict["termination"]
            dict_organism[:objective_value_ch] = dict["objective_value"]
            dict_organism[:time_ch] = dict["time"]
            dict_organism[:feasibility_ch] = dict["thermo_feasible"]
        end 

        if ch_mis
            try 
                dict = JSON.parse(open("json/" * organism * "_constraint_handler_5_mis_" * file))
                if dict["thermo_feasible"] == NaN || isnothing(dict["thermo_feasible"])
                    dict["thermo_feasible"] = missing
                end
            catch e        
                println(e) 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => missing,
                    "time" => missing, 
                    "thermo_feasible" => missing,
                )
            end  
            dict_organism[:termination_ch_mis_5] = dict["termination"]
            dict_organism[:objective_value_ch_mis_5] = dict["objective_value"]
            dict_organism[:time_ch_mis_5] = dict["time"]
            dict_organism[:feasibility_ch_mis_5] = dict["thermo_feasible"]
        end 

        if fba 
            dict = JSON.parse(open("json/" * organism * "_fba_" * solver * ".json"))
            dict_organism[:termination_fba] = dict["termination"]
            dict_organism[:objective_value_fba] = dict["objective_value"]
            dict_organism[:time_fba] = dict["time"]
            dict_organism[:feasibility_fba] = dict["thermo_feasible"]
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
    column_names = names(df)

    for name in column_names 
        if occursin("objective_value", name)
            @show name
            df[!, name] = round.(df[!, name], digits=3)
        end 
    end

    for name in column_names 
        if occursin("time", name) && !occursin("time_limit", name)
            @show name
            df[!, name] = [((item >= 1800) && !isnan(item)) ? 1800 : item for item in df[!, name]]
            df[!, name] = [isnan(item) ? item : Int(round(item, digits=0)) for item in df[!, name]]
            df[!, name] = [isnan(item) ? missing : item for item in df[!, name]]
        end 
    end
    
    if yeast
        type = "results_yeast_ch"
    else
        type = "results_bigg_ch"
    end 
    if solver != "SCIP"
        type = type * "_" * solver
    end
    file_name = type * ".csv"
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
solver_data(organisms, time_limit=1800, yeast=false, cb=false, fba=false, cb_big_m=false, ch=true, ch_mis=true, solver="SCIP")

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