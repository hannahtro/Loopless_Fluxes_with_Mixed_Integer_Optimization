using DataFrames, JSON, CSV

"""
compare ll-FBA, no good cuts and fast combinatorial Benders
"""
function loopless_fba_vs_cb(organisms; cuts=false, ll_formulations=false, yeast=false, time_limit=1800)
    if yeast 
        df = DataFrame(
            organism = String[], 
            time_limit= Int64[],
            termination_ll_fba = String[], 
            objective_value_ll_fba = Float64[], 
            time_ll_fba = Float64[], 
            termination_no_good_cuts = String[], 
            objective_value_no_good_cuts = Float64[], 
            time_no_good_cuts = Float64[], 
            termination_cb = String[], 
            objective_value_cb = Float64[], 
            time_cb = Float64[], 
        )
    else 
        df = DataFrame(
            organism = String[], 
            time_limit= Int64[],
            termination_ll_fba_cobrexa = String[], 
            objective_value_ll_fba_cobrexa = Float64[], 
            time_ll_fba_cobrexa = Float64[], 
            termination_ll_fba_nullspace = String[], 
            objective_value_ll_fba_nullspace = Float64[], 
            time_ll_fba_nullspace = Float64[], 
            termination_ll_fba = String[], 
            objective_value_ll_fba = Float64[], 
            time_ll_fba = Float64[], 
            termination_no_good_cuts = String[], 
            objective_value_no_good_cuts = Float64[], 
            time_no_good_cuts = Float64[], 
            termination_cb = String[], 
            objective_value_cb = Float64[], 
            time_cb = Float64[], 
        )
    end
    
    for organism in organisms
        @show organism
        dict_organism = Dict{Symbol, Any}()
        dict_organism[:organism] = organism

        # read cobrexa data
        if !yeast
            dict = JSON.parse(open("json/" * organism * "_cobrexa_loopless_fba_" * string(time_limit) * ".json"))
            dict_organism[:termination_ll_fba_cobrexa] = dict["termination"]
            dict_organism[:objective_value_ll_fba_cobrexa] = dict["objective_value"]
            dict_organism[:time_ll_fba_cobrexa] = dict["time"]
        end

        # read ll-FBA data
        if !yeast
            dict = JSON.parse(open("json/" * organism * "_loopless_fba_nullspace_" * string(time_limit) * ".json"))
            dict_organism[:termination_ll_fba_nullspace] = dict["termination"]
            dict_organism[:objective_value_ll_fba_nullspace] = dict["objective_value"]
            dict_organism[:time_ll_fba_nullspace] = dict["time"]
        end

        # read ll-FBA without nullspace formulation data
        dict = JSON.parse(open("json/" * organism * "_loopless_fba_" * string(time_limit) * ".json"))
        dict_organism[:termination_ll_fba] = dict["termination"]
        dict_organism[:objective_value_ll_fba] = dict["objective_value"]
        dict_organism[:time_ll_fba] = dict["time"]
        dict_organism[:time_limit] = dict["time_limit"]

        # read no good cuts data 
        if organism == "iAF692"
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_" * string(time_limit) * "_1e-5.json"))
        else
            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_" * string(time_limit) * ".json"))
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

        # read fast CB data
        dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * string(time_limit) * ".json"))
        if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
            if dict["time"] >= dict["time_limit"]
                dict_organism[:termination_cb] = "TIME_LIMIT"
            else 
                dict_organism[:termination_cb] = "INFEASIBLE"
            end
        elseif dict["thermo_feasible"] == true
            dict_organism[:termination_cb] = "OPTIMAL"
        end
        dict_organism[:objective_value_cb] = dict["objective_value"]
        dict_organism[:time_cb] = dict["time"]

        for (key, value) in dict_organism
            if isnothing(value)
                dict_organism[key] = NaN
            end
        end
        @show dict_organism
        push!(df, dict_organism)
    end
    
    # round objective value, cut time above time limit and convert to int
    for row in eachrow(df)
        if !yeast 
            if row.time_ll_fba_cobrexa > 1800
                row.time_ll_fba_cobrexa = 1800
            end
            if row.time_ll_fba_nullspace > 1800
                row.time_ll_fba_nullspace = 1800
            end
        end
        if row.time_ll_fba > 1800
            row.time_ll_fba = 1800
        end
        if row.time_no_good_cuts > 1800
            row.time_no_good_cuts = 1800
        end
        if row.time_cb > 1800
            row.time_cb = 1800
        end
    end

    if !yeast 
        df[!, :objective_value_ll_fba_cobrexa] = round.(df[!, :objective_value_ll_fba_cobrexa], digits=3)
        df[!, :objective_value_ll_fba_nullspace] = round.(df[!, :objective_value_ll_fba_nullspace], digits=3)
        df[!, :time_ll_fba_cobrexa] = round.(df[!, :time_ll_fba_cobrexa], digits=0)
        df[!, :time_ll_fba_cobrexa] = Int.(df[!, :time_ll_fba_cobrexa])
        df[!, :time_ll_fba_nullspace] = round.(df[!, :time_ll_fba_nullspace], digits=0)
        df[!, :time_ll_fba_nullspace] = Int.(df[!, :time_ll_fba_nullspace])
    end 
    df[!, :objective_value_ll_fba] = round.(df[!, :objective_value_ll_fba], digits=3)
    df[!, :objective_value_no_good_cuts] = round.(df[!, :objective_value_no_good_cuts], digits=3)
    df[!, :objective_value_cb] = round.(df[!, :objective_value_cb], digits=3)
    df[!, :time_ll_fba] = round.(df[!, :time_ll_fba], digits=0)
    df[!, :time_ll_fba] = Int.(df[!, :time_ll_fba])
    df[!, :time_no_good_cuts] = round.(df[!, :time_no_good_cuts], digits=0)
    df[!, :time_no_good_cuts] = Int.(df[!, :time_no_good_cuts])
    df[!, :time_cb] = round.(df[!, :time_cb], digits=0)
    df[!, :time_cb] = Int.(df[!, :time_cb])

    @show df[!, [:objective_value_ll_fba, :objective_value_no_good_cuts, :objective_value_cb]]
    @show df[!, [:termination_ll_fba, :termination_no_good_cuts, :termination_cb]]
    @show df[!, [:time_ll_fba, :time_no_good_cuts, :time_cb]]
    
    if cuts 
        df = df[!, [:organism, :time_limit, :termination_no_good_cuts, :objective_value_no_good_cuts, :time_no_good_cuts, :termination_cb, :objective_value_cb, :time_cb]]
        if yeast 
            file_name = "comparison_ll_fba_vs_cb_yeast.csv"
        else 
            file_name = "comparison_ll_fba_vs_cb.csv"
        end
        CSV.write("csv/" * file_name, df, append=false, writeheader=true)
    elseif ll_formulations
        df = df[!, [:organism, :time_limit,  :termination_ll_fba_nullspace, :objective_value_ll_fba_nullspace, :time_ll_fba_nullspace, :termination_ll_fba, :objective_value_ll_fba, :time_ll_fba]]
        if yeast
            file_name = "comparison_ll_formulations_yeast.csv"
        else
            file_name = "comparison_ll_formulations.csv"
        end 
        CSV.write("csv/" * file_name, df, append=false, writeheader=true)
    else 
        if yeast
            file_name = "comparison_ll_fba_yeast.csv"
        else
            file_name = "comparison_ll_fba.csv"
        end 
        CSV.write("csv/" * file_name, df, append=false, writeheader=true)
    end
end

# organisms = ["iAF692", "e_coli_core", "iJR904", "iML1515", "iNF517", "iNJ661", "iCN900"] # "iSB619" not feasible
organisms = [
    # "iAF692", # recompute for 1e-5
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
loopless_fba_vs_cb(organisms, cuts=true, yeast=false, time_limit=1800)

# organisms = [
#     "Alloascoidea_hylecoeti",
#     "Ambrosiozyma_kashinagacola",
#     "Ambrosiozyma_monospora",
#     # "Arthrobotrys_oligospora", # rerun cb fast 
#     "Arxula_adeninivorans",
#     # "Ascoidea_asiatica", # rerun cb
#     # "Ascoidea_rubescens", # rerun cb
#     # "Ashbya_aceri", # rerun cb
#     # "Aspergillus_nidulans", # rerun cb
#     # "Babjeviella_inositovora", # rerun cb
#     #"Botrytis_cinerea" # rerun cb fast
# ]
# loopless_fba_vs_cb(organisms, cuts=true, yeast=true, time_limit=7200)
