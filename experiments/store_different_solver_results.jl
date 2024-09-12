using DataFrames, JSON, CSV
using Statistics
using StatsBase

"""
compare ll-FBA, no good cuts and fast combinatorial Benders
"""
function solver_data(organisms; no_good_cuts=false, no_good_cuts_big_m=false, fba=false, cobrexa=false, cb=false, cb_big_m=false, cb_indicator_and_big_m=false, yeast=false, nullspace=false, mis_indicator=false, mis_big_m=false, ll_fba_indicator=false, mis_numbers=[], cut_densities=[], distinct_cuts=false, time_limit=1800, solver="SCIP")
    @show mis_numbers
    if solver == "" || solver == "SCIP"
        file = string(time_limit) * ".json"
    else 
        file = solver * "_" * string(time_limit) * ".json"
    end 

    df = DataFrame(
            organism = Union{String,Missing}[], 
            time_limit= Union{Int64,Missing}[],
            termination_ll_fba = Union{String,Missing}[], 
            objective_value_ll_fba = Union{Float64,Missing}[], 
            time_ll_fba = Union{Float64,Missing}[],
            feasibility_ll_fba = Union{Bool,Missing}[]
    )

    if no_good_cuts 
        df[!, "termination_no_good_cuts"] = Union{String,Missing}[] 
        df[!, "objective_value_no_good_cuts"] = Union{Float64,Missing}[] 
        df[!, "time_no_good_cuts"] = Union{Float64,Missing}[] 
        df[!, "feasibility_no_good_cuts"] = Union{Bool,Missing}[]
        df[!, "cuts_no_good_cuts"] = Union{Int64,Missing}[]
        df[!, "iter_no_good_iter"] = Union{Int64,Missing}[]
    end 

    if no_good_cuts_big_m 
        df[!, "termination_no_good_cuts_big_m"] = Union{String,Missing}[] 
        df[!, "objective_value_no_good_cuts_big_m"] = Union{Float64,Missing}[] 
        df[!, "time_no_good_cuts_big_m"] = Union{Float64,Missing}[] 
        df[!, "feasibility_no_good_cuts_big_m"] = Union{Bool,Missing}[]
        df[!, "cuts_no_good_cuts_big_m"] = Union{Int64,Missing}[]
        df[!, "iter_no_good_iter_big_m"] = Union{Int64,Missing}[]
    end 

    if cb 
        df[!, "termination_cb"] = Union{String,Missing}[] 
        df[!, "objective_value_cb"] = Union{Float64,Missing}[] 
        df[!, "time_cb"] = Union{Float64,Missing}[] 
        df[!, "feasibility_cb"] = Union{Bool,Missing}[]
        df[!, "cuts_cb"] = Union{Int64,Missing}[]
        df[!, "iter_cb"] = Union{Int64,Missing}[]
        df[!, "times_master_problem_cb"] = Union{Float64,Missing}[]
        df[!, "times_sub_problem_cb"] = Union{Float64,Missing}[]
        df[!, "times_mis_problem_cb"] = Union{Float64,Missing}[]
    end

    if cb_big_m 
        df[!, "termination_cb_big_m"] = Union{String,Missing}[] 
        df[!, "objective_value_cb_big_m"] = Union{Float64,Missing}[] 
        df[!, "time_cb_big_m"] = Union{Float64,Missing}[] 
        df[!, "feasibility_cb_big_m"] = Union{Bool,Missing}[]
        df[!, "cuts_cb_big_m"] = Union{Int64,Missing}[]
        df[!, "iter_cb_big_m"] = Union{Int64,Missing}[]
        df[!, "times_master_problem_cb_big_m"] = Union{Float64,Missing}[]
        df[!, "times_sub_problem_cb_big_m"] = Union{Float64,Missing}[]
        df[!, "times_mis_problem_cb_big_m"] = Union{Float64,Missing}[]
    end

    if cb_indicator_and_big_m
        df[!, "termination_cb_indicator_and_big_m"] = Union{String,Missing}[] 
        df[!, "objective_value_cb_indicator_and_big_m"] = Union{Float64,Missing}[] 
        df[!, "time_cb_indicator_and_big_m"] = Union{Float64,Missing}[] 
        df[!, "feasibility_cb_indicator_and_big_m"] = Union{Bool,Missing}[]
        df[!, "cuts_cb_indicator_and_big_m"] = Union{Int64,Missing}[]
        df[!, "iter_cb_indicator_and_big_m"] = Union{Int64,Missing}[]
        df[!, "times_master_problem_cb_indicator_and_big_m"] = Union{Float64,Missing}[]
        df[!, "times_sub_problem_cb_indicator_and_big_m"] = Union{Float64,Missing}[]
        df[!, "times_mis_problem_cb_indicator_and_big_m"] = Union{Float64,Missing}[]
    end

    if nullspace 
        df[!, "termination_ll_fba_nullspace"] = Union{String,Missing}[] 
        df[!, "objective_value_ll_fba_nullspace"] = Union{Float64,Missing}[] 
        df[!, "time_ll_fba_nullspace"] = Union{Float64,Missing}[]
        df[!, "feasibility_ll_fba_nullspace"] = Union{Bool,Missing}[]
    end

    if ll_fba_indicator 
        df[!, "termination_ll_fba_indicator"] = Union{String,Missing}[] 
        df[!, "objective_value_ll_fba_indicator"] = Union{Float64,Missing}[] 
        df[!, "time_ll_fba_indicator"] = Union{Float64,Missing}[]
        df[!, "feasibility_ll_fba_indicator"] = Union{Bool,Missing}[]
    end

    if fba 
        df[!, "termination_fba"] = Union{String,Missing}[] 
        df[!, "objective_value_fba"] = Union{Float64,Missing}[] 
        df[!, "time_fba"] = Union{Float64,Missing}[]
        df[!, "feasibility_fba"] = Union{Bool,Missing}[]
    end

    if cobrexa
        df[!, "termination_ll_fba_cobrexa"] = Union{String,Missing}[] 
        df[!, "objective_value_ll_fba_cobrexa"] = Union{Float64,Missing}[] 
        df[!, "time_ll_fba_cobrexa"] = Union{Float64,Missing}[] 
        df[!, "feasibility_cobrexa"] = Union{Bool,Missing}[]
    end 

    for mis in mis_numbers
        if mis_indicator
            df[!, "termination_cb_mis_" * string(mis)] = Union{String,Missing}[] 
            df[!, "objective_value_cb_mis_" * string(mis)] = Union{Float64,Missing}[] 
            df[!, "time_cb_mis_" * string(mis)] = Union{Float64,Missing}[] 
            df[!, "feasibility_cb_mis_" * string(mis)] = Union{Bool,Missing}[]
            df[!, "cuts_cb_mis_" * string(mis)] = Union{Int64,Missing}[]
            df[!, "iter_cb_mis_" * string(mis)] = Union{Int64,Missing}[]
            df[!, "times_master_problem_mis_" * string(mis)] = Union{Float64,Missing}[]
            df[!, "times_sub_problem_mis_" * string(mis)] = Union{Float64,Missing}[]
            df[!, "times_mis_problem_mis_" * string(mis)] = Union{Float64,Missing}[]
            if distinct_cuts
                df[!, "termination_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{String,Missing}[] 
                df[!, "objective_value_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[] 
                df[!, "time_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[] 
                df[!, "feasibility_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{Bool,Missing}[]
                df[!, "cuts_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{Int64,Missing}[]
                df[!, "iter_cb_mis_" * string(mis) * "_distinct_cuts"] = Union{Int64,Missing}[]
                df[!, "times_master_problem_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                df[!, "times_sub_problem_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                df[!, "times_mis_problem_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                # df[!, "times_filtering_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
            end 
            if !isempty(cut_densities)
                for density in cut_densities
                    df[!, "termination_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{String,Missing}[] 
                    df[!, "objective_value_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[] 
                    df[!, "time_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[] 
                    df[!, "feasibility_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{Bool,Missing}[]
                    df[!, "cuts_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{Int64,Missing}[]
                    df[!, "iter_cb_mis_" * string(mis) * "_density_" * string(density)] = Union{Int64,Missing}[]
                    df[!, "times_master_problem_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    df[!, "times_sub_problem_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    df[!, "times_mis_problem_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    # df[!, "times_filtering_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                end
            end
        end
        if mis_big_m
            df[!, "termination_cb_big_m_mis_" * string(mis)] = Union{String,Missing}[] 
            df[!, "objective_value_cb_big_m_mis_" * string(mis)] = Union{Float64,Missing}[] 
            df[!, "time_cb_big_m_mis_" * string(mis)] = Union{Float64,Missing}[] 
            df[!, "feasibility_cb_big_m_mis_" * string(mis)] = Union{Bool,Missing}[]
            df[!, "cuts_cb_big_m_mis_" * string(mis)] = Union{Int64,Missing}[]
            df[!, "iter_cb_big_m_mis_" * string(mis)] = Union{Int64,Missing}[]
            df[!, "times_master_problem_big_m_mis_" * string(mis)] = Union{Float64,Missing}[]
            df[!, "times_sub_problem_big_m_mis_" * string(mis)] = Union{Float64,Missing}[]
            df[!, "times_mis_problem_big_m_mis_" * string(mis)] = Union{Float64,Missing}[]
            # df[!, "times_filtering_big_m_mis_" * string(mis)] = Union{Float64,Missing}[]
            if distinct_cuts
                df[!, "termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{String,Missing}[] 
                df[!, "objective_value_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[] 
                df[!, "time_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[] 
                df[!, "feasibility_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Bool,Missing}[]
                df[!, "cuts_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Int64,Missing}[]
                df[!, "iter_cb_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Int64,Missing}[]
                df[!, "times_master_problem_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                df[!, "times_sub_problem_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                df[!, "times_mis_problem_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
                # df[!, "times_filtering_big_m_mis_" * string(mis) * "_distinct_cuts"] = Union{Float64,Missing}[]
            end 
            if !isempty(cut_densities)
                for density in cut_densities
                    df[!, "termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{String,Missing}[] 
                    df[!, "objective_value_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[] 
                    df[!, "time_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[] 
                    df[!, "feasibility_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Bool,Missing}[]
                    df[!, "cuts_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Int64,Missing}[]
                    df[!, "iter_cb_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Int64,Missing}[]
                    df[!, "times_master_problem_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    df[!, "times_sub_problem_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    df[!, "times_mis_problem_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                    # df[!, "times_filtering_big_m_mis_" * string(mis) * "_density_" * string(density)] = Union{Float64,Missing}[]
                end
            end
        end 
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

        if ll_fba_indicator
            # read ll-FBA data
            dict = JSON.parse(open("json/" * organism * "_loopless_indicator_fba_" * file))
            dict_organism[:termination_ll_fba_indicator] = dict["termination"]
            dict_organism[:objective_value_ll_fba_indicator] = dict["objective_value"]
            dict_organism[:time_ll_fba_indicator] = dict["time"]
            dict_organism[:feasibility_ll_fba_indicator] = dict["thermo_feasible"]
        end 

        if no_good_cuts
            # read no good cuts data 
            try 
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_" * file))
            catch e 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => NaN,
                    "time" => NaN, 
                    "time_limit" => time_limit,
                    "thermo_feasible" => false,
                    "iter" => missing,
                    "cuts" => missing,
                    "times_master_problem" => NaN,
                    "times_mis_problem" => NaN,
                    "times_sub_problem" => NaN
                )
            end 

            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >= dict["time_limit"]
                    dict_organism[:termination_no_good_cuts] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_no_good_cuts] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_no_good_cuts] = "OPTIMAL"
            else 
                dict_organism[:termination_no_good_cuts] = dict["termination"]
            end
            
            dict_organism[:objective_value_no_good_cuts] = dict["objective_value"]
            dict_organism[:time_no_good_cuts] = dict["time"]
            dict_organism[:feasibility_no_good_cuts] = dict["thermo_feasible"]
            dict_organism[:cuts_no_good_cuts] = dict["cuts"]
            dict_organism[:iter_no_good_iter] = dict["iter"]
        end 

        if no_good_cuts_big_m
            # read no good cuts data 
            try 
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_big_m_" * file))
            catch e 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => NaN,
                    "time" => NaN, 
                    "time_limit" => time_limit,
                    "thermo_feasible" => false,
                    "iter" => missing,
                    "cuts" => missing,
                    "times_master_problem" => NaN,
                    "times_mis_problem" => NaN,
                    "times_sub_problem" => NaN
                )
            end 

            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >= dict["time_limit"]
                    dict_organism[:termination_no_good_cuts_big_m] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_no_good_cuts_big_m] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_no_good_cuts_big_m] = "OPTIMAL"
            else 
                dict_organism[:termination_no_good_cuts_big_m] = dict["termination"]
            end
            
            dict_organism[:objective_value_no_good_cuts_big_m] = dict["objective_value"]
            dict_organism[:time_no_good_cuts_big_m] = dict["time"]
            dict_organism[:feasibility_no_good_cuts_big_m] = dict["thermo_feasible"]
            dict_organism[:cuts_no_good_cuts_big_m] = dict["cuts"]
            dict_organism[:iter_no_good_iter_big_m] = dict["iter"]
        end 

        if cb 
            # read fast CB data
            try 
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * file))
            catch e 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => NaN,
                    "time" => NaN, 
                    "time_limit" => time_limit,
                    "thermo_feasible" => false,
                    "iter" => missing,
                    "cuts" => missing,
                    "times_master_problem" => NaN,
                    "times_mis_problem" => NaN,
                    "times_sub_problem" => NaN
                )
            end 

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
            dict_organism[:cuts_cb] = dict["cuts"]
            dict_organism[:iter_cb] = dict["iter"]
            dict_organism[Symbol("times_master_problem_cb")] = geomean(dict["times_master_problem"])
            dict_organism[Symbol("times_sub_problem_cb")] = geomean(dict["times_sub_problem"])
            dict_organism[Symbol("times_mis_problem_cb")] = geomean(dict["times_mis_problem"])
        end 

        if cb_big_m
            try 
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * file))
            catch e 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => NaN,
                    "time" => NaN, 
                    "time_limit" => time_limit,
                    "thermo_feasible" => false,
                    "iter" => missing,
                    "cuts" => missing,
                    "times_master_problem" => NaN,
                    "times_mis_problem" => NaN,
                    "times_sub_problem" => NaN
                )
            end 

            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >= dict["time_limit"]
                    dict_organism[:termination_cb_big_m] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_cb_big_m] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_cb_big_m] = "OPTIMAL"
            else 
                dict_organism[:termination_cb_big_m] = dict["termination"]
            end

            dict_organism[:objective_value_cb_big_m] = dict["objective_value"]
            dict_organism[:time_cb_big_m] = dict["time"]
            dict_organism[:feasibility_cb_big_m] = dict["thermo_feasible"]
            dict_organism[:cuts_cb_big_m] = dict["cuts"]
            dict_organism[:iter_cb_big_m] = dict["iter"]
            dict_organism[Symbol("times_master_problem_cb_big_m")] = geomean(dict["times_master_problem"])
            dict_organism[Symbol("times_sub_problem_cb_big_m")] = geomean(dict["times_sub_problem"])
            dict_organism[Symbol("times_mis_problem_cb_big_m")] = geomean(dict["times_mis_problem"])
        end 

        if cb_indicator_and_big_m
            try 
                dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_indicator_and_big_m_" * file))
            catch e 
                dict = Dict{String, Any}(
                    "termination" => "ERROR",
                    "objective_value" => NaN,
                    "time" => NaN, 
                    "time_limit" => time_limit,
                    "thermo_feasible" => false,
                    "iter" => missing,
                    "cuts" => missing,
                    "times_master_problem" => NaN,
                    "times_mis_problem" => NaN,
                    "times_sub_problem" => NaN
                )
            end 

            if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                if dict["time"] >= dict["time_limit"]
                    dict_organism[:termination_cb_indicator_and_big_m] = "TIME_LIMIT"
                else 
                    dict_organism[:termination_cb_indicator_and_big_m] = "INFEASIBLE"
                end
            elseif dict["thermo_feasible"] == true
                dict_organism[:termination_cb_indicator_and_big_m] = "OPTIMAL"
            else 
                dict_organism[:termination_cb_indicator_and_big_m] = dict["termination"]
            end

            dict_organism[:objective_value_cb_indicator_and_big_m] = dict["objective_value"]
            dict_organism[:time_cb_indicator_and_big_m] = dict["time"]
            dict_organism[:feasibility_cb_indicator_and_big_m] = dict["thermo_feasible"]
            dict_organism[:cuts_cb_indicator_and_big_m] = dict["cuts"]
            dict_organism[:iter_cb_indicator_and_big_m] = dict["iter"]
            dict_organism[Symbol("times_master_problem_cb_indicator_and_big_m")] = geomean(dict["times_master_problem"])
            dict_organism[Symbol("times_sub_problem_cb_indicator_and_big_m")] = geomean(dict["times_sub_problem"])
            dict_organism[Symbol("times_mis_problem_cb_indicator_and_big_m")] = geomean(dict["times_mis_problem"])
        end 

        if fba 
            dict = JSON.parse(open("json/" * organism * "_fba.json"))
            dict_organism[:termination_fba] = dict["termination"]
            dict_organism[:objective_value_fba] = dict["objective_value"]
            dict_organism[:time_fba] = dict["time"]
            dict_organism[:feasibility_fba] = dict["thermo_feasible"]
        end 

        for mis in mis_numbers
            if mis_indicator
                try 
                    dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * string(mis) * "_mis_" * file))
                catch e 
                    dict = Dict{String, Any}(
                        "termination" => "ERROR",
                        "objective_value" => missing,
                        "time" => NaN, 
                        "time_limit" => time_limit,
                        "thermo_feasible" => missing,
                        "iter" => missing,
                        "cuts" => missing,
                        "times_master_problem" => NaN,
                        "times_sub_problem" => NaN,
                        "times_mis_problem" => NaN
                    )
                end 

                if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                    if dict["time"] >= dict["time_limit"]
                        dict_organism[Symbol("termination_cb_mis_" * string(mis))] = "TIME_LIMIT"
                    else 
                        dict_organism[Symbol("termination_cb_mis_" * string(mis))] = "INFEASIBLE"
                    end
                elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                    dict_organism[Symbol("termination_cb_mis_" * string(mis))] = "OPTIMAL"
                else 
                    dict_organism[Symbol("termination_cb_mis_" * string(mis))] = dict["termination"]
                end
                
                dict_organism[Symbol("objective_value_cb_mis_" * string(mis))] = dict["objective_value"]
                dict_organism[Symbol("time_cb_mis_" * string(mis))] = dict["time"]
                dict_organism[Symbol("feasibility_cb_mis_" * string(mis))] = dict["thermo_feasible"]
                dict_organism[Symbol("cuts_cb_mis_" * string(mis))] = dict["cuts"]
                dict_organism[Symbol("iter_cb_mis_" * string(mis))] = dict["iter"]
                dict_organism[Symbol("times_master_problem_mis_" * string(mis))] = geomean(dict["times_master_problem"])
                dict_organism[Symbol("times_sub_problem_mis_" * string(mis))] = geomean(dict["times_sub_problem"])
                dict_organism[Symbol("times_mis_problem_mis_" * string(mis))] = geomean(dict["times_mis_problem"])
                if distinct_cuts                    
                    try 
                        dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * string(mis) * "_mis_distinct_cuts_" * file))
                    catch e 
                        dict = Dict{String, Any}(
                            "termination" => "ERROR",
                            "objective_value" => missing,
                            "time" => NaN, 
                            "time_limit" => time_limit,
                            "thermo_feasible" => missing,
                            "iter" => missing,
                            "cuts" => missing,
                            "times_master_problem" => NaN,
                            "times_sub_problem" => NaN,
                            "times_mis_problem" => NaN,
                            # "times_filtering" => NaN
                        )
                    end 
    
                    if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                        if dict["time"] >= dict["time_limit"]
                            dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_distinct_cuts")] = "TIME_LIMIT"
                        else 
                            dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_distinct_cuts")] = "INFEASIBLE"
                        end
                    elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                        dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_distinct_cuts")] = "OPTIMAL"
                    else 
                        dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["termination"]
                    end
                    
                    dict_organism[Symbol("objective_value_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["objective_value"]
                    dict_organism[Symbol("time_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["time"]
                    dict_organism[Symbol("feasibility_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["thermo_feasible"]
                    dict_organism[Symbol("cuts_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["cuts"]
                    dict_organism[Symbol("iter_cb_mis_" * string(mis) * "_distinct_cuts")] = dict["iter"]
                    dict_organism[Symbol("times_master_problem_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_master_problem"])
                    dict_organism[Symbol("times_sub_problem_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_sub_problem"])
                    dict_organism[Symbol("times_mis_problem_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_mis_problem"])
                    # dict_organism[Symbol("times_filtering_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_filtering"])
                end 
                if !isempty(cut_densities)
                    for density in cut_densities
                        try 
                            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_" * string(mis) * "_mis_" * string(density) * "_max_density_" * file))
                        catch e 
                            dict = Dict{String, Any}(
                                "termination" => "ERROR",
                                "objective_value" => missing,
                                "time" => NaN, 
                                "time_limit" => time_limit,
                                "thermo_feasible" => missing,
                                "iter" => missing,
                                "cuts" => missing,
                                "times_master_problem" => NaN,
                                "times_sub_problem" => NaN,
                                "times_mis_problem" => NaN,
                                # "times_filtering" => NaN
                            )
                        end 
        
                        if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                            if dict["time"] >= dict["time_limit"]
                                dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_density_" * string(density))] = "TIME_LIMIT"
                            else 
                                dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_density_" * string(density))] = "INFEASIBLE"
                            end
                        elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                            dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_density_" * string(density))] = "OPTIMAL"
                        else 
                            dict_organism[Symbol("termination_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["termination"]
                        end
                        
                        dict_organism[Symbol("objective_value_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["objective_value"]
                        dict_organism[Symbol("time_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["time"]
                        dict_organism[Symbol("feasibility_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["thermo_feasible"]
                        dict_organism[Symbol("cuts_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["cuts"]
                        dict_organism[Symbol("iter_cb_mis_" * string(mis) * "_density_" * string(density))] = dict["iter"]
                        dict_organism[Symbol("times_master_problem_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_master_problem"])
                        dict_organism[Symbol("times_sub_problem_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_sub_problem"])
                        dict_organism[Symbol("times_mis_problem_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_mis_problem"])
                        # dict_organism[Symbol("times_filtering_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_filtering"])
                    end
                end
            end 
            if mis_big_m
                try
                    dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * string(mis) * "_mis_" * file))
                catch e
                    println(e)
                    dict = Dict{String, Any}(
                        "termination" => "ERROR",
                        "objective_value" => missing,
                        "time" => NaN, 
                        "time_limit" => time_limit,
                        "thermo_feasible" => missing,
                        "iter" => missing,
                        "cuts" => missing,
                        "times_master_problem" => NaN,
                        "times_sub_problem" => NaN,
                        "times_mis_problem" => NaN
                    )
                end 

                if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                    if dict["time"] >= dict["time_limit"]
                        dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis))] = "TIME_LIMIT"
                    else 
                        dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis))] = "INFEASIBLE"
                    end
                elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                    dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis))] = "OPTIMAL"
                else 
                    dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis))] = dict["termination"]
                end

                dict_organism[Symbol("objective_value_cb_big_m_mis_" * string(mis))] = dict["objective_value"]
                dict_organism[Symbol("time_cb_big_m_mis_" * string(mis))] = dict["time"]
                dict_organism[Symbol("feasibility_cb_big_m_mis_" * string(mis))] = dict["thermo_feasible"]
                dict_organism[Symbol("cuts_cb_big_m_mis_" * string(mis))] = dict["cuts"]
                dict_organism[Symbol("iter_cb_big_m_mis_" * string(mis))] = dict["iter"]
                dict_organism[Symbol("times_master_problem_big_m_mis_" * string(mis))] = geomean(dict["times_master_problem"])
                dict_organism[Symbol("times_sub_problem_big_m_mis_" * string(mis))] = geomean(dict["times_sub_problem"])
                dict_organism[Symbol("times_mis_problem_big_m_mis_" * string(mis))] = geomean(dict["times_mis_problem"])
                if distinct_cuts                    
                    try 
                        dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * string(mis) * "_mis_distinct_cuts_" * file))
                    catch e 
                        dict = Dict{String, Any}(
                            "termination" => "ERROR",
                            "objective_value" => missing,
                            "time" => NaN, 
                            "time_limit" => time_limit,
                            "thermo_feasible" => missing,
                            "iter" => missing,
                            "cuts" => missing,
                            "times_master_problem" => NaN,
                            "times_sub_problem" => NaN,
                            "times_mis_problem" => NaN,
                            # "times_filtering" => NaN
                        )
                    end 
    
                    if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                        if dict["time"] >= dict["time_limit"]
                            dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = "TIME_LIMIT"
                        else 
                            dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = "INFEASIBLE"
                        end
                    elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                        dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = "OPTIMAL"
                    else 
                        dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["termination"]
                    end
                    
                    dict_organism[Symbol("objective_value_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["objective_value"]
                    dict_organism[Symbol("time_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["time"]
                    dict_organism[Symbol("feasibility_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["thermo_feasible"]
                    dict_organism[Symbol("cuts_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["cuts"]
                    dict_organism[Symbol("iter_cb_big_m_mis_" * string(mis) * "_distinct_cuts")] = dict["iter"]
                    dict_organism[Symbol("times_master_problem_big_m_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_master_problem"])
                    dict_organism[Symbol("times_sub_problem_big_m_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_sub_problem"])
                    dict_organism[Symbol("times_mis_problem_big_m_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_mis_problem"])
                    # dict_organism[Symbol("times_filtering_big_m_mis_" * string(mis) * "_distinct_cuts")] = geomean(dict["times_filtering"])
                end
                if !isempty(cut_densities)
                    for density in cut_densities
                        try 
                            dict = JSON.parse(open("json/" * organism * "_combinatorial_benders_fast_big_m_" * string(mis) * "_mis_" * string(density) * "_max_density_" * file))
                        catch e 
                            dict = Dict{String, Any}(
                                "termination" => "ERROR",
                                "objective_value" => missing,
                                "time" => NaN, 
                                "time_limit" => time_limit,
                                "thermo_feasible" => missing,
                                "iter" => missing,
                                "cuts" => missing,
                                "times_master_problem" => NaN,
                                "times_sub_problem" => NaN,
                                "times_mis_problem" => NaN,
                                # "times_filtering" => NaN
                            )
                        end 
        
                        if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
                            if dict["time"] >= dict["time_limit"]
                                dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = "TIME_LIMIT"
                            else 
                                dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = "INFEASIBLE"
                            end
                        elseif !ismissing(dict["thermo_feasible"]) && dict["thermo_feasible"] == true
                            dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = "OPTIMAL"
                        else 
                            dict_organism[Symbol("termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["termination"]
                        end
                        
                        dict_organism[Symbol("objective_value_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["objective_value"]
                        dict_organism[Symbol("time_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["time"]
                        dict_organism[Symbol("feasibility_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["thermo_feasible"]
                        dict_organism[Symbol("cuts_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["cuts"]
                        dict_organism[Symbol("iter_cb_big_m_mis_" * string(mis) * "_density_" * string(density))] = dict["iter"]
                        dict_organism[Symbol("times_master_problem_big_m_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_master_problem"])
                        dict_organism[Symbol("times_sub_problem_big_m_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_sub_problem"])
                        dict_organism[Symbol("times_mis_problem_big_m_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_mis_problem"])
                        # dict_organism[Symbol("times_filtering_big_m_mis_" * string(mis) * "_density_" * string(density))] = geomean(dict["times_filtering"])
                    end
                end
            end
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
            df[!, name] = [((item >= time_limit) && !isnan(item)) ? time_limit : item for item in df[!, name]]
            df[!, name] = [isnan(item) ? item : Int(round(item, digits=0)) for item in df[!, name]]
            df[!, name] = [isnan(item) ? missing : item for item in df[!, name]]
        end 
    end

    # @show df[!, [:objective_value_ll_fba, :objective_value_no_good_cuts, :objective_value_cb]]
    # @show df[!, [:termination_ll_fba, :termination_no_good_cuts, :termination_cb]]
    # @show df[!, [:time_ll_fba, :time_no_good_cuts, :time_cb]]
    
    if yeast
        if time_limit == 1800
            file_name = "results_yeast_" * solver * ".csv"
        else 
            file_name = "results_yeast_" * solver * "_" * string(time_limit) * ".csv"
        end
    else
        file_name = "results_bigg_" * solver * ".csv"
    end 
    CSV.write("csv/" * file_name, df, append=false, writeheader=true)
end

function sub_csv(file_out, file_in="results_bigg_SCIP.csv"; nullspace=false, ll_fba=false, cb=false, cb_big_m=false, cb_indicator_and_big_m=false, objective_values=true, mis_indicator=false, mis_big_m=false, mis_numbers=[], no_good_cuts=false, no_good_cuts_big_m=false, ch=false, ch_mis=false, ll_fba_indicator=false, cut_densities=[], distinct_cuts=false)
    df = CSV.read("csv/" * file_in, DataFrame)

    cols = ["organism"]

    if ll_fba 
        append!(cols, ["termination_ll_fba", "time_ll_fba"])
        if objective_values
            append!(cols, ["objective_value_ll_fba"])
        end
    end 

    if nullspace 
        append!(cols, ["termination_ll_fba_nullspace", "time_ll_fba_nullspace"])
        if objective_values
            append!(cols, ["objective_value_ll_fba_nullspace"])
        end
    end 
    
    if ll_fba_indicator 
        append!(cols, ["termination_ll_fba_indicator", "time_ll_fba_indicator"])
        if objective_values
            append!(cols, ["objective_value_ll_fba_indicator"])
        end
    end 

    if cb
        append!(cols, ["termination_cb", "time_cb"])
        if objective_values
            append!(cols, ["objective_value_cb"])
        end
    end 
    
    if cb_big_m 
        append!(cols, ["termination_cb_big_m", "time_cb_big_m"])
        if objective_values
            append!(cols, ["objective_value_cb_big_m"])
        end
    end 

    if cb_indicator_and_big_m 
        append!(cols, ["termination_cb_indicator_and_big_m", "time_cb_indicator_and_big_m"])
        if objective_values
            append!(cols, ["objective_value_cb_indicator_and_big_m"])
        end
    end 
    
    if !isempty(mis_numbers)
        for mis in mis_numbers
            if mis_indicator
                append!(cols, ["termination_cb_mis_" * string(mis), "time_cb_mis_" * string(mis)])
                if objective_values
                    append!(cols, ["objective_value_cb_mis_" * string(mis)])
                end
                if mis == 2.0
                    for density in cut_densities
                        append!(cols, ["termination_cb_mis_" * string(mis) * "_density_" * string(density), "time_cb_mis_" * string(mis) * "_density_" * string(density)])
                    end 
                    if distinct_cuts
                        append!(cols, ["termination_cb_mis_" * string(mis) * "_distinct_cuts", "time_cb_mis_" * string(mis) * "_distinct_cuts"])
                    end
                end
            end 
            if mis_big_m
                append!(cols, ["termination_cb_big_m_mis_" * string(mis), "time_cb_big_m_mis_" * string(mis)])
                if objective_values
                    append!(cols, ["objective_value_cb_big_m_mis_" * string(mis)])
                end
                if mis == 0.5
                    for density in cut_densities
                        append!(cols, ["termination_cb_big_m_mis_" * string(mis) * "_density_" * string(density), "time_cb_big_m_mis_" * string(mis) * "_density_" * string(density)])
                    end 
                    if distinct_cuts
                        append!(cols, ["termination_cb_big_m_mis_" * string(mis) * "_distinct_cuts", "time_cb_big_m_mis_" * string(mis) * "_distinct_cuts"])
                    end
                end
            end 
        end 
    end

    if no_good_cuts
        append!(cols, ["termination_no_good_cuts", "time_no_good_cuts"])
        if objective_values
            append!(cols, ["objective_value_no_good_cuts"])
        end
    end 

    if no_good_cuts_big_m
        append!(cols, ["termination_no_good_cuts_big_m", "time_no_good_cuts_big_m"])
        if objective_values
            append!(cols, ["objective_value_no_good_cuts_big_m"])
        end
    end 

    if ch || ch_mis 
        df_ch = CSV.read("csv/results_bigg_ch.csv", DataFrame)
        select!(df_ch, Not(["termination_ll_fba", "feasibility_ll_fba", "objective_value_ll_fba", "time_ll_fba", "time_limit"]))
        df = innerjoin(df, df_ch, on="organism")
        if ch 
            append!(cols, ["termination_ch", "time_ch"])
            if objective_values
                append!(cols, ["objective_value_ch"])
            end
        end 
        if ch_mis 
            append!(cols, ["termination_ch_mis_5", "time_ch_mis_5"])
            if objective_values
                append!(cols, ["objective_value_ch_mis_5"])
            end
        end
    end

    # get model size 
    df_size = CSV.read("../molecular_models/bigg_model_data.csv", DataFrame)
    df = innerjoin(df, df_size, on="organism")
    sort!(df, [:reactions])
    df = df[!, cols]

    CSV.write("csv/" * file_out, df, append=false, writeheader=true)
end 

organisms = [
    "iAF692", # recompute for 1e-5
    "iJR904", 
    "iML1515", 
    "e_coli_core",
    "iNF517",
    "iSB619", # AssertionError("feasible")
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
    "iSbBS512_1146", # recompute on cluster
    "RECON1",
    "Recon3D",
    "STM_v1_0",
    "iAB_RBC_283",
    "iAPECO1_1312",
    "iECB_1328",
    "iETEC_1333",
    "iHN637",
    "iIS312_Amastigote",
    "iJB785",
    "iJN746",
    "iLB1027_lipid",
    "iMM1415",
    "iND750",
    "iRC1080",
    "iSFxv_1172",
    # "iSynCJ816",
    "iYO844",
    "iYS1720",
    "iZ_1308"
]
mis_numbers = [0.1, 0.5, 1.0, 2.0, 5, 10, 20, 30]
cut_densities = [5, 10, 15, 20]
solver_data(organisms, time_limit=1800, yeast=false, cb=true, fba=true, cb_big_m=true, mis_indicator=true, mis_big_m=true, nullspace=true, mis_numbers=mis_numbers, no_good_cuts=true, no_good_cuts_big_m=true, cb_indicator_and_big_m=true, ll_fba_indicator=true, cut_densities=cut_densities, distinct_cuts=true)

# sub_csv("results_ll_fba_variants.csv", nullspace=true, ll_fba=true)
# sub_csv("results_ll_fba_indicator.csv", ll_fba_indicator=true, ll_fba=true)

# solver_data(organisms, time_limit=1800, yeast=false, cb=true, fba=true, cb_big_m=true, mis_indicator=true, mis_big_m=true, nullspace=true, mis_numbers=mis_numbers, no_good_cuts=true, no_good_cuts_big_m=true)
# sub_csv("results_bigg_SCIP_llfba_vs_nullspace.csv", nullspace=true, ll_fba=true)
# sub_csv("results_bigg_SCIP_cb_vs_nullspace.csv", ll_fba=true, cb=true, cb_big_m=true, objective_values=false)
# sub_csv("results_bigg_SCIP_multiple_mis_indicator.csv", mis_indicator=true, mis_big_m=false, objective_values=false, mis_numbers=mis_numbers)
# sub_csv("results_bigg_SCIP_multiple_mis_big_m.csv", cb=true, mis_indicator=false, mis_big_m=true, objective_values=false, mis_numbers=mis_numbers)
# sub_csv("results_bigg_SCIP_no_good_cuts_big_m.csv", ll_fba=true,no_good_cuts_big_m=true, cb_big_m=false, objective_values=false)
# sub_csv("results_bigg_SCIP_no_good_cuts_indicator.csv", no_good_cuts=true, cb=true, objective_values=true)

# sub_csv("results_bigg_SCIP_indicator_and_big_m.csv", cb=true, cb_big_m=true, cb_indicator_and_big_m=true, objective_values=false)
sub_csv("results_bigg_SCIP_cut_selection.csv", ll_fba=true, cb=true, cb_big_m=true, mis_indicator=true, mis_big_m=true, objective_values=false, mis_numbers=mis_numbers, cut_densities=cut_densities, distinct_cuts=true)

# CONSTRAINT HANDLER
# organisms = [
#     "iAF692", # recompute for 1e-5
#     "iJR904", 
#     "iML1515", 
#     "e_coli_core",
#     "iNF517",
#     "iSB619", # AssertionError("feasible")
#     "iNJ661",
#     "iCN900",
#     "iAF1260",
#     "iEK1008",
#     "iJO1366",
#     "iMM904",
#     "iSDY_1059",
#     "iSFV_1184", # recompute on cluster
#     "iSF_1195",
#     "iS_1188",
#     "iSbBS512_1146"
# ]
# sub_csv("results_bigg_ch_vs_cb_big_m.csv", cb_big_m=true, ch=true, ch_mis=true, objective_values=true)

# println(" ")
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
# # mis_numbers = [0.1, 0.2, 0.5, 2.0]
# # solver_data(organisms, time_limit=1800, fba=false, yeast=true, cb=true, cb_big_m=true, mis_numbers=mis_numbers, mis_big_m=true, mis_indicator=true)

# mis_numbers = [0.1, 0.2]
# solver_data(organisms, time_limit=14440, fba=false, yeast=true, cb=true, cb_big_m=true, mis_numbers=mis_numbers, mis_big_m=true, mis_indicator=true)
