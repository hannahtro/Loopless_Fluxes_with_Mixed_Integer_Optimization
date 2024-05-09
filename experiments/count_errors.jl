using DataFrames, JSON, CSV
using Statistics
using Query

include("store_different_solver_results.jl")

# sub_csv("results_ll_fba_variants.csv", nullspace=true, ll_fba=true)
# sub_csv("results_cb_indicator_vs_big_m.csv", cb=true, cb_big_m=true, cb_indicator_and_big_m=true)
# sub_csv_gecko("results_ll_fba_vs_cb_gecko.csv")

function count_errors(;
    intervals=[0, 10, 60, 300, 600, 1200, 1800], 
    solving_strategies=["loopless_fba", "combinatorial_benders_big_m", "combinatorial_benders_big_m_0_5"], save_as="sum_termination_status_gecko.csv",
    read_from="results_ll_fba_vs_cb_gecko.csv",
    bins_column="time_loopless_fba",
    mis_list=[],
    filter_instances=false
)
    df = CSV.read("csv/" * read_from, DataFrame)

    updated_df_names = [replace(name, "_gecko_tol_1.0e-8" => "") for name in names(df)]
    updated_df_names = [replace(name, "_gecko_fast" => "") for name in updated_df_names]
    updated_df_names = [replace(name, "_gecko" => "") for name in updated_df_names]
    updated_df_names = [replace(name, "_tol_1.0e-8" => "") for name in updated_df_names]
    updated_df_names = [replace(name, "." => "_") for name in updated_df_names]
    rename!(df, updated_df_names)
    # print(names(df))

    if filter_instances
        df = df[(df.termination_loopless_fba .!= "INFEASIBLE"), :]
        df = df[(df.termination_loopless_fba .!= "ERROR"), :]
        df = df[(df.termination_loopless_fba .== "OPTIMAL") .& (df.objective_value_loopless_fba .<=  30) .| (df.termination_loopless_fba .== "TIME_LIMIT") , :]
        df = df[(df.termination_loopless_fba .== "OPTIMAL") .& (df.objective_value_loopless_fba .>  0) .| (df.termination_loopless_fba .== "TIME_LIMIT") , :]
        df = df[(df.termination_combinatorial_benders_big_m .== "OPTIMAL") .& (df.objective_value_combinatorial_benders_big_m .<=  30) .| (df.termination_combinatorial_benders_big_m .== "TIME_LIMIT") .| (df.termination_combinatorial_benders_big_m .== "ERROR") , :]
        df = df[(df.termination_combinatorial_benders_big_m .== "OPTIMAL") .& (df.iter_combinatorial_benders_big_m .<=  70) .| (df.termination_combinatorial_benders_big_m .== "TIME_LIMIT") .| (df.termination_combinatorial_benders_big_m .== "ERROR") , :]
        df = df[(df.termination_combinatorial_benders_big_m_0_5 .== "OPTIMAL") .& (df.objective_value_combinatorial_benders_big_m_0_5 .<=  30) .| (df.termination_combinatorial_benders_big_m_0_5 .== "TIME_LIMIT") .| (df.termination_combinatorial_benders_big_m_0_5 .== "ERROR") , :]
        CSV.write("csv/results_ll_fba_vs_cb_gecko_filtered.csv", df, append=false, writeheader=true)
    end

    if !isempty(intervals)
        intervals_dict = Dict()
        for i in 1:length(intervals)-1
            intervals_dict[i] = string(intervals[i]) * "-" * string(intervals[i+1])
        end
        intervals_dict[length(intervals)] = string(intervals[end]) * "-Inf"

        @show intervals_dict

        # 7 bins
        if bins_column == "time_loopless_fba"
            df[!, "bins"] = zeros(nrow(df))
            for idx in 1:length(intervals)-1
                df.bins = trunc.(Int64, max.(idx * interval.(df.time_loopless_fba, intervals[idx], intervals[idx+1]), df.bins))
            end
            df.bins = trunc.(Int64, max.(length(intervals) * interval.(df.time_loopless_fba, intervals[end], Inf), df.bins))
        elseif bins_column == "time_ll_fba"
            df[!, "bins"] = zeros(nrow(df))
            for idx in 1:length(intervals)-1
                df.bins = trunc.(Int64, max.(idx * interval.(df.time_ll_fba, intervals[idx], intervals[idx+1]), df.bins))
            end
            df.bins = trunc.(Int64, max.(length(intervals) * interval.(df.time_ll_fba, intervals[end], Inf), df.bins))
        elseif bins_column == "time_ll_fba_indicator"
            df[!, "bins"] = zeros(nrow(df))
            for idx in 1:length(intervals)-1
                df.bins = trunc.(Int64, max.(idx * interval.(df.time_ll_fba_indicator, intervals[idx], intervals[idx+1]), df.bins))
            end
            df.bins = trunc.(Int64, max.(length(intervals) * interval.(df.time_ll_fba_indicator, intervals[end], Inf), df.bins))
        end
    else 
        intervals_dict = Dict()
        intervals_dict[1] = "all"
        df[!,"bins"] = ones(nrow(df))
    end 

    for strategy in solving_strategies 
        if occursin("cb_big_m", strategy)
            df[(df.organism .== "iRC1080"), Symbol("termination_" * strategy)] = ["INFEASIBLE"]
            df[(df.organism .== "iSbBS512_1146"), Symbol("termination_" * strategy)] = ["TIME_LIMIT"]
        end 
        if occursin("cb_mis", strategy) || strategy == "cb"     
            df[(df.organism .== "iSB619"), Symbol("termination_" * strategy)] = ["INFEASIBLE"]
        end

        df[!, "optimal_" * strategy] = [termination=="OPTIMAL" ? 1 : 0 for termination in df[!, "termination_" * strategy]]
        df[!, "timelimit_" * strategy] = [termination=="TIME_LIMIT" ? 1 : 0 for termination in df[!, "termination_" * strategy]]
        df[!, "error_" * strategy] = [(termination=="ERROR" || termination=="INFEASIBLE") ? 1 : 0 for termination in df[!, "termination_" * strategy]]
    end

    # print(df[!,:bins])

    # groupby bins and count termination status for different solving strategies
    gdf = groupby(df, [:bins])
    gdf = combine(gdf, nrow => :count, append!(append!([
        Symbol("optimal_" * strategy) for strategy in solving_strategies
        ], [Symbol("timelimit_" * strategy) for strategy in solving_strategies]), Symbol("error_" * strategy) for strategy in solving_strategies) .=> sum)

    gdf[!,"bins_string"] = [intervals_dict[i] for i in gdf[!,"bins"]]

    # print(gdf)
    if isempty(intervals)
        if !contains(read_from, "yeast") && !(contains(read_from, "gecko"))
            new_df = DataFrame(mis=String[], optimal_ll_fba=Int64[], timelimit_ll_fba=Int64[], optimal_cb=Int64[], timelimit_cb=Int64[], error_cb=Int64[])
            for mis in mis_list
                if bins_column=="time_ll_fba_indicator"
                    push!(new_df, (mis, gdf[!, "optimal_ll_fba_indicator_sum"][1], gdf[!, "timelimit_ll_fba_indicator_sum"][1], gdf[!,"optimal_cb_mis_" * string(mis) * "_sum"][1], gdf[!,"timelimit_cb_mis_" * string(mis) * "_sum"][1], gdf[!,"error_cb_mis_" * string(mis) * "_sum"][1]))
                else 
                    push!(new_df, (mis, gdf[!, "optimal_ll_fba_sum"][1], gdf[!, "timelimit_ll_fba_sum"][1], gdf[!,"optimal_cb_big_m_mis_" * string(mis) * "_sum"][1], gdf[!,"timelimit_cb_big_m_mis_" * string(mis) * "_sum"][1], gdf[!,"error_cb_big_m_mis_" * string(mis) * "_sum"][1]))
                end
            end
        elseif contains(read_from, "yeast")
            new_df = DataFrame(solving_strategy=String[], optimal=Int64[], timelimit=Int64[], error=Int64[])
            print(names(gdf))
            push!(new_df, ("ll_fba", gdf[!, "optimal_ll_fba_sum"][1], gdf[!, "timelimit_ll_fba_sum"][1], gdf[!, "error_ll_fba_sum"][1]))
            push!(new_df, ("cb_big_m", gdf[!,"optimal_cb_big_m_sum"][1], gdf[!,"timelimit_cb_big_m_sum"][1], gdf[!,"error_cb_big_m_sum"][1]))
            push!(new_df, ("cb_indicator", gdf[!,"optimal_cb_sum"][1], gdf[!,"timelimit_cb_sum"][1], gdf[!,"error_cb_sum"][1]))
            for mis in mis_list
                push!(new_df, (mis, gdf[!,"optimal_cb_big_m_mis_" * string(mis) * "_sum"][1], gdf[!,"timelimit_cb_big_m_mis_" * string(mis) * "_sum"][1], gdf[!,"error_cb_big_m_mis_" * string(mis) * "_sum"][1]))
            end
            for mis in mis_list
                push!(new_df, (mis, gdf[!,"optimal_cb_mis_" * string(mis) * "_sum"][1], gdf[!,"timelimit_cb_mis_" * string(mis) * "_sum"][1], gdf[!,"error_cb_mis_" * string(mis) * "_sum"][1]))
            end
        elseif contains(read_from, "gecko")
            print(names(gdf))
            new_df = DataFrame(solving_strategy=String[], 
            optimal=Int64[], timelimit=Int64[], error=Int64[])
            print(names(gdf))
            push!(new_df, ("ll_fba", gdf[!, "optimal_loopless_fba_sum"][1], gdf[!, "timelimit_loopless_fba_sum"][1], gdf[!, "error_loopless_fba_sum"][1]))
            push!(new_df, ("cb_big_m", gdf[!,"optimal_combinatorial_benders_big_m_sum"][1], gdf[!,"timelimit_combinatorial_benders_big_m_sum"][1], gdf[!,"error_combinatorial_benders_big_m_sum"][1]))
            push!(new_df, ("cb_big_m_0_5", gdf[!,"optimal_combinatorial_benders_big_m_0_5_sum"][1], gdf[!,"timelimit_combinatorial_benders_big_m_0_5_sum"][1], gdf[!,"error_combinatorial_benders_big_m_0_5_sum"][1]))
        end
        gdf=new_df
    end 

    CSV.write("csv/" * save_as , gdf, append=false, writeheader=true)
end 

function interval(val, lb, ub)
    if (val >= lb) && (val <= ub)
        return true
    else 
        return false 
    end
end

count_errors(
    intervals=[], 
    solving_strategies=[
        "loopless_fba", 
        "combinatorial_benders_big_m", "combinatorial_benders_big_m_0_5"
    ], 
    save_as="sum_termination_status_gecko.csv",
    read_from="results_ll_fba_vs_cb_gecko.csv",
    bins_column="time_loopless_fba",
    mis_list=[],
    filter_instances=true
)

# count_errors(
#     intervals=[0, 10, 600, 1800], 
#     solving_strategies=["ll_fba", "ll_fba_indicator", "cb_big_m", "cb"], 
#     save_as="sum_termination_status_bigg.csv",
#     read_from="results_bigg_SCIP.csv",
#     bins_column="time_ll_fba"
# )

# count_errors(
#     intervals=[0, 10, 600, 1800], 
#     solving_strategies=[
#         "ll_fba", 
#         "cb_big_m", 
#         "cb_big_m_mis_0_1", 
#         "cb_big_m_mis_0_5", 
#         "cb_big_m_mis_1_0", 
#         "cb_big_m_mis_2_0", 
#         "cb_big_m_mis_5_0", 
#         "cb_big_m_mis_10_0", 
#         "cb_big_m_mis_20_0", 
#         "cb_big_m_mis_30_0"], 
#     save_as="sum_termination_status_mis_big_m.csv",
#     read_from="results_bigg_SCIP.csv",
#     bins_column="time_ll_fba"
# )

# count_errors(
#     intervals=[], 
#     solving_strategies=[
#         "ll_fba_indicator", 
#         "cb", 
#         "cb_mis_0_1", 
#         "cb_mis_0_5", 
#         "cb_mis_1_0", 
#         "cb_mis_2_0", 
#         "cb_mis_5_0", 
#         "cb_mis_10_0", 
#         "cb_mis_20_0", 
#         "cb_mis_30_0"], 
#     save_as="sum_termination_status_mis_indicator.csv",
#     read_from="results_bigg_SCIP.csv",
#     bins_column="time_ll_fba_indicator",
#     mis_list=["0_1", "0_5", "1_0", "2_0", "5_0", "10_0", "20_0", "30_0"]
# )

# count_errors(
#     intervals=[], 
#     solving_strategies=[
#         "ll_fba", 
#         "cb", 
#         "cb_big_m_mis_0_1", 
#         "cb_big_m_mis_0_5", 
#         "cb_big_m_mis_1_0", 
#         "cb_big_m_mis_2_0", 
#         "cb_big_m_mis_5_0", 
#         "cb_big_m_mis_10_0", 
#         "cb_big_m_mis_20_0", 
#         "cb_big_m_mis_30_0"], 
#     save_as="sum_termination_status_mis_big_m.csv",
#     read_from="results_bigg_SCIP.csv",
#     bins_column="time_ll_fba",
#     mis_list=["0_1", "0_5", "1_0", "2_0", "5_0", "10_0", "20_0", "30_0"]
# )

# count_errors(
#     intervals=[], 
#     solving_strategies=[
#         "ll_fba", 
#         # "ll_fba_indicator", 
#         "cb_big_m", 
#         "cb",
#         "cb_big_m_mis_0_1", 
#         "cb_big_m_mis_0_2", 
#         "cb_big_m_mis_0_5", 
#         "cb_big_m_mis_2_0", 
#         "cb_mis_0_1", 
#         "cb_mis_0_2", 
#         "cb_mis_0_5", 
#         "cb_mis_2_0", 
#     ], 
#     save_as="sum_termination_status_yeast.csv",
#     read_from="results_yeast_SCIP.csv",
#     bins_column="time_ll_fba",
#     mis_list=["0_1", "0_2", "0_5", "2_0",]
# )