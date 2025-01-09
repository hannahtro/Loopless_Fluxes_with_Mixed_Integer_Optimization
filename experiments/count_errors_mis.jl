using DataFrames, JSON, CSV
using Statistics
# using Query

include("store_different_solver_results.jl")

# sub_csv("results_ll_fba_variants.csv", nullspace=true, ll_fba=true)
# sub_csv("results_cb_indicator_vs_big_m.csv", cb=true, cb_big_m=true, cb_indicator_and_big_m=true)
# sub_csv_gecko("results_ll_fba_vs_cb_gecko.csv")

function count_errors(;
    intervals=[0, 10, 60, 300, 600, 1200, 1800], 
    solving_strategies=["loopless_fba", "combinatorial_benders_big_m", "combinatorial_benders_big_m_0_5"], save_as="sum_termination_status_gecko.csv",
    read_from="results_ll_fba_vs_cb_gecko.csv",
    mis_list=[],
    filter_instances=false
)
    df = CSV.read("csv/" * read_from, DataFrame)

    updated_df_names = [replace(name, "." => "_") for name in names(df)]
    rename!(df, updated_df_names)
    # print(names(df))

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

    # set runtime to 1800 if status not optimal
    df[df.optimal_ll_fba_indicator .!= 1, :time_ll_fba_indicator] .= 1800
    df[df.optimal_ll_fba .!= 1, :time_ll_fba] .= 1800
    df[df.optimal_cb .!= 1, :time_cb] .= 1800
    df[df.optimal_cb_big_m .!= 1, :time_cb_big_m] .= 1800
    for mis in mis_list
        df[df[!,"optimal_cb_mis_" * string(mis)] .!= 1, "time_cb_mis_" * string(mis)] .= 1800
        df[df[!,"optimal_cb_big_m_mis_" * string(mis)] .!= 1, "time_cb_big_m_mis_" * string(mis)] .= 1800
    end

    @infiltrate

    new_df = DataFrame(
        instances=Int[], 
        mis=String[], 
        optimal_ll_fba=Int64[], 
        timelimit_ll_fba=Int64[], 
        optimal_ll_fba_indicator=Int64[], 
        timelimit_ll_fba_indicator=Int64[], 
        error_ll_fba_indicator=Int64[], 
        optimal_cb=Int64[], 
        time_cb_geom_shifted_mean=Float64[],
        optimal_cb_perc=Float64[],
        timelimit_cb=Int64[], 
        error_cb=Int64[], 
        optimal_cb_big_m=Int64[], 
        time_cb_big_m_geom_shifted_mean=Float64[],
        optimal_cb_big_m_perc=Float64[],
        timelimit_cb_big_m=Int64[], 
        error_cb_big_m=Int64[]
    )
    for (idx,mis) in enumerate(mis_list)
        push!(new_df, [NaN for col in names(new_df)], promote=true)

        new_df[idx, :][:instances] = length(df[!, :organism])
        new_df[idx, :][:mis] = mis
        new_df[idx, :]["optimal_ll_fba"] = sum(df[!, "optimal_ll_fba"])
        new_df[idx, :]["timelimit_ll_fba"] = sum(df[!, "timelimit_ll_fba"])
        new_df[idx, :]["optimal_ll_fba_indicator"] = sum(df[!, "optimal_ll_fba_indicator"])
        new_df[idx, :]["timelimit_ll_fba_indicator"] = sum(df[!, "timelimit_ll_fba_indicator"])
        new_df[idx, :]["error_ll_fba_indicator"] = sum(df[!, "error_ll_fba_indicator"])
        new_df[idx, :]["optimal_cb"] = sum(df[!,"optimal_cb_mis_" * string(mis)])
        new_df[idx, :]["time_cb_geom_shifted_mean"] = geom_shifted_mean(df[!, "time_cb_mis_" * string(mis)])
        new_df[idx, :]["optimal_cb_perc"] = round(sum(df[!,"optimal_cb_mis_" * string(mis)]) ./ length(df[!, :organism]) * 100, digits=2)
        new_df[idx, :]["timelimit_cb"] = sum(df[!,"timelimit_cb_mis_" * string(mis)])
        new_df[idx, :]["error_cb"] = sum(df[!,"error_cb_mis_" * string(mis)])
        new_df[idx, :]["optimal_cb_big_m"] = sum(df[!,"optimal_cb_big_m_mis_" * string(mis)])
        new_df[idx, :]["time_cb_big_m_geom_shifted_mean"] = geom_shifted_mean(df[!, "time_cb_big_m_mis_" * string(mis)])
        new_df[idx, :]["optimal_cb_big_m_perc"] = round(sum(df[!,"optimal_cb_big_m_mis_" * string(mis)]) ./ length(df[!, :organism]) * 100, digits=2)
        new_df[idx, :]["timelimit_cb_big_m"] = sum(df[!,"timelimit_cb_big_m_mis_" * string(mis)])
        new_df[idx, :]["error_cb_big_m"] =  sum(df[!,"error_cb_big_m_mis_" * string(mis)])
    end

    # add CB column
    push!(new_df, [NaN for col in names(new_df)], promote=true)
    new_df[length(mis_list)+1, :][:mis] = 0
    # new_df[length(mis_list)+1, :][:optimal_ll_fba_perc] = round(sum(df[!,"optimal_ll_fba"]) ./ length(df[!, :organism]) * 100, digits=2)
    # new_df[length(mis_list)+1, :][:optimal_ll_fba_indicator_perc] = round(sum(df[!,"optimal_ll_fba_indicator"]) ./ length(df[!, :organism]) * 100, digits=2)
    new_df[length(mis_list)+1, :][:optimal_cb_big_m_perc] = round(sum(df[!,"optimal_cb_big_m"]) ./ length(df[!, :organism]) * 100, digits=2)
    new_df[length(mis_list)+1, :][:optimal_cb_perc] = round(sum(df[!,"optimal_cb"]) ./ length(df[!, :organism]) * 100, digits=2)
    # new_df[length(intervals)+1, :][:time_ll_fba_geom_shifted_mean] = geom_shifted_mean(df[!, :time_ll_fba])
    # new_df[length(intervals)+1, :][:time_ll_fba_indicator_geom_shifted_mean] = geom_shifted_mean(df[!, :time_cb_big_m])
    new_df[length(mis_list)+1, :][:time_cb_geom_shifted_mean] = geom_shifted_mean(df[!, :time_cb])
    new_df[length(mis_list)+1, :][:time_cb_big_m_geom_shifted_mean] = geom_shifted_mean(df[!, :time_cb_big_m])

    new_df = new_df[!, [:mis, :optimal_cb_perc, :time_cb_geom_shifted_mean, :optimal_cb_big_m_perc, :time_cb_big_m_geom_shifted_mean]]
    CSV.write("csv/" * save_as , new_df, append=false, writeheader=true)
end 

function interval(val, lb, ub)
    if (val >= lb) && (val <= ub)
        return true
    else 
        return false 
    end
end

function geom_shifted_mean(xs; shift=big"1.0")
    n = length(xs)
    r = prod(xi + shift for xi in xs)
    return round(Float64(r^(1/n) - shift), digits=2)
end

count_errors(
    intervals=[], 
    solving_strategies=[
        "ll_fba",
        "ll_fba_indicator", 
        "cb", 
        "cb_mis_0_1", 
        "cb_mis_0_5", 
        "cb_mis_2_0", 
        "cb_mis_5_0", 
        "cb_mis_10_0", 
        "cb_big_m",
        "cb_big_m_mis_0_1", 
        "cb_big_m_mis_0_5", 
        "cb_big_m_mis_2_0", 
        "cb_big_m_mis_5_0", 
        "cb_big_m_mis_10_0", 
    ], 
    save_as="sum_termination_status_mis.csv",
    read_from="results_bigg_SCIP.csv",
    mis_list=["0_1", "0_5", "2_0", "5_0", "10_0"]
)
