using DataFrames, JSON, CSV
using Statistics
# using Query

# include("store_different_solver_results.jl")
# sub_csv("results_ll_fba_variants.csv", nullspace=true, ll_fba=true)
# sub_csv("results_cb_indicator_vs_big_m.csv", cb=true, cb_big_m=true, cb_indicator_and_big_m=true)
# sub_csv_gecko("results_ll_fba_vs_cb_gecko.csv")

function count_errors(;
    intervals=[0, 10, 60, 300, 600, 1200, 1800], 
    solving_strategies=["loopless_fba", "combinatorial_benders_big_m", "combinatorial_benders_big_m_0_5"], 
    save_as="sum_termination_status_gecko.csv",
    read_from="results_ll_fba_vs_cb_gecko.csv",
    bins_column="time_loopless_fba",
    mis_list=[],
    filter_instances=false
)
    df = CSV.read("csv/" * read_from, DataFrame)

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
    end

    # TODO: recheck instances
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

    # @infiltrate
    # filter(row -> row.termination_loopless_fba == "ERROR", df)[!, [:organism]]
    # debug_df = filter(row -> row.bins == 3, df)[!, [:organism, :optimal_ll_fba, :time_ll_fba, :optimal_ll_fba_indicator, :time_ll_fba_indicator]]
    # CSV.write("csv/debug_bins_data.csv" , debug_df, append=false, writeheader=true)

    # groupby bins and count termination status for different solving strategies
    gdf = groupby(df, [:bins])
    gdf_sum = combine(gdf, nrow => :count, append!(
        append!([Symbol("optimal_" * strategy) for strategy in solving_strategies], [Symbol("timelimit_" * strategy) for strategy in solving_strategies]), Symbol("error_" * strategy) for strategy in solving_strategies) .=> sum)

    # groupby bins and apply geom_shifted_mean
    # TODO: check missing instances
    replace!(df.time_cb_big_m, missing => 1800)
    gdf = groupby(df, [:bins])
    gdf_geom_mean = combine(gdf, nrow, ["time_ll_fba", "time_ll_fba_indicator", "time_cb", "time_cb_big_m"] .=> geom_shifted_mean)

    gdf = innerjoin(gdf_sum, gdf_geom_mean, on = :bins)
    gdf[!,"bins_string"] = [intervals_dict[i] for i in gdf[!,"bins"]]

    # add row with data for all instances
    push!(gdf, [NaN for col in names(gdf)], promote=true)
    gdf[length(intervals)+1, :][:bins] = length(intervals)+1
    gdf[length(intervals)+1, :][:bins_string] = "all"
    gdf[length(intervals)+1, :][:count] = size(df)[1]
    gdf[length(intervals)+1, :][:optimal_ll_fba_sum] = sum(df[!, :optimal_ll_fba])
    gdf[length(intervals)+1, :][:optimal_ll_fba_indicator_sum] = sum(df[!, :optimal_ll_fba_indicator])
    gdf[length(intervals)+1, :][:optimal_cb_sum] = sum(df[!, :optimal_cb])
    gdf[length(intervals)+1, :][:optimal_cb_big_m_sum] = sum(df[!, :optimal_cb_big_m])
    gdf[length(intervals)+1, :][:time_ll_fba_geom_shifted_mean] = geom_shifted_mean(df[!, :time_ll_fba])
    gdf[length(intervals)+1, :][:time_ll_fba_indicator_geom_shifted_mean] = geom_shifted_mean(df[!, :time_ll_fba_indicator])
    gdf[length(intervals)+1, :][:time_cb_geom_shifted_mean] = geom_shifted_mean(df[!, :time_cb])
    gdf[length(intervals)+1, :][:time_cb_big_m_geom_shifted_mean] = geom_shifted_mean(df[!, :time_cb_big_m])

    # calculate percentage of optimally solved instances
    if "optimal_ll_fba_sum" in names(gdf)
        gdf[!, "optimal_ll_fba_perc"] = round.(gdf[!, "optimal_ll_fba_sum"] ./ gdf[!, :count], digits=2) * 100
        gdf[!, "optimal_cb_big_m_perc"] = round.(gdf[!, "optimal_cb_big_m_sum"] ./ gdf[!, :count], digits=2) * 100
        for mis in mis_list
            gdf[!, "optimal_cb_big_m_mis_" * string(mis) * "_perc"] = round.(gdf[!, "optimal_cb_big_m_mis_" * string(mis) * "_sum"] ./ gdf[!, :count], digits=2) * 100 
        end
    end 
    if "optimal_ll_fba_indicator_sum" in names(gdf)
        gdf[!, "optimal_ll_fba_indicator_perc"] = round.(gdf[!, "optimal_ll_fba_indicator_sum"] ./ gdf[!, :count], digits=2) * 100
        gdf[!, "optimal_cb_perc"] = round.(gdf[!, "optimal_cb_sum"] ./ gdf[!, :count], digits=2) * 100
        for mis in mis_list
            gdf[!, "optimal_cb_mis_" * string(mis) * "_perc"] = round.(gdf[!, "optimal_cb_mis_" * string(mis) * "_sum"] ./ gdf[!, :count], digits=2) * 100
        end
    end

    # @infiltrate
    # CSV.write("csv/debug_bins_data.csv" , gdf, append=false, writeheader=true)

    gdf = gdf[!, [:bins_string, :count, :optimal_ll_fba_perc, :time_ll_fba_geom_shifted_mean, :optimal_ll_fba_indicator_perc, :time_ll_fba_indicator_geom_shifted_mean, :optimal_cb_big_m_perc, :time_cb_big_m_geom_shifted_mean, :optimal_cb_perc, :time_cb_geom_shifted_mean]]

    CSV.write("csv/" * save_as , gdf, append=false, writeheader=true)
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
    intervals=[0, 10, 600, 1800], 
    solving_strategies=["ll_fba", "ll_fba_indicator", "cb_big_m", "cb"], 
    save_as="sum_termination_status_bigg.csv",
    read_from="results_bigg_SCIP.csv",
    bins_column="time_ll_fba"
)
