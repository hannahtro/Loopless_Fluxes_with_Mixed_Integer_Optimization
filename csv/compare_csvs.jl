using DataFrames
using CSV

function compare_shortest_cycles(;file_names, organism)
    df = DataFrame(type = String[], objective_value = Float64[], termination = String[], num_blocked_cycles = Int64[], time = Float64[], nodes = Int64[])

    for name in file_names
        df_temp = first(CSV.read(organism * "/server/" * name * ".csv", DataFrame),1)
        df_temp = df_temp[!,[:objective_value, :time, :termination, :nodes, :num_blocked_cycles]]
        name = replace(name, organism * "_" => "")
        name = replace(name, "_50_1800_10000" => "")
        name = replace(name, "_100_1800_10000" => "")
        df_temp[!, :type] = [name]
        append!(df, df_temp)
    end 

    df[!,:objective_value] = round.(df[!,:objective_value], digits=4)
    df[!,:time] = round.(df[!,:time], digits=2)
    rename!(df,:objective_value => :objective)
    rename!(df,:time => :time_s)
    rename!(df,:num_blocked_cycles => :blocked_cycles)

    file_name = organism * "/server/comparison_shortest_cycles.csv"
    CSV.write(file_name, df, append=false, writeheader=true)
end

# block_limit, time_limit, ceiling
# file_names = ["iAF692_loopless_fba_blocked_50_1800_10000", "iAF692_loopless_fba_blocked_100_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000"]
# file_names = ["iJR904_loopless_fba_blocked_50_1800_10000", "iJR904_loopless_fba_blocked_100_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_100_1800_10000"]

organism = "iJR904"
file_names = ["iJR904_loopless_fba_1800", "iJR904_loopless_indicator_fba_mu_1800", "iJR904_loopless_indicator_fba_blocked_mu_1800_10_same_objective", "iJR904_loopless_fba_blocked_1800_50", "iJR904_loopless_fba_blocked_100_1800_100", "iJR904_loopless_fba_blocked_100_1800_200", "iJR904_loopless_fba_blocked_100_1800_500", "iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50", "iJR904_loopless_fba_blocked_1800_50_same_objective", "iJR904_loopless_fba_blocked_100_1800_100_same_objective", "iJR904_loopless_fba_blocked_100_1800_200_same_objective", "iJR904_loopless_fba_blocked_100_1800_500_same_objective"]

df = DataFrame(type = String[], objective_value = Float64[], termination = String[], num_blocked_cycles = Union{Missing, Int64}[], time = Float64[], nodes = Int64[])

for name in file_names
    df_temp = first(CSV.read(organism * "/server/" * name * ".csv", DataFrame),1)
    try 
        df_temp = df_temp[!,[:objective_value, :time, :termination, :nodes, :num_blocked_cycles]]
    catch 
        println("num_blocked_cycles not in " * name)
        df_temp = df_temp[!,[:objective_value, :time, :termination, :nodes]]
        df_temp[!, :num_blocked_cycles] = missings(Int, 1)
    end
    name = replace(name, organism * "_" => "")
    # name = replace(name, "_50_1800_10000" => "")
    # name = replace(name, "_100_1800_10000" => "")
    name = replace(name, "_1800" => "")
    df_temp[!, :type] = [name]
    append!(df, df_temp)
end 

df[!,:objective_value] = round.(df[!,:objective_value], digits=4)
df[!,:time] = round.(df[!,:time], digits=2)
rename!(df,:objective_value => :objective)
rename!(df,:time => :time_s)
rename!(df,:num_blocked_cycles => :blocked_cycles)

@show df
file_name = organism * "/server/comparison_blocked_cycles.csv"
CSV.write(file_name, df, append=false, writeheader=true)