using DataFrames
using CSV

df = DataFrame(type = String[], objective_value = Float64[], termination = String[], num_blocked_cycles = Int64[], time = Float64[], nodes = Int64[])

file_names = ["iAF692_loopless_fba_blocked_50_1800_10000", "iAF692_loopless_fba_blocked_100_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000"]
for name in file_names
    df_temp = first(CSV.read("iAF692/server/" * name * ".csv", DataFrame),1)
    df_temp = df_temp[!,[:objective_value, :time, :termination, :nodes, :num_blocked_cycles]]
    name = replace(name, "iAF692_" => "")
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

file_name = "iAF692/server/comparison_shortest_cycles.csv"
CSV.write(file_name, df, append=false, writeheader=true)


df = DataFrame(type = String[], objective_value = Float64[], termination = String[], num_blocked_cycles = Int64[], time = Float64[], nodes = Int64[])
file_names = ["iJR904_loopless_fba_blocked_50_1800_10000", "iJR904_loopless_fba_blocked_100_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_100_1800_10000"]
for name in file_names
    df_temp = first(CSV.read("iJR904/server/" * name * ".csv", DataFrame),1)
    df_temp = df_temp[!,[:objective_value, :time, :termination, :nodes, :num_blocked_cycles]]
    name = replace(name, "iJR904_" => "")
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

file_name = "iJR904/server/comparison_shortest_cycles.csv"
CSV.write(file_name, df, append=false, writeheader=true)