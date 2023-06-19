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

# file_names = ["iAF692_loopless_fba_blocked_50_1800_10000", "iAF692_loopless_fba_blocked_100_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000"]
# file_names = ["iJR904_loopless_fba_blocked_50_1800_10000", "iJR904_loopless_fba_blocked_100_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000", "iJR904_loopless_fba_blocked_shortest_cycles_100_1800_10000"]

function compare_blocked_cycles(;file_names, organism)
    df = DataFrame(type = String[], objective_value = Float64[], dual_bound = Float64[], termination = String[], num_blocked_cycles = Union{Missing, Int64}[], ceiling = Union{Missing, Int64}[], block_limit = Union{Missing, Int64}[], time_limit= Int64[], time = Float64[], nodes = Int64[])

    for name in file_names
        @show name
        df_temp = first(CSV.read(organism * "/server/" * name * ".csv", DataFrame),1)
        try 
            df_temp = df_temp[!,[:objective_value, :dual_bound, :time, :termination, :nodes, :num_blocked_cycles, :block_limit, :ceiling, :time_limit]]
        catch 
            println("num_blocked_cycles not in " * name)
            df_temp = df_temp[!,[:objective_value, :dual_bound, :time, :termination, :nodes, :time_limit]]
            df_temp[!, :num_blocked_cycles] = missings(Int, 1)
            df_temp[!, :ceiling] = missings(Int, 1)
            df_temp[!, :block_limit] = missings(Int, 1)
        end
        name = replace(name, organism * "_" => "")
        name = replace(name, "_500" => "")
        name = replace(name, "_50" => "")
        name = replace(name, "_100" => "")
        name = replace(name, "_200" => "")
        name = replace(name, "_1800" => "")
        df_temp[!, :type] = [name]
        append!(df, df_temp)
    end 

    df[!,:objective_value] = round.(df[!,:objective_value], digits=4)
    df[!,:dual_bound] = round.(df[!,:dual_bound], digits=4)
    df[!,:time] = round.(df[!,:time], digits=2)
    rename!(df,:objective_value => :objective)
    rename!(df,:time => :time_s)
    rename!(df,:num_blocked_cycles => :blocked_cycles)

    @show df
    file_name = organism * "/server/comparison_blocked_cycles.csv"
    CSV.write(file_name, df, append=false, writeheader=true)
end

function compare_loopless_formulation(; file_names, organisms=["iAF692","iJR904","iML1515"])
    df = DataFrame(organism = String[], type = String[], objective_value = Float64[], dual_bound = Float64[], termination = String[], time_limit= Int64[], time = Float64[], nodes = Int64[])

    for organism in organisms
        for name in file_names
            @show name
            df_temp = first(CSV.read(organism * "/server/" * organism * "_" * name * ".csv", DataFrame),1)
            try 
                df_temp = df_temp[!,[:objective_value, :dual_bound, :time, :termination, :nodes, :num_blocked_cycles, :block_limit, :ceiling, :time_limit]]
            catch 
                println("num_blocked_cycles not in " * name)
                df_temp = df_temp[!,[:objective_value, :dual_bound, :time, :termination, :nodes, :time_limit]]
            end
            name = replace(name, organism * "_" => "")
            name = replace(name, "_1800" => "")
            df_temp[!, :type] = [name]
            df_temp[!, :organism] = [organism]
            append!(df, df_temp)
        end 
    end

    df[!,:objective_value] = round.(df[!,:objective_value], digits=4)
    df[!,:dual_bound] = round.(df[!,:dual_bound], digits=4)
    df[!,:time] = round.(df[!,:time], digits=2)
    rename!(df,:objective_value => :objective)
    rename!(df,:time => :time_s)

    @show df
    file_name = "comparison_loopless_formulation.csv"
    CSV.write(file_name, df, append=false, writeheader=true)
end

compare_loopless_formulation(file_names=["loopless_fba_1800", "loopless_fba_nullspace_1800", "loopless_indicator_fba_1800"])

# organism = "iAF692"
# file_names = [
#     "iAF692_loopless_fba_1800", 
#     "iAF692_loopless_fba_nullspace_1800", 
#     "iAF692_loopless_indicator_fba_1800", 
#     "iAF692_loopless_indicator_fba_blocked_500_1800_50_same_objective", 
#     "iAF692_loopless_indicator_fba_blocked_500_1800_100_same_objective", #TODO: differet objective
#     "iAF692_loopless_fba_blocked_100_1800_50_same_objective", 
#     "iAF692_loopless_fba_blocked_100_1800_100_same_objective", 
#     "iAF692_loopless_fba_blocked_100_1800_200_same_objective", 
#     "iAF692_loopless_fba_blocked_100_1800_500_same_objective",
#     "iAF692_loopless_fba_blocked_100_1800_50", 
#     "iAF692_loopless_fba_blocked_100_1800_100", 
#     "iAF692_loopless_fba_blocked_100_1800_200", 
#     "iAF692_loopless_fba_blocked_100_1800_500", 
#     "iAF692_loopless_fba_blocked_50_1800_10000", 
#     "iAF692_loopless_fba_blocked_100_1800_10000", 
#     "iAF692_loopless_fba_blocked_shortest_cycles_50_1800_10000", 
#     "iAF692_loopless_fba_blocked_shortest_cycles_100_1800_10000"
# ]

# organism = "iML1515"
# file_names = [
#     "iML1515_loopless_fba_1800", 
#     "iML1515_loopless_fba_nullspace_1800", 
#     "iML1515_loopless_indicator_fba_1800", 
#     "iML1515_loopless_indicator_fba_blocked_500_1800_50_same_objective", #TODO: different objective
#     "iML1515_loopless_indicator_fba_blocked_500_1800_100_same_objective", 
#     "iML1515_loopless_fba_blocked_100_1800_50_same_objective", 
#     "iML1515_loopless_fba_blocked_100_1800_100_same_objective", 
#     "iML1515_loopless_fba_blocked_100_1800_200_same_objective", 
#     "iML1515_loopless_fba_blocked_100_1800_500_same_objective",
#     "iML1515_loopless_fba_blocked_100_1800_50", 
#     "iML1515_loopless_fba_blocked_100_1800_100", 
#     "iML1515_loopless_fba_blocked_100_1800_200", 
#     "iML1515_loopless_fba_blocked_100_1800_500", 
#     #"iML1515_loopless_fba_blocked_50_1800_10000", 
#     #"iML1515_loopless_fba_blocked_100_1800_10000", 
#     #"iML1515_loopless_fba_blocked_shortest_cycles_50_1800_10000", 
#     "iML1515_loopless_fba_blocked_shortest_cycles_100_1800_50"
# ]

# compare_blocked_cycles(file_names=file_names, organism=organism)

# organism = "iJR904"
# file_names = [
#     "iJR904_loopless_fba_1800", 
#     "iJR904_loopless_fba_nullspace_1800", 
#     "iJR904_loopless_indicator_fba_1800", 
#     "iJR904_loopless_indicator_fba_blocked_500_1800_50_same_objective", #TODO: different objective
#     "iJR904_loopless_indicator_fba_blocked_500_1800_100_same_objective", 
#     "iJR904_loopless_fba_blocked_100_1800_50_same_objective", 
#     "iJR904_loopless_fba_blocked_100_1800_100_same_objective", 
#     "iJR904_loopless_fba_blocked_100_1800_200_same_objective", 
#     "iJR904_loopless_fba_blocked_100_1800_500_same_objective",
#     "iJR904_loopless_fba_blocked_100_1800_50", 
#     "iJR904_loopless_fba_blocked_100_1800_100", 
#     "iJR904_loopless_fba_blocked_100_1800_200", 
#     "iJR904_loopless_fba_blocked_100_1800_500", 
#     "iJR904_loopless_fba_blocked_50_1800_10000", 
#     "iJR904_loopless_fba_blocked_100_1800_10000", 
#     "iJR904_loopless_fba_blocked_shortest_cycles_50_1800_10000", 
#     "iJR904_loopless_fba_blocked_shortest_cycles_100_1800_50"
# ]

# compare_blocked_cycles(file_names=file_names, organism=organism)