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

function compare_loopless_formulation_blocked_cycles(; file_names, organisms=["iAF692","iJR904","iML1515"])
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

"""
compare ll-FBA, no good cuts and fast combinatorial Benders
"""
function loopless_fba_vs_cb(organisms; cuts=true)
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
    
    for organism in organisms
        @show organism
        dict_organism = Dict{Symbol, Any}()
        dict_organism[:organism] = organism

        # read cobrexa data
        df_temp = first(CSV.read(organism * "/server/" * organism * "_cobrexa_loopless_fba_1800.csv", DataFrame),1)
        dict_organism[:termination_ll_fba_cobrexa] = df_temp[!, [:termination][1]][1]
        dict_organism[:objective_value_ll_fba_cobrexa] = df_temp[!, [:objective_value][1]][1]
        dict_organism[:time_ll_fba_cobrexa] = df_temp[!, [:time][1]][1]

        # read ll-FBA data
        df_temp = first(CSV.read(organism * "/server/" * organism * "_loopless_fba_nullspace_1800.csv", DataFrame),1)
        dict_organism[:time_limit] = df_temp[!, [:time_limit][1]][1]
        dict_organism[:termination_ll_fba_nullspace] = df_temp[!, [:termination][1]][1]
        dict_organism[:objective_value_ll_fba_nullspace] = df_temp[!, [:objective_value][1]][1]
        dict_organism[:time_ll_fba_nullspace] = df_temp[!, [:time][1]][1]

        # read ll-FBA without nullspace formulation data
        df_temp = first(CSV.read(organism * "/server/" * organism * "_loopless_fba_1800.csv", DataFrame),1)
        dict_organism[:termination_ll_fba] = df_temp[!, [:termination][1]][1]
        dict_organism[:objective_value_ll_fba] = df_temp[!, [:objective_value][1]][1]
        dict_organism[:time_ll_fba] = df_temp[!, [:time][1]][1]

        # read no good cuts data 
        if organism == "iAF692"
            df_temp = first(CSV.read(organism * "/server/" * organism * "_combinatorial_benders_1800_1e-5.csv", DataFrame),1)
        else 
            df_temp = first(CSV.read(organism * "/server/" * organism * "_combinatorial_benders_1800.csv", DataFrame),1)
        end
        # set termination status to TIME_LIMIT, INFEASIBLE, OPTIMAL
        if df_temp[!, [:termination][1]][1] == "INFEASIBLE"
            if df_temp[!, [:time][1]][1] >=  df_temp[!, [:time_limit][1]][1]
                dict_organism[:termination_no_good_cuts] = "TIME_LIMIT"
            else 
                dict_organism[:termination_no_good_cuts] = "INFEASIBLE"
            end
        elseif df_temp[!, [:thermo_feasible][1]][1] == true
            dict_organism[:termination_no_good_cuts] = "OPTIMAL"
        end
        dict_organism[:objective_value_no_good_cuts] = df_temp[!, [:objective_value][1]][1]
        dict_organism[:time_no_good_cuts] = df_temp[!, [:time][1]][1]

        # read fast CB data
        df_temp = first(CSV.read(organism * "/server/" * organism * "_combinatorial_benders_fast_1800.csv", DataFrame),1)
        if df_temp[!, [:termination][1]][1] == "INFEASIBLE"
            if df_temp[!, [:time][1]][1] >=  df_temp[!, [:time_limit][1]][1]
                dict_organism[:termination_cb] = "TIME_LIMIT"
            else 
                dict_organism[:termination_cb] = "INFEASIBLE"
            end
        elseif df_temp[!, [:thermo_feasible][1]][1] == true
            dict_organism[:termination_cb] = "OPTIMAL"
        end
        dict_organism[:objective_value_cb] = df_temp[!, [:objective_value][1]][1]
        dict_organism[:time_cb] = df_temp[!, [:time][1]][1]

        push!(df, dict_organism)
    end
    
    # round objective value, cut time above time limit and convert to int
    for row in eachrow(df)
        if row.time_ll_fba > 1800
            row.time_ll_fba = 1800
        end
        if row.time_ll_fba_cobrexa > 1800
            row.time_ll_fba_cobrexa = 1800
        end
        if row.time_ll_fba_nullspace > 1800
            row.time_ll_fba_nullspace = 1800
        end
        if row.time_no_good_cuts > 1800
            row.time_no_good_cuts = 1800
        end
        if row.time_cb > 1800
            row.time_cb = 1800
        end
    end

    df[!, :objective_value_ll_fba_cobrexa] = round.(df[!, :objective_value_ll_fba_cobrexa], digits=3)
    df[!, :objective_value_ll_fba] = round.(df[!, :objective_value_ll_fba], digits=3)
    df[!, :objective_value_ll_fba_nullspace] = round.(df[!, :objective_value_ll_fba_nullspace], digits=3)
    df[!, :objective_value_no_good_cuts] = round.(df[!, :objective_value_no_good_cuts], digits=3)
    df[!, :objective_value_cb] = round.(df[!, :objective_value_cb], digits=3)
    df[!, :time_ll_fba_cobrexa] = round.(df[!, :time_ll_fba_cobrexa], digits=0)
    df[!, :time_ll_fba_cobrexa] = Int.(df[!, :time_ll_fba_cobrexa])
    df[!, :time_ll_fba] = round.(df[!, :time_ll_fba], digits=0)
    df[!, :time_ll_fba] = Int.(df[!, :time_ll_fba])
    df[!, :time_ll_fba_nullspace] = round.(df[!, :time_ll_fba_nullspace], digits=0)
    df[!, :time_ll_fba_nullspace] = Int.(df[!, :time_ll_fba_nullspace])
    df[!, :time_no_good_cuts] = round.(df[!, :time_no_good_cuts], digits=0)
    df[!, :time_no_good_cuts] = Int.(df[!, :time_no_good_cuts])
    df[!, :time_cb] = round.(df[!, :time_cb], digits=0)
    df[!, :time_cb] = Int.(df[!, :time_cb])

    @show df[!, [:objective_value_ll_fba_nullspace, :objective_value_ll_fba, :objective_value_no_good_cuts, :objective_value_cb]]
    @show df[!, [:termination_ll_fba_nullspace, :termination_ll_fba, :termination_no_good_cuts, :termination_cb]]
    @show df[!, [:time_ll_fba_nullspace, :time_ll_fba, :time_no_good_cuts, :time_cb]]
    
    if cuts 
        file_name = "comparison_ll_fba_vs_cb.csv"
        CSV.write(file_name, df, append=false, writeheader=true)
    else 
        df = df[!, [:time_ll_fba, :objective_value_ll_fba, :termination_ll_fba, :time_ll_fba_nullspace, :objective_value_ll_fba_nullspace, :termination_ll_fba_nullspace, :time_ll_fba_cobrexa, :objective_value_ll_fba_cobrexa, :termination_ll_fba_cobrexa]]
        file_name = "comparison_ll_fba.csv"
        CSV.write(file_name, df, append=false, writeheader=true)
    end
end

organisms = ["iAF692", "e_coli_core", "iJR904", "iML1515", "iNF517", "iNJ661", "iCN900"] # "iSB619" not feasible
loopless_fba_vs_cb(organisms, cuts=false)


# compare_loopless_formulation(file_names=["loopless_fba_1800", "loopless_fba_nullspace_1800", "loopless_indicator_fba_1800"])

# organism = "iAF692"
# file_names = [
#     "iAF692_loopless_fba_1800", 
#     "iAF692_loopless_fba_nullspace_1800", 
#     "iAF692_loopless_indicator_fba_1800", 
#     "iAF692_loopless_indicator_fba_blocked_500_1800_50_same_objective", 
#     "iAF692_loopless_indicator_fba_blocked_500_1800_100_same_objective",
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
#     "iML1515_loopless_indicator_fba_blocked_500_1800_50_same_objective",
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
#     "iJR904_loopless_indicator_fba_blocked_500_1800_50_same_objective",
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