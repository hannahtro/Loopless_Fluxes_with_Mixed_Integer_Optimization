using DataFrames, JSON, CSV
using Statistics

"""
generate filenames or column names
"""
function generate_setups(
    ll_fba,
    block_limits,        
    cycles_blocked,
    shortest_cycles_blocked,
    seeds
)

    setups = String[]
    if ll_fba
        push!(setups, "loopless_fba")
        for seed in seeds 
            if seed != 0
                push!(setups, "loopless_fba_1800_scip_seed_" * string(seed))
            end
        end
    end
    if cycles_blocked
        for limit in block_limits
            if !isempty(seeds)
                for seed in seeds
                    if seed == 0
                        push!(setups, "loopless_fba_blocked_" * string(limit))
                    else 
                        push!(setups, "loopless_fba_blocked_" * string(limit) * "_1800_scip_seed_" * string(seed))
                    end
                end 
            else
                push!(setups, "loopless_fba_blocked_" * string(limit))
            end
        end
    end 
    if shortest_cycles_blocked
        for limit in block_limits
            if !isempty(seeds)
                for seed in seeds
                    if seed == 0
                        push!(setups, "loopless_fba_blocked_shortest_cycles_" * string(limit))
                    else 
                        push!(setups, "loopless_fba_blocked_shortest_cycles_" * string(limit) * "_1800_scip_seed_" * string(seed))
                    end
                end
            else
                push!(setups, "loopless_fba_blocked_shortest_cycles_" * string(limit))
            end
        end
    end

    # if enzyme_data
    #     if set_seed 
    #         if mean == 1.0
    #             if ll_fba
    #                 push!(setups, "loopless_fba_gecko_" * string(seed) * "_tol_" * string(tol))
    #             end
    #             if cb_big_m 
    #                 push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_fast_big_m_tol_" * string(tol))
    #                 if mis_big_m
    #                     for mis in mis_numbers
    #                         push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_fast_big_m_" * string(mis) * "_mis_tol_" * string(tol))
    #                     end
    #                 end
    #             end 
    #         else 
    #             if ll_fba
    #                 push!(setups, "loopless_fba_gecko_" * string(seed) * "_" * string(mean) * "_tol_" * string(tol))
    #             end
    #             if cb_big_m
    #                 push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_" * string(mean) * "_fast_big_m_tol_" * string(tol))
    #                 if mis_big_m
    #                     for mis in mis_numbers
    #                         push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_" * string(mean) * "_fast_big_m_" * string(mis) * "_mis_tol_" * string(tol))
    #                     end
    #                 end
    #             end
    #         end
    #     else 
    #         if ll_fba
    #             push!(setups, "loopless_fba_gecko_tol_" * string(tol))
    #         end
    #         if cb_big_m 
    #             push!(setups, "combinatorial_benders_gecko_fast_big_m_tol_" * string(tol))
    #             if mis_big_m
    #                 for mis in mis_numbers
    #                     push!(setups, "combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_" * string(tol))
    #                 end
    #             end
    #         end
    #     end
    # end

    return setups
end 

# read and clean up instance
function build_dict(
    organism,
    file_name,
    setup;
    num_cycles=true,
    termination=true, 
    objective_value=true, 
    time=true,
    feasibility=true,
    time_limit,
    )

    cleaned_up_dict = Dict{Symbol, Any}()
    dict = Dict()
    try 
        if occursin("blocked", file_name)
            if occursin("seed", file_name)
                dict = JSON.parse(open("json/" * organism * "_" * file_name * "_1000.json"))
            else
                dict = JSON.parse(open("json/" * organism * "_" * file_name * "_" * string(time_limit) * "_1000.json"))
            end
        else 
            if occursin("seed", file_name)
                dict = JSON.parse(open("json/" * organism * "_" * file_name * ".json"))
            else
                dict = JSON.parse(open("json/" * organism * "_" * file_name * "_" * string(time_limit) * ".json"))
            end
        end
    catch e 
        println(e)
        dict = Dict{String, Any}(
            "termination" => "ERROR",
            "objective_value" => missing,
            "time" => NaN, 
            "thermo_feasible" => false,
            "iter" => missing,
            # "cuts" => missing
            num_cycles => missing
        )
    end 

    if dict["termination"] == "INFEASIBLE" || dict["termination"] == "TIME_LIMIT"
        if dict["time"] >= dict["time_limit"]
            cleaned_up_dict[Symbol("termination_" * setup)] = "TIME_LIMIT"
        else 
            cleaned_up_dict[Symbol("termination_" * setup)] = "INFEASIBLE"
        end
    elseif dict["thermo_feasible"] == true
        cleaned_up_dict[Symbol("termination_" * setup)] = "OPTIMAL"
    else 
        cleaned_up_dict[Symbol("termination_" * setup)] = dict["termination"]
    end

    cleaned_up_dict[Symbol("objective_value_" * setup)] = dict["objective_value"]
    if isnothing(cleaned_up_dict[Symbol("objective_value_" * setup)])
        cleaned_up_dict[Symbol("objective_value_" * setup)] = missing
    end
    cleaned_up_dict[Symbol("time_" * setup)] = dict["time"]
    cleaned_up_dict[Symbol("feasibility_" * setup)] = dict["thermo_feasible"]
    if occursin("blocked", file_name)
        cleaned_up_dict[Symbol("num_cycles_" * setup)] = dict["num_blocked_cycles"]
    else 
        cleaned_up_dict[Symbol("num_cycles_" * setup)] = missing
    end

    try 
        cleaned_up_dict[Symbol("iter_" * setup)] = dict["iter"]
    catch e 
        cleaned_up_dict[Symbol("iter_" * setup)] = missing
    end 

    # cleaned_up_dict[:cuts_cb] = dict["cuts"]
    # cleaned_up_dict[Symbol("times_master_problem_cb")] = mean(dict["times_master_problem"])
    # cleaned_up_dict[Symbol("times_sub_problem_cb")] = mean(dict["times_sub_problem"])
    # cleaned_up_dict[Symbol("times_mis_problem_cb")] = mean(dict["times_mis_problem"])

    return cleaned_up_dict
end 

"""
build CSV to compare the results of different setups
"""
function build_csv(
    save_as, 
    organisms;
    time_limit=1800,
    block_limits = [10,20,50,100],
    seeds=[],
    ll_fba=true,
    cycles_blocked=true,
    shortest_cycles_blocked=true
    )
    
    df = DataFrame(
        organism = Union{String,Missing}[], 
        time_limit= Union{Int64,Missing}[],
    )

    setups = generate_setups(
        ll_fba,
        block_limits,        
        cycles_blocked,
        shortest_cycles_blocked,
        seeds
    )
    # @show setups

    for setup in setups
        df[!, "termination_" * setup] = Union{String,Missing}[] 
        df[!, "objective_value_" * setup] = Union{Float64,Missing}[] 
        df[!, "time_" * setup] = Union{Float64,Missing}[]
        df[!, "feasibility_" * setup] = Union{Bool,Missing}[]
        df[!, "iter_" * setup] = Union{Int64,Missing}[]
        df[!, "num_cycles_" * setup] = Union{Int64,Missing}[]
    end

    for organism in organisms
        dict_organism = Dict{Symbol, Any}()
        dict_organism[:organism] = organism
        dict_organism[:time_limit] = time_limit 
        
        # fill df
        file_names = generate_setups(
            ll_fba,
            block_limits,        
            cycles_blocked,
            shortest_cycles_blocked,
            seeds
        )

        for (idx, file_name) in enumerate(file_names)
            dict = build_dict(
                organism,
                file_name,
                setups[idx];
                termination=true, 
                objective_value=true, 
                time=true,
                feasibility=true,
                time_limit=time_limit,
                num_cycles=true
            )
            merge!(dict_organism, dict)
        end 
        push!(df, dict_organism)
    end
    # println(df)
    CSV.write("csv/" * save_as, df, append=false, writeheader=true)
end

function build_sub_csv(
    input_file="results_cff.csv",
    output_file="results_cff_reduced.csv",
    characteristics=["time", "num_cycles"]
)
    df = CSV.read("csv/" * input_file, DataFrame)

    columns_to_select = []
    for char in characteristics
        for col in names(df)
            if occursin(char, col)
                push!(columns_to_select, col)
            end
        end
    end

    for col in columns_to_select
        df[!, col] = round.(df[!, col])
    end

    push!(columns_to_select, "organism")
    select!(df, columns_to_select)
    select!(df, Not("time_limit"))

    CSV.write("csv/" * output_file, df, append=false, writeheader=true)
end

organisms = [
    "iAF692", 
    "iJR904", 
    "iML1515", 
    "e_coli_core",
    # "iNF517",
    "iNJ661",
    "iCN900",
    "iEK1008",
    "iMM904",
    "iSFV_1184",
]

build_csv(
    "results_cff.csv", 
    organisms,
    time_limit=1800,
    block_limits = [10,20,50,100],
    ll_fba=true,
    cycles_blocked=true,
    shortest_cycles_blocked=true,
    seeds=0:4
)

# build_sub_csv()
