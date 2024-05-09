using DataFrames, JSON, CSV
using Statistics

"""
generate filenames or column names
"""
function generate_setups(
    enzyme_data,
    tol, 
    ll_fba,
    cb_big_m, 
    mis_big_m, 
    mis_numbers;
    set_seed=false,
    seed=1,
    mean=1.0
)

    setups = String[]
    if enzyme_data
        if set_seed 
            if mean == 1.0
                if ll_fba
                    push!(setups, "loopless_fba_gecko_" * string(seed) * "_tol_" * string(tol))
                end
                if cb_big_m 
                    push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_fast_big_m_tol_" * string(tol))
                    if mis_big_m
                        for mis in mis_numbers
                            push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_fast_big_m_" * string(mis) * "_mis_tol_" * string(tol))
                        end
                    end
                end 
            else 
                if ll_fba
                    push!(setups, "loopless_fba_gecko_" * string(seed) * "_" * string(mean) * "_tol_" * string(tol))
                end
                if cb_big_m
                    push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_" * string(mean) * "_fast_big_m_tol_" * string(tol))
                    if mis_big_m
                        for mis in mis_numbers
                            push!(setups, "combinatorial_benders_gecko_" * string(seed) * "_" * string(mean) * "_fast_big_m_" * string(mis) * "_mis_tol_" * string(tol))
                        end
                    end
                end
            end
        else 
            if ll_fba
                push!(setups, "loopless_fba_gecko_tol_" * string(tol))
            end
            if cb_big_m 
                push!(setups, "combinatorial_benders_gecko_fast_big_m_tol_" * string(tol))
                if mis_big_m
                    for mis in mis_numbers
                        push!(setups, "combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_" * string(tol))
                    end
                end
            end
        end
    end

    return setups
end 

# read and clean up instance
function build_dict(
    organism,
    file_name,
    setup,
    enzyme_data;
    termination=true, 
    objective_value=true, 
    time=true,
    feasibility=true,
    time_limit,
    )

    cleaned_up_dict = Dict{Symbol, Any}()
    dict = Dict()
    try 
        if enzyme_data
            dict = JSON.parse(open("json/" * organism * "_" * file_name * "_" * string(time_limit) * ".json"))
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
    enzyme_data=true,
    seeds=[10],
    tol=1.0e-8, 
    ll_fba=false,
    cb_big_m=false, 
    μ_mean=[], 
    mis_big_m=false, 
    mis_numbers=[]
    )
    
    df = DataFrame(
        organism = Union{String,Missing}[], 
        time_limit= Union{Int64,Missing}[],
        seed = Int64[],
        mean = Float64[],
        tolerance = Float64[]
    )

    setups = generate_setups(
        enzyme_data,
        tol, 
        ll_fba,
        cb_big_m, 
        mis_big_m, 
        mis_numbers
    )
    @show setups

    for setup in setups
        df[!, "termination_" * setup] = Union{String,Missing}[] 
        df[!, "objective_value_" * setup] = Union{Float64,Missing}[] 
        df[!, "time_" * setup] = Union{Float64,Missing}[]
        df[!, "feasibility_" * setup] = Union{Bool,Missing}[]
        df[!, "iter_" * setup] = Union{Int64,Missing}[]
    end

    for organism in organisms
        for seed in seeds 
            dict_organism = Dict{Symbol, Any}()
            dict_organism[:organism] = organism
            dict_organism[:time_limit] = time_limit
            dict_organism[:seed] = seed
            dict_organism[:mean] = 1.0 
            dict_organism[:tolerance] = tol 
            
            # fill df
            file_names = generate_setups(
                enzyme_data,
                tol, 
                ll_fba,
                cb_big_m, 
                mis_big_m, 
                mis_numbers;
                seed=seed,
                mean=1.0,
                set_seed=true
            )

            for (idx, file_name) in enumerate(file_names)
                dict = build_dict(
                    organism,
                    file_name,
                    setups[idx],
                    enzyme_data;
                    termination=true, 
                    objective_value=true, 
                    time=true,
                    feasibility=true,
                    time_limit=time_limit
                )
                merge!(dict_organism, dict)
            end 
            push!(df, dict_organism)
            for mean in μ_mean
                dict_organism = Dict{Symbol, Any}()
                dict_organism[:organism] = organism
                dict_organism[:time_limit] = time_limit
                dict_organism[:seed] = seed
                dict_organism[:mean] = mean 
                dict_organism[:tolerance] = tol 
                
                # fill df
                file_names = generate_setups(
                    enzyme_data,
                    tol, 
                    ll_fba,
                    cb_big_m, 
                    mis_big_m, 
                    mis_numbers;
                    seed=seed,
                    mean=mean,
                    set_seed=true
                )

                for (idx, file_name) in enumerate(file_names)
                    dict = build_dict(
                        organism,
                        file_name,
                        setups[idx],
                        enzyme_data;
                        termination=true, 
                        objective_value=true, 
                        time=true,
                        feasibility=true,
                        time_limit=time_limit
                    )
                    merge!(dict_organism, dict)
                end 
                push!(df, dict_organism)
            end
        end
    end
    # println(df)
    CSV.write("csv/" * save_as, df, append=false, writeheader=true)
end

organisms = [
    "iAF692", 
    # "iJR904", 
    "iML1515", 
    "e_coli_core",
    "iNF517",
    "iSB619",
    #  "iNJ661",
    # "iCN900",   
    "iAF1260",
    "iEK1008",
    "iJO1366",
    "iMM904",
    "iSDY_1059",
    "iSFV_1184",
    "iSF_1195",
    "iS_1188",
    "iSbBS512_1146",
    # "RECON1",
    "Recon3D",
    "STM_v1_0",
    # "iAB_RBC_283",
    "iAPECO1_1312",
    "iECB_1328",
    "iETEC_1333",
    "iHN637",
    # "iIS312_Amastigote",
    "iJB785",
    # "iJN746",
    "iLB1027_lipid",
    # "iMM1415",
    "iND750",
    # "iRC1080",
    "iSFxv_1172",
    # "iSynCJ816",
    "iYO844",
    "iYS1720",
    "iZ_1308"
]

build_csv(
    "results_bigg_SCIP_gecko_1.0e-8.csv", 
    organisms,
    time_limit=1800,
    enzyme_data=true,
    seeds=[10],
    tol=1.0e-8, 
    ll_fba=true,
    cb_big_m=true, 
    μ_mean=[1.2, 1.5], 
    mis_big_m=true, 
    mis_numbers=[0.5]
)

# solver_data(organisms, save_as="results_bigg_gecko_1.0e-8.csv", time_limit=1800, enzyme_data=true, seeds=[10], tol=1.0e-8, cb_big_m=true, μ_mean=[1.2, 1.5], mis_big_m=true, mis_numbers=[0.5]) 
# mis_indicator_and_big_m=true,