using DataFrames, JSON, CSV
using Statistics

function sub_csv(file_out, file_in="results_bigg_SCIP.csv"; nullspace=false, ll_fba=false, cb=false, cb_big_m=false, cb_indicator_and_big_m=false, objective_values=true, mis_indicator=false, mis_big_m=false, mis_numbers=[], no_good_cuts=false, no_good_cuts_big_m=false, ch=false, ch_mis=false, ll_fba_indicator=false)
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
            end 
            if mis_big_m
                append!(cols, ["termination_cb_big_m_mis_" * string(mis), "time_cb_big_m_mis_" * string(mis)])
                if objective_values
                    append!(cols, ["objective_value_cb_big_m_mis_" * string(mis)])
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

function sub_csv_gecko(file_out, file_in="results_bigg_SCIP_gecko_1.0e-8.csv"; cb_big_m=true, objective_values=true, iter=false, mis_big_m=true, mis_numbers=[0.5], ll_fba=true)
    df = CSV.read("csv/" * file_in, DataFrame)

    cols = ["organism", "mean"]

    if ll_fba 
        append!(cols, ["termination_loopless_fba_gecko_tol_1.0e-8", "time_loopless_fba_gecko_tol_1.0e-8"])
        if objective_values
            append!(cols, ["objective_value_loopless_fba_gecko_tol_1.0e-8"])
        end
    end 
    
    if cb_big_m 
        append!(cols, ["termination_combinatorial_benders_gecko_fast_big_m_tol_1.0e-8", "time_combinatorial_benders_gecko_fast_big_m_tol_1.0e-8"])
        if objective_values
            append!(cols, ["objective_value_combinatorial_benders_gecko_fast_big_m_tol_1.0e-8"])
        end
        if iter 
            append!(cols, ["iter_combinatorial_benders_gecko_fast_big_m_tol_1.0e-8"])
        end
    end 
    
    if !isempty(mis_numbers)
        for mis in mis_numbers
            if mis_big_m
                append!(cols, ["termination_combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_1.0e-8", "time_combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_1.0e-8"])
                if objective_values
                    append!(cols, ["objective_value_combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_1.0e-8"])
                end
                if iter 
                    append!(cols, ["iter_combinatorial_benders_gecko_fast_big_m_" * string(mis) * "_tol_1.0e-8"])
                end
            end 
        end 
    end

    # get model size 
    df_size = CSV.read("../molecular_models/bigg_model_data.csv", DataFrame)
    df = innerjoin(df, df_size, on="organism")
    sort!(df, [:reactions])
    df = df[!, cols]

    CSV.write("csv/" * file_out, df, append=false, writeheader=true)
    return file_out
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
    "iSynCJ816",
    "iYO844",
    "iYS1720",
    "iZ_1308"
]

# sub_csv("results_ll_fba_variants.csv", nullspace=true, ll_fba=true)
# sub_csv("results_cb_indicator_vs_big_m.csv", cb=true, cb_big_m=true, cb_indicator_and_big_m=true)

sub_csv_gecko("results_ll_fba_vs_cb_gecko.csv", iter=true)
