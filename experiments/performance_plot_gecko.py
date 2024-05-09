import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(filename='csv/results_bigg_SCIP.csv', indicator_and_big_m=False):
    # load data 
    df_data = pd.read_csv(filename)
    # print(df_data)
    data = {}

    # ll fba data
    time_ll_fba = df_data[df_data["termination_loopless_fba"] == "OPTIMAL"]["time_loopless_fba"].to_list()
    time_ll_fba.sort()
    time_ll_fba.append(1800)
    print(time_ll_fba)
    instances_ll_fba = len(time_ll_fba) + 1
    instances_ll_fba = np.arange(1, instances_ll_fba)
    data["time_ll_fba"] = time_ll_fba
    data["instances_ll_fba"] = instances_ll_fba

    # # cb data indicator
    # time_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["time_cb"].to_list()
    # time_cb.sort()
    # time_cb.append(1800)
    # instances_cb = len(time_cb) + 1
    # instances_cb = np.arange(1, instances_cb)
    # iter_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["iter_cb"].to_list()
    # iter_cb.sort()
    # # time_mp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_master_problem_cb"].to_list()
    # # time_mp_cb.sort()
    # # time_sp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_sub_problem_cb"].to_list()
    # # time_sp_cb.sort()
    # # time_mis_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_mis_problem_cb"].to_list()
    # # time_mis_cb.sort()
    # data["time_cb"] = time_cb
    # data["instances_cb"] = instances_cb
    # data["iter_cb"] = iter_cb
    # # data["time_mp_cb"] = time_mp_cb
    # # data["time_sp_cb"] = time_sp_cb
    # # data["time_mis_cb"] = time_mis_cb

    # cb data big m
    df_cb_big_m = df_data[df_data["termination_combinatorial_benders_big_m"] == "OPTIMAL"]
    df_cb_big_m = df_cb_big_m[df_cb_big_m["time_combinatorial_benders_big_m"] <= 1800]
    time_cb_big_m = df_cb_big_m["time_combinatorial_benders_big_m"].to_list()
    time_cb_big_m.sort()
    time_cb_big_m.append(1800)
    instances_cb_big_m = len(time_cb_big_m) + 1
    instances_cb_big_m = np.arange(1, instances_cb_big_m)
    iter_cb_big_m = df_cb_big_m["iter_combinatorial_benders_big_m"].to_list()
    iter_cb_big_m.sort()
    # time_mp_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_master_problem_cb"].to_list()
    # time_mp_cb.sort()
    # time_sp_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_sub_problem_cb"].to_list()
    # time_sp_cb.sort()
    # time_mis_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_mis_problem_cb"].to_list()
    # time_mis_cb.sort()
    data["time_cb_big_m"] = time_cb_big_m
    data["instances_cb_big_m"] = instances_cb_big_m
    data["iter_cb_big_m"] = iter_cb_big_m
    # data["time_mp_cb_big_m"] = time_mp_cb
    # data["time_sp_cb_big_m"] = time_sp_cb
    # data["time_mis_cb_big_m"] = time_mis_cb

    # cb big m mis 0.5
    df_cb_big_m_mis = df_data[df_data["termination_combinatorial_benders_big_m_0_5"] == "OPTIMAL"]
    df_cb_big_m_mis = df_cb_big_m_mis[df_cb_big_m_mis["time_combinatorial_benders_big_m_0_5"] <= 1800]
    time_cb_big_m_mis = df_cb_big_m_mis["time_combinatorial_benders_big_m_0_5"].to_list()
    time_cb_big_m_mis.sort()
    time_cb_big_m_mis.append(1800)
    instances_cb_big_m_mis = len(time_cb_big_m_mis) + 1
    instances_cb_big_m_mis = np.arange(1, instances_cb_big_m_mis)
    iter_cb_big_m_mis = df_cb_big_m_mis["iter_combinatorial_benders_big_m_0_5"].to_list()
    iter_cb_big_m_mis.sort()
    data["time_cb_big_m_mis"] = time_cb_big_m_mis
    data["instances_cb_big_m_mis"] = instances_cb_big_m_mis
    data["iter_cb_big_m_mis"] = iter_cb_big_m_mis

    if indicator_and_big_m:
        df_cb = df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]
        df_cb = df_cb[df_cb["time_cb_indicator_and_big_m"] <= 1800]
        time_cb = df_cb["time_cb_indicator_and_big_m"].to_list()
        time_cb.sort()
        time_cb.append(1800)
        instances_cb = len(time_cb) + 1
        instances_cb = np.arange(1, instances_cb)
        iter_cb = df_cb["iter_cb_indicator_and_big_m"].to_list()
        iter_cb.sort()
        data["time_cb_indicator_and_big_m"] = time_cb
        data["instances_cb_indicator_and_big_m"] = instances_cb
        data["iter_cb_indicator_and_big_m"] = iter_cb

    return data 


def build_solved_instances_plot(filename, colors, linestyles, save_as="", indicator_and_big_m=False):
    data = load_data(filename, indicator_and_big_m=indicator_and_big_m)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })

    print(data.keys())
    # create solved instances plot
    fig, axs = plt.subplots(2, 1, figsize=[6.8, 4.8]) #, layout='constrained')

    axs[0].plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], label="ll-FBA (big-M)", linestyle=linestyles[0])
    axs[0].plot(data["time_cb_big_m"], data["instances_cb_big_m"], color=colors[1], linestyle=linestyles[1])
    axs[0].plot(data["time_cb_big_m_mis"], data["instances_cb_big_m_mis"], color=colors[2], linestyle=linestyles[2])
    # if indicator_and_big_m:
    #     axs[0].plot(data["time_cb_indicator_and_big_m"], data["instances_cb_indicator_and_big_m"], color=colors[3], linestyle=linestyles[3])
    axs[0].set_ylabel(r"$\textbf{solved instances}$")
    axs[0].set_xlabel(r"$\textbf{time (s)}$")
    axs[0].grid(True)

    axs[1].plot(data["iter_cb_big_m"], np.arange(1, len(data["iter_cb_big_m"])+1), label="CB (big-M)", color=colors[1], linestyle=linestyles[1])
    axs[1].plot(data["iter_cb_big_m_mis"], np.arange(1, len(data["iter_cb_big_m_mis"])+1), label="CB (big-M MIS 0.5\%)", color=colors[2], linestyle=linestyles[2])
    # if indicator_and_big_m:
    #         axs[1].plot(data["iter_cb_indicator_and_big_m"], np.arange(1, len(data["iter_cb_indicator_and_big_m"])+1), label="CB (indicator + big M)", color=colors[3], linestyle=linestyles[3])
    axs[1].set_ylabel(r"$\textbf{solved instances}$")
    axs[1].set_xlabel(r"$\textbf{iterations}$")
    axs[1].grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.05), ncol=3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    # plt.show()
    if save_as == "":
        plt.savefig("plots/comparison_solved_instances_gecko.pdf", format="pdf", bbox_inches="tight")
    else:
        plt.savefig(save_as, format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
linestyles=["-.", "--", ":", "-", "--", "-"]

# filename = 'csv/results_bigg_SCIP_gecko.csv'
# build_solved_instances_plot(filename, colors, linestyles, save_as="plots/comparison_solved_instances_gecko.pdf")

filename = 'csv/results_ll_fba_vs_cb_gecko_filtered.csv'
build_solved_instances_plot(filename, colors, linestyles, save_as="plots/comparison_solved_instances_gecko_1.0e-8.pdf")

# filename = 'csv/results_bigg_gecko_indicator_and_big_m_1.0e-8.csv'
# build_solved_instances_plot(filename, colors, linestyles, indicator_and_big_m=True, save_as="plots/comparison_solved_instances_gecko_1.0e-8_indicator_and_big_m.pdf")
