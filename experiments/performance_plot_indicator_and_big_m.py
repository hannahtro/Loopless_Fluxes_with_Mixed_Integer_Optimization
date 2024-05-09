import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_SCIP.csv', big_m=False, indicator=True, indicator_and_big_m=False, ll_fba=False, ll_fba_indicator=False, extended=False):
    # load data 
    df_data = pd.read_csv(file)
    # print(df_data)
    data = {}

    # ll fba data
    if ll_fba:
        time_ll_fba = df_data[df_data["termination_ll_fba"] == "OPTIMAL"]["time_ll_fba"].to_list()
        time_ll_fba.sort()
        time_ll_fba.append(1800)
        # print(time_ll_fba)
        instances_ll_fba = len(time_ll_fba) + 1
        instances_ll_fba = np.arange(1, instances_ll_fba)
        data["time_ll_fba"] = time_ll_fba
        data["instances_ll_fba"] = instances_ll_fba

    # ll fba indicator data
    if ll_fba_indicator:
        time_ll_fba = df_data[df_data["termination_ll_fba_indicator"] == "OPTIMAL"]["time_ll_fba_indicator"].to_list()
        time_ll_fba.sort()
        time_ll_fba.append(1800)
        # print(time_ll_fba)
        instances_ll_fba = len(time_ll_fba) + 1
        instances_ll_fba = np.arange(1, instances_ll_fba)
        data["time_ll_fba_indicator"] = time_ll_fba
        data["instances_ll_fba_indicator"] = instances_ll_fba

    # cb data
    if indicator:
        time_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["time_cb"].to_list()
        time_cb.sort()
        time_cb.append(1800)
        instances_cb = len(time_cb) + 1
        instances_cb = np.arange(1, instances_cb)
        data["time_cb"] = time_cb
        data["instances_cb"] = instances_cb

        if extended:
            iter_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["iter_cb"].to_list()
            iter_cb.sort()
            data["iter_cb"] = iter_cb
            time_mp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_master_problem_cb"].to_list()
            time_mp_cb.sort()
            time_sp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_sub_problem_cb"].to_list()
            time_sp_cb.sort()
            time_mis_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_mis_problem_cb"].to_list()
            time_mis_cb.sort()
            data["time_mp_cb"] = time_mp_cb
            data["time_sp_cb"] = time_sp_cb
            data["time_mis_cb"] = time_mis_cb


    if big_m:
        time_cb_big_m = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["time_cb_big_m"].to_list()
        time_cb_big_m.sort()
        time_cb_big_m.append(1800)
        instances_cb_big_m = len(time_cb_big_m) + 1
        instances_cb_big_m = np.arange(1, instances_cb_big_m)
        data["time_cb_big_m"] = time_cb_big_m
        data["instances_cb_big_m"] = instances_cb_big_m
        
        if extended:
            iter_cb_big_m = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["iter_cb_big_m"].to_list()
            iter_cb_big_m.sort()
            time_mp_cb_big_m= df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_master_problem_cb_big_m"].to_list()
            time_mp_cb_big_m.sort()
            time_sp_cb_big_m = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_sub_problem_cb_big_m"].to_list()
            time_sp_cb_big_m.sort()
            time_mis_cb_big_m = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_mis_problem_cb_big_m"].to_list()
            time_mis_cb_big_m.sort()
            data["iter_cb_big_m"] = iter_cb_big_m
            data["time_mp_cb_big_m"] = time_mp_cb_big_m
            data["time_sp_cb_big_m"] = time_sp_cb_big_m
            data["time_mis_cb_big_m"] = time_mis_cb_big_m

    if indicator_and_big_m:
        time_cb_indicator_and_big_m = df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]["time_cb_indicator_and_big_m"].to_list()
        time_cb_indicator_and_big_m.sort()
        time_cb_indicator_and_big_m.append(1800)
        instances_cb_indicator_and_big_m = len(time_cb_indicator_and_big_m) + 1
        instances_cb_indicator_and_big_m = np.arange(1, instances_cb_indicator_and_big_m)
        data["time_cb_indicator_and_big_m"] = time_cb_indicator_and_big_m
        data["instances_cb_indicator_and_big_m"] = instances_cb_indicator_and_big_m

        if extended:
            iter_cb_indicator_and_big_m = df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]["iter_cb_indicator_and_big_m"].to_list()
            iter_cb_indicator_and_big_m.sort()
            time_mp_cb_indicator_and_big_m= df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]["times_master_problem_cb_indicator_and_big_m"].to_list()
            time_mp_cb_indicator_and_big_m.sort()
            time_sp_cb_indicator_and_big_m = df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]["times_sub_problem_cb_indicator_and_big_m"].to_list()
            time_sp_cb_indicator_and_big_m.sort()
            time_mis_cb_indicator_and_big_m = df_data[df_data["termination_cb_indicator_and_big_m"] == "OPTIMAL"]["times_mis_problem_cb_indicator_and_big_m"].to_list()
            time_mis_cb_indicator_and_big_m.sort()
            data["iter_cb_indicator_and_big_m"] = iter_cb_indicator_and_big_m
            data["time_mp_cb_indicator_and_big_m"] = time_mp_cb_indicator_and_big_m
            data["time_sp_cb_indicator_and_big_m"] = time_sp_cb_indicator_and_big_m
            data["time_mis_cb_indicator_and_big_m"] = time_mis_cb_indicator_and_big_m

    return data 

def build_solved_instances_plot(colors, linestyles, markerstyles, ll_fba=False, ll_fba_indicator=False, big_m=False, file='csv/results_bigg_SCIP.csv', indicator=True, indicator_and_big_m=False, all_subplots=False, extended=False):
    data = load_data(big_m=big_m, indicator=indicator, indicator_and_big_m=indicator_and_big_m, file=file, ll_fba=ll_fba, ll_fba_indicator=ll_fba_indicator, extended=extended)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })
    
    # print(data)
    # create solved instances plot
    if all_subplots:
        fig, axs = plt.subplots(5, 1, figsize=[6.8, 9.6]) #, layout='constrained')
    else: 
        fig, axs = plt.subplots(2, 1, figsize=[6.8, 4.8]) #, layout='constrained')

    if ll_fba:
        axs[0].plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], label="ll FBA (big-M)", linestyle=linestyles[0])
    if ll_fba_indicator:
        axs[0].plot(data["time_ll_fba_indicator"], data["instances_ll_fba_indicator"], color=colors[5], label="ll FBA (indicator)", linestyle=linestyles[1])
    if indicator:
        axs[0].plot(data["time_cb"], data["instances_cb"], color=colors[2], linestyle=linestyles[2], label="CB (indicator)")
    if big_m:
        axs[0].plot(data["time_cb_big_m"], data["instances_cb_big_m"], color=colors[1], linestyle=linestyles[3], label="CB (big-M)")
    if indicator_and_big_m:    
        axs[0].plot(data["time_cb_indicator_and_big_m"], data["instances_cb_indicator_and_big_m"], color=colors[3], linestyle=linestyles[4], label="CB (indicator + big-M)")

    axs[0].set_ylabel(r"$\textbf{solved instances}$")
    axs[0].set_xlabel(r"$\textbf{time (s)}$")
    axs[0].grid(True)

    if indicator:
        axs[1].plot(data["iter_cb"], np.arange(1, len(data["iter_cb"])+1), color=colors[2], linestyle=linestyles[2])
    if big_m:
        axs[1].plot(data["iter_cb_big_m"], np.arange(1, len(data["iter_cb_big_m"])+1), color=colors[1], linestyle=linestyles[3])
    if indicator_and_big_m:    
        axs[1].plot(data["iter_cb_indicator_and_big_m"], np.arange(1, len(data["iter_cb_indicator_and_big_m"])+1), color=colors[3], linestyle=linestyles[4])
    
    axs[1].set_ylabel(r"$\textbf{solved instances}$")
    axs[1].set_xlabel(r"$\textbf{iterations}$")
    axs[1].grid(True)

    # if all_subplots:
        # axs[2].scatter(data["time_mp_cb"], np.arange(1, len(data["time_mp_cb"])+1), color=colors[1], marker=markerstyles[1], s=10, alpha=0.7)

        # axs[2].set_ylabel('solved instances')
        # axs[2].set_xlabel('average time (s) in master problem')
        # axs[2].grid(True)

        # axs[3].scatter(data["time_sp_cb"], np.arange(1, len(data["time_sp_cb"])+1), color=colors[0], marker=markerstyles[0], s=10, alpha=0.7)

        # axs[3].set_ylabel('solved instances')
        # axs[3].set_xlabel('average time (s) in sub problem')
        # axs[3].grid(True)

        # axs[4].scatter(data["time_mis_cb"], np.arange(1, len(data["time_mis_cb"])+1), color=colors[1], marker=markerstyles[1], s=10, alpha=0.7)

        # axs[4].set_ylabel('solved instances')
        # axs[4].set_xlabel('average time (s) in MIS search')
        # axs[4].grid(True)

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
    if indicator_and_big_m:
        plt.savefig("plots/comparison_solved_instances_indicator_and_big_m.pdf", format="pdf", bbox_inches="tight")
    if indicator:
        plt.savefig("plots/comparison_solved_instances.pdf", format="pdf", bbox_inches="tight")
    if big_m:
        plt.savefig("plots/comparison_solved_instances_big_m.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

# def build_time_vs_iterations_plot(colors, markerstyles, big_m=False, indicator=True, mis_list=[]):
#     data = load_data(big_m=big_m, indicator=indicator, mis_list=mis_list)

#     # create average time vs number of iterations plot
#     fig, axs = plt.subplots(2, 1, figsize=[6.8, 4.8])
#     axs[0].scatter(data["iter_cb"], data["time_mis_cb"], color=colors[0], marker=markerstyles[0], alpha=0.7)

#     for idx, mis in enumerate(mis_list):
#         axs[0].scatter(data["iter_cb_mis_" + str(mis)],data["time_mis_mis_" + str(mis)], color=colors[idx+2], marker=markerstyles[idx+2], alpha=0.7)

#     axs[0].set_xlabel('iterations')
#     axs[0].set_ylabel('average time (s) in MIS search')
#     axs[0].grid(True)

#     axs[1].scatter(data["iter_cb"], data["time_mp_cb"], label="cb", color=colors[1], marker=markerstyles[1], alpha=0.7)
    
#     for idx, mis in enumerate(mis_list):
#         axs[1].scatter(data["iter_cb_mis_" + str(mis)], data["time_mp_mis_" + str(mis)], label="mis_" + str(mis), color=colors[idx+2], marker=markerstyles[idx+2], alpha=0.7)

#     axs[1].set_xlabel('iterations')
#     axs[1].set_ylabel('average time (s) in master problem')
#     axs[1].grid(True)

#     lines = []
#     labels = []
#     for ax in fig.axes:
#         Line, Label = ax.get_legend_handles_labels()
#         # print(Label)
#         lines.extend(Line)
#         labels.extend(Label)

#     fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.05), ncol=3)
#     plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

#     # plt.show()
#     if indicator:
#         plt.savefig("plots/mis_comparison_time_vs_iterations.pdf", format="pdf", bbox_inches="tight")
#     if big_m:
#         plt.savefig("plots/mis_comparison_time_vs_iterations_big_m.pdf", format="pdf", bbox_inches="tight")
#     plt.clf()

# make plots
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#999999', '#984ea3', '#e41a1c', '#dede00', '#f781bf', '#a65628']
linestyles = ["dashdot", "dashed", "dotted", (0, (3, 1, 1, 1, 1, 1)), (0, (3, 5, 1, 5, 1, 5)), (0, (3, 10, 1, 10)), (0, (1, 10)), (0, (5, 10)), (0, (3, 1, 1, 1)), (0, (3, 5, 1, 5))]
markerstyles = ['o', 'v', '^', 's', 'p', 'D', 'd', 'p', 'D']

build_solved_instances_plot(colors, linestyles, markerstyles, big_m=True, indicator=True, indicator_and_big_m=True, all_subplots=False, extended=True, ll_fba=True, ll_fba_indicator=True)
