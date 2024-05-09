import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_dp_cleaned_up.csv', ll_fba=True, hull=True, big_m_1000=True, big_m_10000=True, dp_solvers=["HiGHS"], mis_numbers=[0.5], mis_big_m=True):
    # load data 
    df_data = pd.read_csv(file)
    # print(df_data)
    data = {}

    # ll fba data
    if ll_fba:
        time_ll_fba = df_data[df_data["termination_ll_fba"] == "OPTIMAL"]["time_ll_fba"].to_list()
        time_ll_fba.sort()
        time_ll_fba.append(1800)
        print(time_ll_fba)
        instances_ll_fba = len(time_ll_fba) + 1
        instances_ll_fba = np.arange(1, instances_ll_fba)
        data["time_ll_fba"] = time_ll_fba
        data["instances_ll_fba"] = instances_ll_fba

    for solver in dp_solvers:
        if hull:
            if solver == "HiGHS":
                time_dp_hull = df_data[df_data["termination_dp_Hull" ] == "OPTIMAL"]["time_dp_Hull"].to_list()
                time_dp_hull.sort()
                time_dp_hull.append(1800)
                instances_dp_hull = len(time_dp_hull) + 1
                instances_dp_hull = np.arange(1, instances_dp_hull)
                data["time_dp_hull"] = time_dp_hull
                data["instances_dp_hull"] = instances_dp_hull
            else:
                time_dp_hull_GLPK = df_data[df_data["termination_dp_Hull_GLPK"] == "OPTIMAL"]["time_dp_Hull_GLPK"].to_list()
                time_dp_hull_GLPK.sort()
                time_dp_hull_GLPK.append(1800)
                instances_dp_hull_GLPK = len(time_dp_hull_GLPK) + 1
                instances_dp_hull_GLPK = np.arange(1, instances_dp_hull_GLPK)
                data["time_dp_hull_GLPK"] = time_dp_hull_GLPK
                data["instances_dp_hull_GLPK"] = instances_dp_hull_GLPK
        
        if big_m_1000:
            if solver == "HiGHS":
                time_dp_big_m = df_data[df_data["termination_dp_BigM_1000"] == "OPTIMAL"]["time_dp_BigM_1000"].to_list()
                time_dp_big_m.sort()
                time_dp_big_m.append(1800)
                instances_dp_big_m = len(time_dp_big_m) + 1
                instances_dp_big_m = np.arange(1, instances_dp_big_m)
                data["time_dp_big_m"] = time_dp_big_m
                data["instances_dp_big_m"] = instances_dp_big_m
            else:
                time_dp_big_m_GLPK = df_data[df_data["termination_dp_BigM_1000_GLPK"] == "OPTIMAL"]["time_dp_BigM_1000_GLPK"].to_list()
                time_dp_big_m_GLPK.sort()
                time_dp_big_m_GLPK.append(1800)
                instances_dp_big_m_GLPK = len(time_dp_big_m_GLPK) + 1
                instances_dp_big_m_GLPK = np.arange(1, instances_dp_big_m_GLPK)
                data["time_dp_big_m_GLPK"] = time_dp_big_m_GLPK
                data["instances_dp_big_m_GLPK"] = instances_dp_big_m_GLPK

        if big_m_10000:
            if solver == "HiGHS":
                time_dp_big_m_10000 = df_data[df_data["termination_dp_BigM_10000"] == "OPTIMAL"]["time_dp_BigM_10000"].to_list()
                time_dp_big_m_10000.sort()
                time_dp_big_m_10000.append(1800)
                instances_dp_big_m_10000 = len(time_dp_big_m_10000) + 1
                instances_dp_big_m_10000 = np.arange(1, instances_dp_big_m_10000)
                data["time_dp_big_m_10000"] = time_dp_big_m_10000
                data["instances_dp_big_m_10000"] = instances_dp_big_m_10000
            else:
                time_dp_big_m_GLPK_10000 = df_data[df_data["termination_dp_BigM_10000_GLPK"] == "OPTIMAL"]["time_dp_BigM_10000"].to_list()
                time_dp_big_m_GLPK_10000.sort()
                time_dp_big_m_GLPK_10000.append(1800)
                instances_dp_big_m_GLPK_10000 = len(time_dp_big_m_GLPK_10000) + 1
                instances_dp_big_m_GLPK_10000 = np.arange(1, instances_dp_big_m_GLPK_10000)
                data["time_dp_big_m_GLPK_10000"] = time_dp_big_m_GLPK_10000
                data["instances_dp_big_m_GLPK_10000"] = instances_dp_big_m_GLPK_10000

        if mis_big_m:
            time_cb_big_m_mis = df_data[df_data["termination_cb_big_m_mis_0.5"] == "OPTIMAL"]["time_cb_big_m_mis_0.5"].to_list()
            time_cb_big_m_mis.sort()
            time_cb_big_m_mis.append(1800)
            print(time_cb_big_m_mis)
            instances_cb_big_m_mis = len(time_cb_big_m_mis) + 1
            instances_cb_big_m_mis = np.arange(1, instances_cb_big_m_mis)
            data["time_cb_big_m_mis"] = time_cb_big_m_mis
            data["instances_cb_big_m_mis"] = instances_cb_big_m_mis
    return data 

def build_solved_instances_plot(colors, linestyles, markerstyles, ll_fba=True, hull=True, big_m_1000=True, big_m_10000=False, mis_big_m=True, file='csv/results_bigg_dp_cleaned_up.csv', save_as='plots/comparison_dp_cleaned_up.pdf'):
    data = load_data(file=file, ll_fba=ll_fba, hull=hull, big_m_1000=big_m_1000, big_m_10000=big_m_10000)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })

    print(data)
    # create solved instances plot
    fig, ax = plt.subplots(figsize=[6.8, 2.4]) #, layout='constrained')

    if ll_fba:
        ax.plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], label="ll-FBA (big-M)", linestyle=linestyles[0])
    if mis_big_m:
        ax.plot(data["time_cb_big_m_mis"], data["instances_cb_big_m_mis"], color=colors[1], label="CB (big-M MIS 0.5 \%)", linestyle=(0, (3, 1, 1, 1, 1, 1)))
    if hull:
        ax.plot(data["time_dp_hull"], data["instances_dp_hull"], color=colors[4], linestyle=linestyles[4], label="DP (convex-hull)")
    if big_m_1000:
        ax.plot(data["time_dp_big_m"], data["instances_dp_big_m"], color=colors[2], linestyle='dashed', label="DP (big-M)")
    if big_m_10000:    
        ax.plot(data["time_dp_big_m_10000"], data["instances_dp_big_m_10000"], color=colors[3], linestyle='dotted', label="DP (big M 10000)")
    

    ax.set_ylabel(r"$\textbf{solved instances}$")
    ax.set_xlabel(r"$\textbf{time (s)}$")
    ax.grid(True)

    # if ll_fba:
    #     axs[1].plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], linestyle=linestyles[0])
    # if hull:
    #     axs[1].plot(data["time_dp_hull_GLPK"], data["instances_dp_hull_GLPK"], color=colors[1], linestyle=linestyles[1], label="DP (hull) GLPK")
    # if big_m_1000:
    #     axs[1].plot(data["time_dp_big_m_GLPK"], data["instances_dp_big_m_GLPK"], color=colors[2], linestyle=linestyles[2], label="DP (big M 1000) GLPK")
    # if big_m_10000:    
    #     axs[1].plot(data["time_dp_big_m_GLPK_10000"], data["instances_dp_big_m_GLPK_10000"], color=colors[3], linestyle=linestyles[3], label="DP (big M 10000) GLPK")
    
    # axs[1].set_ylabel('solved instances')
    # axs[1].set_xlabel('time (s)')
    # axs[1].grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.05), ncol=2)
    plt.tight_layout(pad=1.4, w_pad=0.5, h_pad=1.0)

    # plt.show()
    plt.savefig(save_as, format="pdf", bbox_inches="tight")
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
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#984ea3', '#999999', '#e41a1c', '#dede00', '#f781bf', '#a65628']
linestyles = ["-.", "--", "-.", ":", "-", "--", "-", "-", "--"]
markerstyles = ['o', 'v', '^', 's', 'p', 'D', 'd', 'p', 'D']

build_solved_instances_plot(colors, linestyles, markerstyles)
build_solved_instances_plot(colors, linestyles, markerstyles, save_as='plots/comparison_dp.pdf', file='csv/results_bigg_dp.csv')

# build_time_vs_iterations_plot(colors, markerstyles, big_m=False, indicator=True, mis_list=[5.0, 10.0, 20.0, 30.0])
# build_solved_instances_plot(colors, linestyles, markerstyles, big_m=True, indicator=False, all_subplots=False, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0])
# build_time_vs_iterations_plot(colors, markerstyles, big_m=True, indicator=False, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0])