import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_SCIP_cut_selection.csv', big_m=False, indicator=True, mis_list=[], remove_easy_instances=False, cut_densities=[], distinct_cuts=False):
    assert big_m != indicator
    # load data 
    df_data = pd.read_csv(file)
    # print(df_data)
    if remove_easy_instances:
        df_data[(df_data.time_ll_fba <= 15)]

    data = {}

    # ll fba data
    time_ll_fba = df_data[df_data["termination_ll_fba"] == "OPTIMAL"]["time_ll_fba"].to_list()
    time_ll_fba.sort()
    time_ll_fba.append(1800)
    print(time_ll_fba)
    instances_ll_fba = len(time_ll_fba) + 1
    instances_ll_fba = np.arange(1, instances_ll_fba)
    data["time_ll_fba"] = time_ll_fba
    data["instances_ll_fba"] = instances_ll_fba

    # ll fba data indicator
    time_ll_fba_indicator = df_data[df_data["termination_ll_fba_indicator"] == "OPTIMAL"]["time_ll_fba_indicator"].to_list()
    time_ll_fba_indicator.sort()
    time_ll_fba_indicator.append(1800)
    print(time_ll_fba_indicator)
    instances_ll_fba_indicator = len(time_ll_fba_indicator) + 1
    instances_ll_fba_indicator = np.arange(1, instances_ll_fba_indicator)
    data["time_ll_fba_indicator"] = time_ll_fba_indicator
    data["instances_ll_fba_indicator"] = instances_ll_fba_indicator

    # cb data
    if indicator:
        time_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["time_cb"].to_list()
    if big_m:
        time_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["time_cb_big_m"].to_list()
    time_cb.sort()
    time_cb.append(1800)
    instances_cb = len(time_cb) + 1
    instances_cb = np.arange(1, instances_cb)
    
    if indicator:
        iter_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["iter_cb"].to_list()
        iter_cb.sort()
        time_mp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_master_problem_cb"].to_list()
        time_mp_cb.sort()
        time_sp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_sub_problem_cb"].to_list()
        time_sp_cb.sort()
        time_mis_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_mis_problem_cb"].to_list()

    if big_m:
        iter_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["iter_cb_big_m"].to_list()
        iter_cb.sort()
        time_mp_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_master_problem_cb_big_m"].to_list()
        time_mp_cb.sort()
        time_sp_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_sub_problem_cb_big_m"].to_list()
        time_sp_cb.sort()
        time_mis_cb = df_data[df_data["termination_cb_big_m"] == "OPTIMAL"]["times_mis_problem_cb_big_m"].to_list()
    
    time_mis_cb.sort()
    data["time_cb"] = time_cb
    data["instances_cb"] = instances_cb
    data["iter_cb"] = iter_cb
    data["time_mp_cb"] = time_mp_cb
    data["time_sp_cb"] = time_sp_cb
    data["time_mis_cb"] = time_mis_cb

    # cb mis data
    for mis in mis_list:
        if indicator:
            time_cb_mis_temp = df_data[df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL"]["time_cb_mis_" + str(mis)].to_list()
            if distinct_cuts:
                time_cb_mis_distinct_cuts_temp = df_data[df_data["termination_cb_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL"]["time_cb_mis_" + str(mis) + "_distinct_cuts"].to_list()
        
        if big_m:
            time_cb_mis_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL"]["time_cb_big_m_mis_" + str(mis)].to_list()
            if distinct_cuts:
                time_cb_mis_distinct_cuts_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL"]["time_cb_big_m_mis_" + str(mis) + "_distinct_cuts"].to_list()
        
        time_cb_mis_temp.sort()
        time_cb_mis_temp.append(1800)
        instances_cb_mis_temp = len(time_cb_mis_temp) + 1
        instances_cb_mis_temp = np.arange(1, instances_cb_mis_temp)

        if distinct_cuts:
            time_cb_mis_distinct_cuts_temp.sort()
            time_cb_mis_distinct_cuts_temp.append(1800)
            instances_cb_mis_distinct_cuts_temp = len(time_cb_mis_temp) + 1
            instances_cb_mis_distinct_cuts_temp = np.arange(1, instances_cb_mis_distinct_cuts_temp)

        if indicator:
            iter_cb_mis_temp = df_data[df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL"]["iter_cb_mis_" + str(mis)].to_list()
            iter_cb_mis_temp.sort()
            time_mp_mis_temp = df_data[df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL"]["times_master_problem_mis_" + str(mis)].to_list()
            time_mp_mis_temp.sort()
            time_sp_mis_temp = df_data[df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL"]["times_sub_problem_mis_" + str(mis)].to_list()
            time_sp_mis_temp.sort()
            time_mis_mis_temp = df_data[df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL"]["times_mis_problem_mis_" + str(mis)].to_list()
            if distinct_cuts:
                iter_cb_mis_distinct_cuts_temp = df_data[df_data["termination_cb_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL"]["iter_cb_mis_" + str(mis) + "_distinct_cuts"].to_list()
                iter_cb_mis_distinct_cuts_temp.sort()
        if big_m: 
            iter_cb_mis_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL"]["iter_cb_big_m_mis_" + str(mis)].to_list()
            iter_cb_mis_temp.sort()
            time_mp_mis_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL"]["times_master_problem_big_m_mis_" + str(mis)].to_list()
            time_mp_mis_temp.sort()
            time_sp_mis_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL"]["times_sub_problem_big_m_mis_" + str(mis)].to_list()
            time_sp_mis_temp.sort()
            time_mis_mis_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL"]["times_mis_problem_big_m_mis_" + str(mis)].to_list()
            if distinct_cuts:
                iter_cb_mis_distinct_cuts_temp = df_data[df_data["termination_cb_big_m_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL"]["iter_cb_big_mmis_" + str(mis) + "_distinct_cuts"].to_list()
                iter_cb_mis_distinct_cuts_temp.sort()
        
        time_mis_mis_temp.sort()
        data["time_cb_mis_" + str(mis)] = time_cb_mis_temp
        data["instances_cb_mis_" + str(mis)] = instances_cb_mis_temp
        data["iter_cb_mis_" + str(mis)] = iter_cb_mis_temp
        data["time_mp_mis_" + str(mis)] = time_mp_mis_temp
        data["time_sp_mis_" + str(mis)] = time_sp_mis_temp
        data["time_mis_mis_" + str(mis)] = time_mis_mis_temp
        data["time_cb_mis_" + str(mis) + "_distinct_cuts"] = time_cb_mis_distinct_cuts_temp
        data["iter_cb_mis_" + str(mis) + "_distinct_cuts"] = iter_cb_mis_distinct_cuts_temp

    return data 


def build_solved_instances_plot(colors, linestyles, markerstyles, big_m=False, indicator=True, all_subplots=False, mis_list=[], distinct_cuts=False, cut_densities=[]):
    data = load_data(big_m=big_m, indicator=indicator, mis_list=mis_list, distinct_cuts=distinct_cuts, cut_densities=cut_densities)
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })

    if big_m: 
        method = "big-M"
    elif indicator:
        method = "indicator"


    # create solved instances plot
    if all_subplots:
        fig, axs = plt.subplots(5, 1, figsize=[6.8, 9.6]) #, layout='constrained')
    else: 
        fig, axs = plt.subplots(2, 1, figsize=[6.8, 4.8]) #, layout='constrained')

    if indicator:
        axs[0].plot(data["time_ll_fba_indicator"], data["instances_ll_fba_indicator"], color='red', label="ll-FBA (indicator)", linestyle='dashed', linewidth=2.0)
        axs[0].plot(data["time_cb"], data["instances_cb"], color='green', linestyle='dotted', linewidth=2.0)
    else:
        axs[0].plot(data["time_ll_fba"], data["instances_ll_fba"], color='#377eb8', label="ll-FBA (big-M)", linestyle='dashdot', linewidth=2.0)
        axs[0].plot(data["time_cb"], data["instances_cb"], color='orange', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=2.0)

    for idx, mis in enumerate(mis_list):
        if big_m:
            axs[0].plot(data["time_cb_mis_" + str(mis)], data["instances_cb_mis_" + str(mis)], color=colors[idx], linestyle=linestyles[idx])
        if indicator:
            axs[0].plot(data["time_cb_mis_" + str(mis)], data["instances_cb_mis_" + str(mis)], color=colors[idx], linestyle=linestyles[idx])

    axs[0].set_ylabel(r"$\textbf{solved instances}$")
    axs[0].set_xlabel(r"$\textbf{time (s)}$")
    axs[0].grid(True)

    if indicator: 
        axs[1].plot(data["iter_cb"], np.arange(1, len(data["iter_cb"])+1), label="CB (" + method + ")", color='green', linestyle='dotted', linewidth=2.0)
    if big_m:
        axs[1].plot(data["iter_cb"], np.arange(1, len(data["iter_cb"])+1), label="CB (" + method + ")", color='orange', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=2.0)

    
    for idx, mis in enumerate(mis_list):
        if big_m:
            if mis == 0.5:
                axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[idx], linestyle=linestyles[idx], linewidth=2.0)
            else:
                axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[idx], linestyle=linestyles[idx])
        if indicator:
            if mis == 2.0:
                axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[idx], linestyle=linestyles[idx], linewidth=2.0)
            else:
                axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[idx], linestyle=linestyles[idx])

    axs[1].set_ylabel(r"$\textbf{solved instances}$")
    axs[1].set_xlabel(r"$\textbf{iterations}$")
    axs[1].grid(True)

    # if all_subplots:
    #     axs[2].scatter(data["time_mp_cb"], np.arange(1, len(data["time_mp_cb"])+1), color=colors[5], marker=markerstyles[1], s=10, alpha=0.6)

    #     for idx, mis in enumerate(mis_list):
    #         axs[2].scatter(data["time_mp_mis_" + str(mis)], np.arange(1, len(data["time_mp_mis_" + str(mis)])+1), color=colors[idx+2], marker=markerstyles[idx+2], s=10, alpha=0.6)

    #     axs[2].set_ylabel('solved instances')
    #     axs[2].set_xlabel('average time (s) in master problem')
    #     axs[2].grid(True)

    #     axs[3].scatter(data["time_sp_cb"], np.arange(1, len(data["time_sp_cb"])+1), color=colors[0], marker=markerstyles[0], s=10, alpha=0.6)

    #     for idx, mis in enumerate(mis_list):
    #         axs[3].scatter(data["time_sp_mis_" + str(mis)], np.arange(1, len(data["time_sp_mis_" + str(mis)])+1), color=colors[idx+2], marker=markerstyles[idx+2], s=10, alpha=0.6)

    #     axs[3].set_ylabel('solved instances')
    #     axs[3].set_xlabel('average time (s) in sub problem')
    #     axs[3].grid(True)

    #     axs[4].scatter(data["time_mis_cb"], np.arange(1, len(data["time_mis_cb"])+1), color=colors[1], marker=markerstyles[1], s=10, alpha=0.6)
    #     for idx, mis in enumerate(mis_list):
    #         axs[4].scatter(data["time_mis_mis_" + str(mis)], np.arange(1, len(data["time_mis_mis_" + str(mis)])+1), color=colors[idx+2], marker=markerstyles[idx+2], s=10, alpha=0.6)

    #     axs[4].set_ylabel('solved instances')
    #     axs[4].set_xlabel('average time (s) in MIS search')
    #     axs[4].grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.1), ncol=3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    # plt.show()
    if indicator:
        plt.savefig("plots/mis_comparison_solved_instances.pdf", format="pdf", bbox_inches="tight")
    if big_m:
        plt.savefig("plots/mis_comparison_solved_instances_big_m.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

def build_time_vs_iterations_plot(colors, markerstyles, big_m=False, indicator=True, mis_list=[]):
    data = load_data(big_m=big_m, indicator=indicator, mis_list=mis_list, remove_easy_instances=True)

    print(data["time_ll_fba"])

    # create average time vs number of iterations plot
    fig, axs = plt.subplots(2, 1, figsize=[6.8, 4.8])
    if indicator:
        axs[0].scatter(data["iter_cb"], data["time_mis_cb"], color='green', marker='v', alpha=0.6)
    if big_m:
        axs[0].scatter(data["iter_cb"], data["time_mis_cb"], color='orange', marker='^', alpha=0.6)

    for idx, mis in enumerate(mis_list):
        axs[0].scatter(data["iter_cb_mis_" + str(mis)],data["time_mis_mis_" + str(mis)], color=colors[idx], marker=markerstyles[idx], alpha=0.6)

    axs[0].set_xlabel(r"$\textbf{total number of iterations}$")
    axs[0].set_ylabel(r"$\textbf{fraction of average time (s)}$"
                    "\n"
                    r"$\textbf{in MIS search}$")
    axs[0].grid(True)

    if indicator:
        axs[1].scatter(data["iter_cb"], [a/b for a,b in zip(data["time_mp_cb"],data["time_cb"])], label="CB (indicator)", color='green', marker='v', alpha=0.6)
    else:    
        axs[1].scatter(data["iter_cb"], [a/b for a,b in zip(data["time_mp_cb"],data["time_cb"])], label="CB (big-M)", color='orange', marker='^', alpha=0.6)
    
    for idx, mis in enumerate(mis_list):
        if indicator:
            axs[1].scatter(data["iter_cb_mis_" + str(mis)], [a/b for a,b in zip(data["time_mp_mis_" + str(mis)], data["time_cb_mis_" + str(mis)])], label="CB (indicator MIS " + str(mis) + "\%)", color=colors[idx], marker=markerstyles[idx], alpha=0.6)
        if big_m:
            axs[1].scatter(data["iter_cb_mis_" + str(mis)], [a/b for a,b in zip(data["time_mp_mis_" + str(mis)],data["time_cb_mis_" + str(mis)])], label="CB (big-M MIS " + str(mis) + "\%)", color=colors[idx], marker=markerstyles[idx], alpha=0.6)

    axs[1].set_xlabel(r"$\textbf{total number of iterations}$")
    axs[1].set_ylabel(r"$\textbf{fraction of average time (s)}$"
                    "\n"
                    r"$\textbf{in MP}$")
    axs[1].grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.5, -0.1), ncol=3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    # plt.show()
    if indicator:
        plt.savefig("plots/mis_comparison_time_vs_iterations.pdf", format="pdf", bbox_inches="tight")
    if big_m:
        plt.savefig("plots/mis_comparison_time_vs_iterations_big_m.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors = [
        # '#377eb8', # blue
        # '#ff7f00', # orange
        # '#4daf4a', # green
        #   'indigo', # gray
        # '#e41a1c', # red
        #   'goldenrod', # yellow
        #   '#f781bf', # pink
          'black', 
        #   'salmon',
        #   'turquoise',
          'brown',
          'olive',
          'plum',
        #   '#984ea3', # purple
        #   'turquoise',
          'darkmagenta',
          'darkturquoise',
          ]
linestyles = [
            # "dashdot", 
            # "dashed", 
            # "dotted", 
            # (0, (3, 1, 1, 1, 1, 1)), 
            # (0, (1, 1)), 
            (0, (3, 3, 1, 3)),
            # (0, (3, 2, 1, 2, 1, 2)), 
            (0, (1, 3, 1, 2, 1, 3)),
            (0, (3, 2, 1, 2)), 
            # (0, (1, 2)), 
            (0, (4, 4)), 
            (0, (2, 1, 1, 1)), 
            (0, (3, 1, 1, 1)), 
            ]
markerstyles = ['o', 's', 'p', 'D', 'd', 'o', 'o']

### indicator mis plots 
build_solved_instances_plot(colors, linestyles, markerstyles, big_m=False, indicator=True, all_subplots=False, mis_list=[2.0], distinct_cuts=True, cut_densities=[5,10,15])
# build_time_vs_iterations_plot(colors, markerstyles, big_m=False, indicator=True, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0, 20.0])

build_solved_instances_plot(colors, linestyles, markerstyles, big_m=True, indicator=False, all_subplots=False, mis_list=[0.5], distinct_cuts=True, cut_densities=[5,10,15])
# build_time_vs_iterations_plot(colors, markerstyles, big_m=True, indicator=False, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0])