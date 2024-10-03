import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_SCIP.csv', big_m=False, indicator=True, mis_list=[], remove_easy_instances=False, cut_densities=[], distinct_cuts=False, time_limit=1800, max_cuts=[]):
    assert big_m != indicator
    # load data 
    df_data = pd.read_csv(file)
    print("termination_ll_fba" in list(df_data.columns.values))
    print("termination_cb_mis_2.0_distinct_cuts" in list(df_data.columns.values))
    if remove_easy_instances:
        df_data = df_data[(df_data.time_ll_fba >= 15)]

    data = {}

    # ll fba data
    time_ll_fba = df_data[(df_data["termination_ll_fba"] == "OPTIMAL") & (df_data["time_ll_fba"] <= time_limit)]["time_ll_fba"].to_list()
    time_ll_fba.sort()
    time_ll_fba.append(time_limit)
    print(time_ll_fba)
    instances_ll_fba = len(time_ll_fba) + 1
    instances_ll_fba = np.arange(1, instances_ll_fba)
    data["time_ll_fba"] = time_ll_fba
    data["instances_ll_fba"] = instances_ll_fba

    # ll fba data indicator
    time_ll_fba_indicator = df_data[(df_data["termination_ll_fba_indicator"] == "OPTIMAL") & (df_data["time_ll_fba_indicator"] <= time_limit)]["time_ll_fba_indicator"].to_list()
    time_ll_fba_indicator.sort()
    time_ll_fba_indicator.append(time_limit)
    print(time_ll_fba_indicator)
    instances_ll_fba_indicator = len(time_ll_fba_indicator) + 1
    instances_ll_fba_indicator = np.arange(1, instances_ll_fba_indicator)
    data["time_ll_fba_indicator"] = time_ll_fba_indicator
    data["instances_ll_fba_indicator"] = instances_ll_fba_indicator

    # cb data
    if indicator:
        time_cb = df_data[(df_data["termination_cb"] == "OPTIMAL") & (df_data["time_cb"] <= time_limit)]["time_cb"].to_list()
    if big_m:
        time_cb = df_data[(df_data["termination_cb_big_m"] == "OPTIMAL") & (df_data["time_cb_big_m"] <= time_limit)]["time_cb_big_m"].to_list()
    time_cb.sort()
    time_cb.append(time_limit)
    instances_cb = len(time_cb) + 1
    instances_cb = np.arange(1, instances_cb)
    
    if indicator:
        iter_cb = df_data[(df_data["termination_cb"] == "OPTIMAL") & (df_data["time_cb"] <= time_limit)]["iter_cb"].to_list()
        iter_cb.sort()
        time_mp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_master_problem_cb"].to_list()
        time_mp_cb.sort()
        time_sp_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_sub_problem_cb"].to_list()
        time_sp_cb.sort()
        time_mis_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["times_mis_problem_cb"].to_list()

    if big_m:
        iter_cb = df_data[(df_data["termination_cb_big_m"] == "OPTIMAL") & (df_data["time_cb"] <= time_limit)]["iter_cb_big_m"].to_list()
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

    # cb mis data for indicator mis=2.0 and big m mis = 0.5 and distinct cuts and density
    if indicator:
        mis = 2.0
        time_cb_mis_temp = df_data[(df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis)] <= time_limit)]["time_cb_mis_" + str(mis)].to_list()
        time_cb_mis_temp.sort()
        time_cb_mis_temp.append(time_limit)
        instances_cb_mis_temp = len(time_cb_mis_temp) + 1
        instances_cb_mis_temp = np.arange(1, instances_cb_mis_temp)
        iter_cb_mis_temp = df_data[(df_data["termination_cb_mis_" + str(mis)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis)] <= time_limit)]["iter_cb_mis_" + str(mis)].to_list()
        iter_cb_mis_temp.sort()

        if distinct_cuts:
            time_cb_mis_distinct_cuts_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_distinct_cuts"] <= time_limit)]["time_cb_mis_" + str(mis) + "_distinct_cuts"].to_list()
            time_cb_mis_distinct_cuts_temp.sort()
            time_cb_mis_distinct_cuts_temp.append(time_limit)
            instances_cb_mis_distinct_cuts_temp = len(time_cb_mis_distinct_cuts_temp) + 1
            instances_cb_mis_distinct_cuts_temp = np.arange(1, instances_cb_mis_distinct_cuts_temp)
            iter_cb_mis_distinct_cuts_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_distinct_cuts"] <= time_limit)]["iter_cb_mis_" + str(mis) + "_distinct_cuts"].to_list()
            iter_cb_mis_distinct_cuts_temp.sort()

        for density in cut_densities:
            iter_cb_mis_density_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_density_" + str(density)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_density_" + str(density)] <= time_limit)]["iter_cb_mis_" + str(mis) + "_density_" + str(density)].to_list()
            iter_cb_mis_density_temp.sort()
            time_cb_mis_density_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_density_" + str(density)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_density_" + str(density)] <= time_limit) ]["time_cb_mis_" + str(mis) + "_density_" + str(density)].to_list()

            time_cb_mis_density_temp.sort()
            time_cb_mis_density_temp.append(time_limit)
            instances_cb_mis_density_temp = len(time_cb_mis_density_temp) + 1
            instances_cb_mis_density_temp = np.arange(1, instances_cb_mis_density_temp)

            data["time_cb_mis_" + str(mis) + "_density_" + str(density)] = time_cb_mis_density_temp
            data["instances_cb_mis_" + str(mis) + "_density_" + str(density)] = instances_cb_mis_density_temp
            data["iter_cb_mis_" + str(mis) + "_density_" + str(density)] = iter_cb_mis_density_temp

    if big_m:
        mis = 0.5
        time_cb_mis_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis)] <= time_limit)]["time_cb_big_m_mis_" + str(mis)].to_list()
        time_cb_mis_temp.sort()
        time_cb_mis_temp.append(time_limit)
        instances_cb_mis_temp = len(time_cb_mis_temp) + 1
        instances_cb_mis_temp = np.arange(1, instances_cb_mis_temp)
        iter_cb_mis_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis)] <= time_limit)]["iter_cb_big_m_mis_" + str(mis)].to_list()
        iter_cb_mis_temp.sort()

        if distinct_cuts:
            time_cb_mis_distinct_cuts_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis) + "_distinct_cuts"] <= time_limit)]["time_cb_big_m_mis_" + str(mis) + "_distinct_cuts"].to_list()
            time_cb_mis_distinct_cuts_temp.sort()
            time_cb_mis_distinct_cuts_temp.append(time_limit)
            instances_cb_mis_distinct_cuts_temp = len(time_cb_mis_distinct_cuts_temp) + 1
            instances_cb_mis_distinct_cuts_temp = np.arange(1, instances_cb_mis_distinct_cuts_temp)
            iter_cb_mis_distinct_cuts_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_distinct_cuts"] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_distinct_cuts"] <= time_limit)]["iter_cb_mis_" + str(mis) + "_distinct_cuts"].to_list()
            iter_cb_mis_distinct_cuts_temp.sort()

        for density in cut_densities:
            iter_cb_mis_density_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis) + "_density_" + str(density)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis) + "_density_" + str(density)] <= time_limit)]["iter_cb_big_m_mis_" + str(mis) + "_density_" + str(density)].to_list()
            iter_cb_mis_density_temp.sort()
            time_cb_mis_density_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis) + "_density_" + str(density)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis) + "_density_" + str(density)] <= time_limit) ]["time_cb_big_m_mis_" + str(mis) + "_density_" + str(density)].to_list()

            time_cb_mis_density_temp.sort()
            time_cb_mis_density_temp.append(time_limit)
            instances_cb_mis_density_temp = len(time_cb_mis_density_temp) + 1
            instances_cb_mis_density_temp = np.arange(1, instances_cb_mis_density_temp)

            data["time_cb_mis_" + str(mis) + "_density_" + str(density)] = time_cb_mis_density_temp
            data["instances_cb_mis_" + str(mis) + "_density_" + str(density)] = instances_cb_mis_density_temp
            data["iter_cb_mis_" + str(mis) + "_density_" + str(density)] = iter_cb_mis_density_temp

    data["time_cb_mis_" + str(mis)] = time_cb_mis_temp
    data["instances_cb_mis_" + str(mis)] = instances_cb_mis_temp
    data["iter_cb_mis_" + str(mis)] = iter_cb_mis_temp

    if distinct_cuts:
        data["time_cb_mis_" + str(mis) + "_distinct_cuts"] = time_cb_mis_distinct_cuts_temp
        data["iter_cb_mis_" + str(mis) + "_distinct_cuts"] = iter_cb_mis_distinct_cuts_temp
        data["instances_cb_mis_" + str(mis) + "_distinct_cuts"] = instances_cb_mis_distinct_cuts_temp

    # cb data for different mis and max cut values
    if indicator:
        for mis in mis_list:
            for m in max_cuts:
                try:
                    time_cb_mis_max_cut_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] <= time_limit)]["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)].to_list()
                    time_cb_mis_max_cut_temp.sort()
                    time_cb_mis_max_cut_temp.append(time_limit)
                    instances_cb_mis_max_cuts_temp = len(time_cb_mis_max_cut_temp) + 1
                    instances_cb_mis_max_cuts_temp = np.arange(1, instances_cb_mis_max_cuts_temp)
                    iter_cb_mis_max_cuts_temp = df_data[(df_data["termination_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] == "OPTIMAL") & (df_data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] <= time_limit)]["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)].to_list()
                    iter_cb_mis_max_cuts_temp.sort()

                    data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = time_cb_mis_max_cut_temp
                    data["instances_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = instances_cb_mis_max_cuts_temp
                    data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = iter_cb_mis_max_cuts_temp
                except:
                    print("mis ", mis, "max cuts ", m)
                    pass

    if big_m:
        for mis in mis_list:
            for m in max_cuts:
                try:
                    time_cb_mis_max_cut_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)] <= time_limit)]["time_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)].to_list()
                    time_cb_mis_max_cut_temp.sort()
                    time_cb_mis_max_cut_temp.append(time_limit)
                    instances_cb_mis_max_cuts_temp = len(time_cb_mis_max_cut_temp) + 1
                    instances_cb_mis_max_cuts_temp = np.arange(1, instances_cb_mis_max_cuts_temp)
                    iter_cb_mis_max_cuts_temp = df_data[(df_data["termination_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)] == "OPTIMAL") & (df_data["time_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)] <= time_limit)]["iter_cb_big_m_mis_" + str(mis) + "_max_cuts_" + str(m)].to_list()
                    iter_cb_mis_max_cuts_temp.sort()

                    data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = time_cb_mis_max_cut_temp
                    data["instances_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = instances_cb_mis_max_cuts_temp
                    data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)] = iter_cb_mis_max_cuts_temp
                except:                                        
                    print("mis ", mis, "max cuts ", m)
                    pass

    return data 


def build_solved_instances_plot(colors, linestyles, markerstyles, big_m=False, indicator=True, all_subplots=False, mis_list=[], distinct_cuts=False, cut_densities=[], time_limit=1800, max_cuts=[], remove_easy_instances=False):
    print("LOAD DATA")
    data = load_data(big_m=big_m, indicator=indicator, mis_list=mis_list, distinct_cuts=distinct_cuts, cut_densities=cut_densities, time_limit=time_limit, max_cuts=max_cuts, remove_easy_instances=remove_easy_instances)
    print(data.keys())

    print("BUILD PLOT")
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

    if big_m:
        mis = 0.5
        axs[0].plot(data["time_cb_mis_" + str(mis)], data["instances_cb_mis_" + str(mis)], color=colors[0], linestyle=linestyles[0], linewidth=2.0)

        if distinct_cuts:
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_distinct_cuts"], data["instances_cb_mis_" + str(mis) + "_distinct_cuts"], color=colors[1], linestyle=linestyles[1])

        for didx, density in enumerate(cut_densities):
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_density_" + str(density)], data["instances_cb_mis_" + str(mis) + "_density_" + str(density)], color=colors[didx+2], linestyle=linestyles[2])

        for idx, (mis, m) in enumerate([(0.5, 0.2), (0.5, 0.3), (0.5, 0.4), (1.0, 0.5), (2.0, 0.5), (5.0, 0.5)]):
            print(mis, m)         
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], data["instances_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], color=colors[idx+6], linestyle=linestyles[6])  

    if indicator: 
        mis = 2.0
        axs[0].plot(data["time_cb_mis_" + str(mis)], data["instances_cb_mis_" + str(mis)], color=colors[0], linestyle=linestyles[0], linewidth=2.0)
        
        if distinct_cuts:
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_distinct_cuts"], data["instances_cb_mis_" + str(mis) + "_distinct_cuts"], color=colors[1], linestyle=linestyles[1])

        for didx, density in enumerate(cut_densities):
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_density_" + str(density)], data["instances_cb_mis_" + str(mis) + "_density_" + str(density)], color=colors[didx+2], linestyle=linestyles[2])       

        for idx, (mis, m) in enumerate([(2.0, 0.5), (2.0, 1.0), (2.0, 1.5), (3.0, 2.0), (4.0, 2.0), (5.0, 2.0)]):
            print(mis, m)         
            axs[0].plot(data["time_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], data["instances_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], color=colors[idx+6], linestyle=linestyles[6])  


    axs[0].set_ylabel(r"$\textbf{solved instances}$")
    axs[0].set_xlabel(r"$\textbf{time (s)}$")
    axs[0].grid(True)

    if indicator: 
        axs[1].plot(data["iter_cb"], np.arange(1, len(data["iter_cb"])+1), label="CB (" + method + ")", color='green', linestyle='dotted', linewidth=2.0)
    if big_m:
        axs[1].plot(data["iter_cb"], np.arange(1, len(data["iter_cb"])+1), label="CB (" + method + ")", color='orange', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=2.0)

    
    if big_m:
        mis = 0.5
        axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[0], linestyle=linestyles[0], linewidth=2.0)

        if distinct_cuts:
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_distinct_cuts"], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_distinct_cuts"])+1), label="CB distinct cuts (" + method + " " + str(mis) + "\%)", color=colors[1], linestyle=linestyles[1])

        for didx, density in enumerate(cut_densities):
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_density_" + str(density)], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_density_" + str(density)])+1), label="CB max density " + str(density) + " (" + method + " " + str(mis) + "\%)", color=colors[didx+2], linestyle=linestyles[2])

        for idx, (mis, m) in enumerate([(0.5, 0.2), (0.5, 0.3), (0.5, 0.4), (1.0, 0.5), (2.0, 0.5), (5.0, 0.5)]):
            print(mis, m)
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)])+1), label="CB max cuts " + str(m) + " (" + method + " " + str(mis) + "\%)", color=colors[idx+6], linestyle=linestyles[6])
        
    if indicator:
        mis = 2.0
        axs[1].plot(data["iter_cb_mis_" + str(mis)], np.arange(1, len(data["iter_cb_mis_" + str(mis)])+1), label="CB (" + method + " " + str(mis) + "\%)", color=colors[0], linestyle=linestyles[0], linewidth=2.0)
        
        if distinct_cuts:
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_distinct_cuts"], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_distinct_cuts"])+1), label="CB distinct cuts (" + method + " " + str(mis) + "\%)", color=colors[1], linestyle=linestyles[1])

        for didx, density in enumerate(cut_densities):
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_density_" + str(density)], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_density_" + str(density)])+1), label="CB max density " + str(density) + " (" + method + " " + str(mis) + "\%)", color=colors[didx+2], linestyle=linestyles[2])

        for idx, (mis, m) in enumerate([(2.0, 0.5), (2.0, 1.0), (2.0, 1.5), (3.0, 2.0), (4.0, 2.0), (5.0, 2.0)]):
            print(mis, m)
            axs[1].plot(data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)], np.arange(1, len(data["iter_cb_mis_" + str(mis) + "_max_cuts_" + str(m)])+1), label="CB max cuts " + str(m) + " (" + method + " " + str(mis) + "\%)", color=colors[idx+6], linestyle=linestyles[6])
       
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

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.1), ncol=3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    # plt.show()
    if indicator:
        plt.savefig("plots/cut_selection_solved_instances.pdf", format="pdf", bbox_inches="tight")
    if big_m:
        plt.savefig("plots/cut_selection_solved_instances_big_m.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors = [
            'black', 
            'brown',
            'olive',
            'darkturquoise',
            'plum',
            'darkmagenta',
            '#377eb8', # blue
            '#ff7f00', # orange
            '#4daf4a', # green
            'indigo', # gray
            '#e41a1c', # red
            'goldenrod', # yellow
            '#f781bf', # pink
            'salmon',
            'turquoise',
            '#984ea3', # purple
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
            (0, (4, 4)), 
            (0, (2, 1, 1, 1)), 
            (0, (3, 1, 1, 1)),
            (0, (4, 4)), 
            (0, (2, 1, 1, 1)), 
            (0, (3, 1, 1, 1)),
            ]
markerstyles = ['o', 's', 'p', 'D', 'd', 'o', 'o']

### indicator mis plots 
build_solved_instances_plot(
    colors, 
    linestyles, 
    markerstyles, 
    big_m=False, 
    indicator=True, 
    all_subplots=False, 
    mis_list=[0.5,1.0,2.0,3.0,4.0,5.0], 
    distinct_cuts=True, 
    cut_densities=[5,10,15,20], 
    time_limit=300, 
    max_cuts=[0.5,1.0,1.5,2.0],
    remove_easy_instances=True
    )
# build_time_vs_iterations_plot(colors, markerstyles, big_m=False, indicator=True, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0, 20.0])

build_solved_instances_plot(
    colors, 
    linestyles, 
    markerstyles, 
    big_m=True, 
    indicator=False, 
    all_subplots=False, 
    mis_list=[0.5,1.0,2.0,5.0], 
    distinct_cuts=True, 
    cut_densities=[5,10,15,20], 
    time_limit=300, 
    max_cuts=[0.2,0.3,0.4,0.5],
    remove_easy_instances=True
    )
# build_time_vs_iterations_plot(colors, markerstyles, big_m=True, indicator=False, mis_list=[0.1, 0.5, 2.0, 5.0, 10.0])