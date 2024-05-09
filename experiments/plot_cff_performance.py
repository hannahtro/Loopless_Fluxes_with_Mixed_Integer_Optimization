import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(filename='csv/results_cff.csv', indicator_and_big_m=False):
    # load df 
    df_data = pd.read_csv(filename)
    df_data = df_data.sort_values('time_loopless_fba')
    # print(df_data)
    data = {}
    data["organism"] = df_data["organism"]

    # ll fba data
    data["time_ll_fba"] = df_data["time_loopless_fba"].to_list()

    # cycles blocked 
    data["time_ll_fba_blocked_10"] = df_data["time_loopless_fba_blocked_10"].to_list()
    data["time_ll_fba_blocked_20"] = df_data["time_loopless_fba_blocked_20"].to_list()
    data["time_ll_fba_blocked_50"] = df_data["time_loopless_fba_blocked_50"].to_list()
    data["time_ll_fba_blocked_100"] = df_data["time_loopless_fba_blocked_100"].to_list()

    # shortest cycles blocked
    data["time_ll_fba_blocked_shortest_cycles_10"] = df_data["time_loopless_fba_blocked_shortest_cycles_10"].to_list()
    data["time_ll_fba_blocked_shortest_cycles_20"] = df_data["time_loopless_fba_blocked_shortest_cycles_20"].to_list()
    data["time_ll_fba_blocked_shortest_cycles_50"] = df_data["time_loopless_fba_blocked_shortest_cycles_50"].to_list()
    data["time_ll_fba_blocked_shortest_cycles_100"] = df_data["time_loopless_fba_blocked_shortest_cycles_100"].to_list()

    return data 


def build_solved_instances_plot(filename, colors, markerstyles, save_as="", indicator_and_big_m=False):
    data = load_data(filename, indicator_and_big_m=indicator_and_big_m)

    print(data.keys())
    # create solved instances plot
    fig, axs = plt.subplots(1, 1, figsize=[8.0, 3.8]) #, layout='constrained')

    axs.scatter(data["organism"], data["time_ll_fba_blocked_10"], color=colors[1], marker=markerstyles[1], label="cycles blocked (limit 10)", s=80)     
    axs.scatter(data["organism"], data["time_ll_fba_blocked_20"], color=colors[2], marker=markerstyles[2], label="cycles blocked (limit 20)", s=80)
    axs.scatter(data["organism"], data["time_ll_fba_blocked_50"], color=colors[3], marker=markerstyles[3], label="cycles blocked (limit 50)", s=80)
    axs.scatter(data["organism"], data["time_ll_fba_blocked_100"], color=colors[4], marker=markerstyles[4], label="cycles blocked (limit 100)", s=80)

    axs.scatter(data["organism"], data["time_ll_fba_blocked_shortest_cycles_10"], color=colors[5], marker=markerstyles[5], label="shortest cycles blocked (limit 10)", s=80)
    # axs.scatter(data["organism"], data["time_ll_fba_blocked_shortest_cycles_20"], color=colors[6], marker=markerstyles[6], label="shortest cycles blocked (limit 20)")
    # axs.scatter(data["organism"], data["time_ll_fba_blocked_shortest_cycles_50"], color=colors[7], marker=markerstyles[7], label="shortest cycles blocked (limit 50)")
    # axs.scatter(data["organism"], data["time_ll_fba_blocked_shortest_cycles_100"], color=colors[8], marker=markerstyles[8], label="shortest cycles blocked (limit 100)")

    axs.scatter(data["organism"], data["time_ll_fba"], color=colors[0], label="ll FBA", marker=markerstyles[0], s=80)

    axs.set_ylabel('time (s)')
    axs.grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.25), ncol=2)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.xticks(rotation=60)

    # plt.show()
    if save_as == "":
        plt.savefig("plots/comparison_solved_instances_gecko.pdf", format="pdf", bbox_inches="tight")
    else:
        plt.savefig(save_as, format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
markerstyles=["s", "D", "o", "^", "v", "h", "star"]

# filename = 'csv/results_bigg_SCIP_gecko.csv'
# build_solved_instances_plot(filename, colors, linestyles, save_as="plots/comparison_solved_instances_gecko.pdf")

build_solved_instances_plot("csv/results_cff.csv", colors, markerstyles, save_as="plots/cff_comparison.pdf")

# filename = 'csv/results_bigg_gecko_indicator_and_big_m_1.0e-8.csv'
# build_solved_instances_plot(filename, colors, linestyles, indicator_and_big_m=True, save_as="plots/comparison_solved_instances_gecko_1.0e-8_indicator_and_big_m.pdf")