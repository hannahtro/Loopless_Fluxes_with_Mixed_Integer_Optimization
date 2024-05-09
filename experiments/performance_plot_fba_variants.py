import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_SCIP.csv', ll_fba=True, indicator=False, nullspace=True):
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

    if indicator:
        time_indicator = df_data[df_data["termination_ll_fba_indicator"] == "OPTIMAL"]["time_ll_fba_indicator"].to_list()
        time_indicator.sort()
        time_indicator.append(1800)
        print(time_indicator)
        instances_indicator = len(time_indicator) + 1
        instances_indicator = np.arange(1, instances_indicator)
        data["time_indicator"] = time_indicator
        data["instances_indicator"] = instances_indicator

    if nullspace:
        time_nullspace = df_data[df_data["termination_ll_fba_nullspace"] == "OPTIMAL"]["time_ll_fba_nullspace"].to_list()
        time_nullspace.sort()
        time_nullspace.append(1800)
        print(time_nullspace)
        instances_nullspace = len(time_nullspace) + 1
        instances_nullspace = np.arange(1, instances_nullspace)
        data["time_nullspace"] = time_nullspace
        data["instances_nullspace"] = instances_nullspace

    return data 


def build_solved_instances_plot(colors, linestyles, markerstyles, ll_fba=False, nullspace=True, indicator=False):
    data = load_data(ll_fba=ll_fba, indicator=indicator, nullspace=nullspace)

    plt.rcParams.update({
            "text.usetex": True,
            "font.family": "lmodern"
        })

    # create solved instances plot 
    fig, axs = plt.subplots(1, 1, figsize=[6.8, 2.4]) #, layout='constrained')

    if ll_fba:
        axs.plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], label="ll-FBA (big-M)", linestyle=linestyles[0])

    if indicator:
        axs.plot(data["time_indicator"], data["instances_indicator"], color=colors[5], label="ll-FBA (indicator)", linestyle=linestyles[1])

    if nullspace:
        axs.plot(data["time_nullspace"], data["instances_nullspace"], color=colors[2], label="ll-FBA (nullspace)", linestyle=linestyles[2], linewidth=2.0)

    axs.set_ylabel(r"$\textbf{solved instances}$")
    axs.set_xlabel(r"$\textbf{time (s)}$")
    axs.grid(True)

    lines = []
    labels = []
    for ax in fig.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)

    fig.legend(lines, labels, loc='center', bbox_to_anchor=(0.52, -0.05), ncol=3)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    plt.savefig("plots/fba_variants_comparison_plot.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#999999', '#984ea3', '#e41a1c', '#dede00', '#f781bf', '#a65628']
linestyles = ["dashdot", "dashed", "dotted", (0, (1, 10)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 5, 1, 5)), (0, (3, 1, 1, 1)), (0, (3, 5, 1, 5, 1, 5))]
markerstyles = ['o', 'v', '^', 's', 'p', 'D', 'd', 'p', 'D']

build_solved_instances_plot(colors, linestyles, markerstyles, nullspace=True, indicator=True, ll_fba=True)