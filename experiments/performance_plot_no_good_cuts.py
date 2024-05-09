import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_data(file='csv/results_bigg_SCIP_no_good_cuts_big_m.csv', ll_fba=True, no_good_cut=True):
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


    if no_good_cut:
        time_no_good_cut = df_data[df_data["termination_no_good_cuts_big_m"] == "OPTIMAL"]["time_no_good_cuts_big_m"].to_list()
        time_no_good_cut.sort()
        time_no_good_cut.append(1800)
        print(time_no_good_cut)
        instances_no_good_cut = len(time_no_good_cut) + 1
        instances_no_good_cut = np.arange(1, instances_no_good_cut)
        data["time_no_good_cut"] = time_no_good_cut
        data["instances_no_good_cut"] = instances_no_good_cut

    return data 


def build_solved_instances_plot(colors, linestyles, markerstyles, ll_fba=False, no_good_cut=False):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern"
    })

    data = load_data(ll_fba=ll_fba, no_good_cut=no_good_cut)

    # create solved instances plot 
    fig, axs = plt.subplots(1, 1, figsize=[6.8, 2.4]) #, layout='constrained')

    if ll_fba:
        axs.plot(data["time_ll_fba"], data["instances_ll_fba"], color=colors[0], label="ll-FBA (big-M)", linestyle=linestyles[0])

    if no_good_cut:
        axs.plot(data["time_no_good_cut"], data["instances_no_good_cut"], color=colors[4], label="no-good cut", linestyle=linestyles[5])

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

    plt.savefig("plots/no_good_cuts_comparison_plot.pdf", format="pdf", bbox_inches="tight")
    plt.clf()

# make plots
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#984ea3', '#999999', '#e41a1c', '#dede00', '#f781bf', '#a65628']
linestyles = ["-.", "--", "-.", ":", "-", "--", "-", "-", "--"]
markerstyles = ['o', 'v', '^', 's', 'p', 'D', 'd', 'p', 'D']

build_solved_instances_plot(colors, linestyles, markerstyles, no_good_cut=True, ll_fba=True)