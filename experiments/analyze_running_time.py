import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

organisms = ["Ashbya_aceri", 
            "yHMPu5000035696_Hanseniaspora_singularis",
            "yHMPu5000035695_Hanseniaspora_pseudoguilliermondii",
            "yHMPu5000035684_Kloeckera_hatyaiensis",
            "yHMPu5000035659_Saturnispora_dispora",
            # "yHMPu5000034963_Hanseniaspora_clermontiae",
            "Tortispora_caseinolytica",
            "Starmerella_bombicola_JCM9596",
            "Hanseniaspora_uvarum",
            "Eremothecium_sinecaudum",
            "Eremothecium_gossypii"
]

# num_rows = 3
# num_cols= 4
plt.rcParams.update({
        "text.usetex": True,
        "font.family": "lmodern",
        "axes.titlesize" : "x-small"
})
fig, axs = plt.subplots(10, 1, figsize=[4.8, 12.8], layout='constrained')

for idx in range(0, len(organisms)):
    print(idx)
    file = 'json/' + organisms[idx] + '_combinatorial_benders_fast_big_m_14440.json'
    with open(file) as json_file:
        data = json.load(json_file)
        # print(data["cuts"])
        # print(data["times_sub_problem"])
        # print(len(data["times_sub_problem"]))

        # plt.plot(range(0,len(data["times_sub_problem"])), data['times_sub_problem'], label="subproblem")
        axs[idx].plot(range(0,len(data["times_master_problem"])), data['times_master_problem'], label="master problem")
        axs[idx].plot(range(0,len(data["times_mis_problem"])), data['times_mis_problem'], label="mis problem")
        axs[idx].xaxis.set_tick_params(labelsize=5)
        axs[idx].yaxis.set_tick_params(labelsize=5)
        axs[idx].title.set_text(organisms[idx])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
# plt.show()
plt.savefig("plots/analyze_running_time.pdf", format="pdf", bbox_inches="tight")