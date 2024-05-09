import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# load data 
df_data = pd.read_csv('csv/results_bigg_SCIP.csv')
print(df_data)

# ll fba data
time_ll_fba = df_data[df_data["termination_ll_fba"] == "OPTIMAL"]["time_ll_fba"].to_list()
time_ll_fba.sort()
time_ll_fba.append(1800)
print(time_ll_fba)
instances_ll_fba = len(time_ll_fba) + 1
instances_ll_fba = np.arange(1, instances_ll_fba)

# cb data
time_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["time_cb"].to_list()
time_cb.sort()
time_cb.append(1800)
instances_cb = len(time_cb) + 1
instances_cb = np.arange(1, instances_cb)
iter_cb = df_data[df_data["termination_cb"] == "OPTIMAL"]["iter_cb"].to_list()
iter_cb.sort()

fig, axs = plt.subplots(2, 1, layout='constrained')

axs[0].plot(time_ll_fba, instances_ll_fba)
axs[0].plot(time_cb, instances_cb)
# axs[0].set_xlim(0, 2)
axs[0].set_ylabel('solved instances')
axs[0].set_xlabel('time (s)')
axs[0].grid(True)

axs[1].plot(iter_cb, np.arange(1, len(iter_cb)+1))
# # axs[1].plot(iterations_cb, instances_cb)
# # axs[0].set_xlim(0, 2)
axs[1].set_ylabel('solved instances')
axs[1].set_xlabel('iterations')
axs[1].grid(True)

plt.show()