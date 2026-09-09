# %%
import pandas as pd
import sys, os
from dotenv import load_dotenv
import matplotlib.pyplot as plt

sys.path.append("..")
load_dotenv()

# %%
Q = 0.5  # ml/min
col_vol = 4.7  # mL

# %%
list_file_cobalt = [
    "cobalt_test_001_10-4_100microL.csv",
    "cobalt_test_002_10-4_100microL.csv",
]

data_cobalt = [f"{os.getenv('data')}/col_test/{i}" for i in list_file_cobalt]

hiscreen_data_cobalt = [
    pd.read_csv(data_file, encoding="utf-16", skiprows=2, sep="\t")
    for data_file in data_cobalt
]

data_cobalt_graph = [i.iloc[:, [0, 1]] for i in hiscreen_data_cobalt]
data_cobalt_graph = [i[i["ml"] >= 0] for i in data_cobalt_graph]

xcobalt = [i["ml"].to_list() for i in data_cobalt_graph]
ycobalt = [i["mAU"].to_list() for i in data_cobalt_graph]

idx_max_cobalt = [i.index(max(i)) for i in ycobalt]
dead_volumes_cobalt = [x_list[idx] for x_list, idx in zip(xcobalt, idx_max_cobalt)]
dead_times_cobalt = [i / Q for i in dead_volumes_cobalt]

print(f"Dead Volumes: {dead_volumes_cobalt}\n")
print(f"Dead Times: {dead_times_cobalt}\n")

dead_volume_cobalt = sum(dead_volumes_cobalt) / len(dead_volumes_cobalt)
dead_time_cobalt = sum(dead_times_cobalt) / len(dead_times_cobalt)

extra_porosity_cobalt = dead_volume_cobalt / col_vol

print(f"Dead Volume: {round(dead_volume_cobalt, 2)}\n")
print(f"Dead Time: {round(dead_time_cobalt, 2)}\n")
print(f"Extra porosity: {round(extra_porosity_cobalt, 2)}\n")

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8), sharex=False)
ax1.scatter(0, 0, label=f"[Co]2+ extra porosity : {round(extra_porosity_cobalt, 2)}")

colors = ["b", "g", "r"]
for i in range(len(xcobalt)):
    ax1.scatter(
        xcobalt[i],
        ycobalt[i],
        marker="x",
        color=colors[i],
        label=f"Dead volume with [Co]2+: {round(dead_volumes_cobalt[i], 2)} mL",
    )

# Total  %%
# --------------------------------------------------------
# --------------------------------------------------------
# --------------------------------------------------------
# --------------------------------------------------------

list_file_bteac = [
    "benzyl_10-4_001.csv",
    "benzyl_10-4_002.csv",
]

data_bteac = [f"{os.getenv('data')}/col_test/{i}" for i in list_file_bteac]

hiscreen_data_bteac = [
    pd.read_csv(data_file, encoding="utf-16", skiprows=2, sep="\t")
    for data_file in data_bteac
]

data_bteac_graph = [i.iloc[:, [0, 1]] for i in hiscreen_data_bteac]
data_bteac_graph = [i[i["ml"] >= 0] for i in data_bteac_graph]

xbteac = [i["ml"].to_list() for i in data_bteac_graph]
ybteac = [i["mAU"].to_list() for i in data_bteac_graph]

idx_max_bteac = [i.index(max(i)) for i in ybteac]
dead_volumes_bteac = [x_list[idx] for x_list, idx in zip(xbteac, idx_max_bteac)]
dead_times_bteac = [i / Q for i in dead_volumes_bteac]

print(f"Dead Volumes: {dead_volumes_bteac}\n")
print(f"Dead Times: {dead_times_bteac}\n")

dead_volume_bteac = sum(dead_volumes_bteac) / len(dead_volumes_bteac)
dead_time = sum(dead_times_bteac) / len(dead_times_bteac)

total_porosity = dead_volume_bteac / col_vol

print(f"Dead Volume: {round(dead_volume_bteac, 2)}\n")
print(f"Dead Time: {round(dead_time, 2)}\n")
print(f"Total porosity: {round(total_porosity, 2)}\n")

ax1.scatter(0, 0, label=f"BTEAC extra porosity : {round(total_porosity, 2)}")

colors = ["r", "orange", "p"]
for i in range(len(xbteac)):
    ax1.scatter(
        xbteac[i],
        ybteac[i],
        marker="x",
        color=colors[i],
        label=f"Dead volume with BTEAC : {round(dead_volumes_bteac[i], 2)} mL",
    )
    ax1.set_xlabel("Volume in mL")
    ax1.set_ylabel("Intensity in mAU")
    ax1.set_title("Porosity measurements with Cobalt and BTEAC")
    ax1.grid(True)
    ax1.legend()

lt.savefig("../graph/extra_chromatogram.png")
