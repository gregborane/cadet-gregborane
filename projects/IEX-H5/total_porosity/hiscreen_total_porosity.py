# %%
import pandas as pd
import numpy as np
import sys, os
from dotenv import load_dotenv
import matplotlib.pyplot as plt

sys.path.append("..")
load_dotenv()

# %%
list_file_aceton = [
    "Acetone_test_002-100microL.csv",
    "Acetone_test_003_100microL.csv",
    "Acetone_test_004_100microL.csv",
]

data_aceton = [f"{os.getenv('data')}/col_test/{i}" for i in list_file_aceton]

hiscreen_data_aceton = [
    pd.read_csv(data_file, encoding="utf-16", skiprows=2, sep="\t")
    for data_file in data_aceton
]

# %%
Q = 0.5  # ml/min
col_vol = 4.7  # mL

# %%
data_aceton_graph = [i.iloc[:, [0, 1]] for i in hiscreen_data_aceton]
data_aceton_graph = [i[i["ml"] >= 0] for i in data_aceton_graph]

# %%
x = [i["ml"].to_list() for i in data_aceton_graph]
y = [i["mAU"].to_list() for i in data_aceton_graph]

idx_max = [i.index(max(i)) for i in y]
dead_volumes = [x_list[idx] for x_list, idx in zip(x, idx_max)]
dead_times = [i / Q for i in dead_volumes]

print(f"Dead Volumes: {dead_volumes}\n")
print(f"Dead Times: {dead_times}\n")

dead_volume = sum(dead_volumes) / len(dead_volumes)
dead_time = sum(dead_times) / len(dead_times)

total_porosity = dead_volume / col_vol

print(f"Dead Volume: {round(dead_volume, 2)}\n")
print(f"Dead Time: {round(dead_time, 2)}\n")
print(f"Total porosity: {round(total_porosity, 2)}\n")

# %%
fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8), sharex=False)
ax1.scatter(0, 0, label=f"total porosity : {round(total_porosity, 2)}")

colors = ["b", "g", "r"]
for i in range(len(x)):
    ax1.scatter(
        x[i],
        y[i],
        marker="x",
        color=colors[i],
        label=f"dead volume : {round(dead_volumes[i], 2)} mL",
    )
    ax1.set_xlabel("Volume in mL")
    ax1.set_ylabel("Intensity in mAU")
    ax1.set_title("Porosity measurements with acetone")
    ax1.grid(True)
    ax1.legend()
plt.savefig("../graph/acetone_chromatogram.png")
