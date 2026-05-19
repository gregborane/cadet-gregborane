# %%
import pandas as pd
import numpy as np
import sys, os
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

sys.path.append("..")
import subprocess

load_dotenv()

# %%
data = f"{os.getenv('data')}/IEX-HYDRO/Hitrap_Capto_Q_for_HA_29-01-2026_pH_8_2_001.xlsx"

hitrap_data = pd.read_excel(data)

# %%

Q = 5  # ml/min
hitrap_data.head(10)
hitrap_data.columns
x = hitrap_data["ml"].to_list()
y = hitrap_data["mAU"].to_list()
dy_dt = np.gradient(y)
dx_dt = np.gradient(x)

with np.errstate(divide="ignore", invalid="ignore"):
    dy_dx = dy_dt / dx_dt

val_max = 0
imax = 0
for i, val in enumerate(dy_dt):
    if val > val_max:
        val_max = val
        imax = i

print(f"max derivate {val_max} \n max indice {imax} \n dead vol {x[imax]}")

# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=False)

# Original Data Plot
ax1.plot(x, y, "b.-", label="Original Data (y)")
ax1.set_ylabel("y")
ax1.set_title("Original Data vs. Computed Derivative")
ax1.grid(True)
ax1.legend()

# Derivative Plot
ax2.plot(x, dy_dx, "r.-", label="Derivative (dy/dx)")
ax2.set_xlabel("x")
ax2.set_ylabel("dy/dx")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
