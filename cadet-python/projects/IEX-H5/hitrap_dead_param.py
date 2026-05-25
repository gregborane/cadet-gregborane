# %%
import pandas as pd
import numpy as np
import sys, os
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_prominences

sys.path.append("..")

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

print(f"Dead Volume: {x[imax]}\n")
print(f"Dead Time: {(x[imax] / 5)}\n")

croissant, decroissant = [], []
peaks, _ = find_peaks(y)

prominences, left_bases, right_bases = peak_prominences(y, peaks)
y_peaks = [y[i] for i in peaks.tolist()]
x_peaks = [x[i] for i in peaks.tolist()]

x_right = [x[i] for i in right_bases.tolist()]
x_left = [x[i] for i in left_bases.tolist()]

y_right = [y[i] for i in right_bases.tolist()]
y_left = [y[i] for i in left_bases.tolist()]


def one_main_peak(xl: list, xr: list, yl: list, yr: list):
    global y_peaks, x_peaks

    xmax, ymax, xl_t, xr_t, yl_t, yr_t, i_l, i_r = 0, 0, 0, 0, 0, 0, 0, 0

    for i in range(len(xl)):
        if 100 < xr[i] < 122:
            xl_t = xr[i]
            yl_t = yr[i]
            i_l = i

        if 130 < xr[i] < 150:
            xr_t = xr[i]
            yr_t = yr[i]
            i_r = i + 1

    for j in range(len(y_peaks)):
        if 1500 < y_peaks[j] < 2000:
            xmax = x_peaks[j]
            ymax = y_peaks[j]

    return xl_t, xr_t, yl_t, yr_t, i_l, i_r, xmax, ymax


main_peak = one_main_peak(x_left, x_right, y_left, y_right)
W = main_peak[1] - main_peak[0]
sigma = W / 4
print(sigma)

tau = main_peak[-2] / Q
print(tau)

# %%
fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8), sharex=False)

# Original Data Plot
ax1.plot(x, y, "b.-", label="Original Data (y)")
ax1.set_xlabel("Volume in mL")
ax1.set_ylabel("Intensity in mAU")
ax1.plot(x_peaks, y_peaks, "x", markersize=16, label="peaks from scipy")
ax1.plot(main_peak[0], main_peak[2], "o", color="magenta", label="selected_points")
ax1.plot(main_peak[1], main_peak[3], "o", color="magenta")
ax1.plot(x, dy_dx, "r.-", label="Derivative (dy/dx)")
ax1.set_xlabel("Volume in mL")
ax1.set_ylabel("d(intensentiy)/d(volume)")
ax1.set_title("Original Data vs. Computed Derivative")
ax1.grid(True)
ax1.legend()
