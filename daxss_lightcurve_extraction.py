# daxss_lightcurve_extraction.py
# Loads the DAXSS level-2 1-hour-average mission file, saves a JD-vs-count-rate
# CSV (used downstream by goes_daxss_lc_plot.py), and plots the light curve
# with the active spectrum's index marked.

import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SPECTRUM_INDEX = {"flare": 127, "quiescent": 167}
ACTIVE = "flare"  # change this to mark a different spectrum on the plot

DAXSS_LEVEL2_FILE = "data/spectra/daxss_solarSXR_level2_1hour_average_2022-02-14-mission_v2.1.0.ncdf"

daxss = nc.Dataset(DAXSS_LEVEL2_FILE)
x_plot = daxss["TIME_JD"][:]
y_plot = daxss["X123_SLOW_CORRECTED"][:]

# Save JD vs. counts/sec CSV — consumed by goes_daxss_lc_plot.py
combined = np.asarray([x_plot, y_plot])
pd.DataFrame(np.transpose(combined), columns=["JD", "cts/sec"]).to_csv(
    "data/spectra/daxss_solarSXR_level2_1hour_average_2022-02-14-mission_v2.1.0.csv"
)

spectrum_index = SPECTRUM_INDEX[ACTIVE]

plt.rcParams["figure.figsize"] = [20, 10]
fig, ax = plt.subplots()
plt.xlabel("Time in JD")
plt.yscale("log")
plt.ylabel("Count Per Second (CPS)")
plt.suptitle("DAXSS Mission CPS vs Time")

plt.scatter(x_plot, y_plot, color="red", label="CPS")
plt.plot(x_plot[spectrum_index], y_plot[spectrum_index], "d", color="blue",
         label=f"Selected spectrum ({ACTIVE})")
plt.legend()
plt.show()