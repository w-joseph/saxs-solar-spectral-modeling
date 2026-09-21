'''
Produces Fig. 3 of the paper: GOES X-ray flux alongside DAXSS SXR count rate,
for either the flaring-Sun or quiescent-Sun observation.

NOTE: figure is saved manually, not via savefig(). Based on paper_plots/:
  "flare" -> goes_daxss_lc_flare.pdf
  "qui"   -> goes_daxss_lc_qui.pdf

Must be run from the repo root, since data files below are referenced
relative to it.
'''
import netCDF4 as nc
import numpy as np
import cftime
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter
from datetime import datetime, timedelta


def jd_to_time(jd, jd_ref=2451545.0):
    '''Formats a Julian Date as an "HH:MM" string, for the x-axis.'''
    dt = datetime(2000, 1, 1, 12) + timedelta(days=jd - jd_ref)
    return dt.strftime('%H:%M')


daxss_file = "data/spectra/daxss_solarSXR_level2_1hour_average_2022-02-14-mission_v2.1.0.csv"

# Each entry: (GOES file, Hinode observation time, date label for the plot
# title, and whether the marked DAXSS point is the count-rate max or min).
CONFIGS = {
    "flare": {
        "goes_file": "data/goes/sci_xrsf-l2-avg1m_g16_d20220425_v2-2-0.ncdf",
        "hinode_time": "2022-04-25T02:06:55",
        "date_string": "April 25, 2022",
        "marker_label": "DAXSS Peak",
        "use_max": True,
    },
    "qui": {
        "goes_file": "data/goes/sci_xrsf-l2-avg1m_g16_d20220629_v2-2-0.ncdf",
        "hinode_time": "2022-06-29T18:03:39",
        "date_string": "June 29, 2022",
        "marker_label": "DAXSS Min",
        "use_max": False,
    },
}

# Change ONLY this line to switch between the flare and quiescent panel.
ACTIVE_CONFIG = "qui"

config = CONFIGS[ACTIVE_CONFIG]
goes_file = config["goes_file"]

daxss = pd.read_csv(daxss_file)

if config["use_max"]:
    daxs_datapt_time = daxss['JD'][daxss['cts/sec'].idxmax()]
    daxs_datapt_cts = daxss['cts/sec'].max()
else:
    daxs_datapt_time = daxss['JD'][daxss['cts/sec'].idxmin()]
    daxs_datapt_cts = daxss['cts/sec'].min()

hinode_time_jd = (datetime.strptime(config["hinode_time"], "%Y-%m-%dT%H:%M:%S")
                   - datetime(2000, 1, 1, 12)).total_seconds() / (60 * 60 * 24) + 2451545.0

ff = nc.Dataset(goes_file)

time_jd = (ff.variables["time"][:] / (60 * 60 * 24)) + 2451545.0  # convert to JD

# GOES provides two energy bands (0.5-4A and 1-8A); only the 1-8A band is
# used here, since that's the band used for GOES flare classification
# (see paper Sect. 3.1) and is what's plotted in Fig. 3.
band_name = "1-8\u212B"
var_name = "xrsb_flux"
goes_classes = ["A", "B", "C", "M", "X"]
goes_values = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]

fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()

ax1.plot(
    time_jd,
    ff.variables[var_name][:],
    linewidth=1,
    color="green",
    label="GOES16 " + band_name,
)
ax2.scatter(
    daxss['JD'],
    daxss['cts/sec'],
    linewidth=1,
    color='red',
    label="DAXSS SXR",
)
ax2.scatter(
    daxs_datapt_time,
    daxs_datapt_cts,
    s=70,
    linewidth=3,
    color='blue',
    marker="x",
    label=config["marker_label"],
)
ax2.axvline(x=daxs_datapt_time - (0.5 / 24), color='gray', linestyle='--')
ax2.axvline(x=daxs_datapt_time + (0.5 / 24), color='gray', linestyle='--')

# vertical span marking the Hinode observation time +/- 5 minutes
ax2.axvspan(hinode_time_jd - (5 / 60 / 24), hinode_time_jd + (5 / 60 / 24), color='pink', alpha=0.5, label="Hinode +/- 5 min")

for y, label in zip(goes_values, goes_classes):
    ax1.axhline(y=y, color='gray', linewidth=1)
    ax1.text(x=time_jd[0] - (0.7 / 24), y=y, s=label, color='black', fontsize=14, ha='left', va='bottom', fontweight='bold')

ax1.set_yscale("log")
ax2.set_yscale("log")
ax2.legend(loc="upper right", prop={"size": 16})
ax1.legend(loc="upper center", prop={"size": 16})
ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: jd_to_time(x)))
ax1.set_xlim(time_jd[0] - (0.7 / 24), time_jd[-1])
ax1.tick_params(axis='x', which='major', labelsize=18)
ax2.tick_params(axis='y', which='major', labelsize=18)
ax1.tick_params(axis='y', which='major', labelsize=18)
ax1.set_xlabel("Time [hours]", fontsize=20)
ax1.set_ylabel("GOES X-ray Flux [{}]".format(ff[var_name].units), fontsize=20)
ax2.set_ylabel("DAXSS counts/sec", fontsize=20)
plt.title("GOES16 flux and DAXSS SXR cts/sec for " + config["date_string"], fontsize=20)
plt.show()