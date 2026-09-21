'''
Produces Appendix Fig. A.1 of the paper: for each coronal region type
(BKC, AR, CO, FL-AVG), the composite SaXS spectral model (colored) together
with its individual isothermal apec components (grey plot), using the FULL,
unrestricted EMD and abundance fixed at 1.0 throughout (per the paper's caption).

This reads EMD data directly from the same CSV files already used in
coronal_regions_cleaned.py, via that file's own (unrestricted, default iter=4)
region objects: core, active_region, bkc, fl_all_fl_wt_avg_red_full

coronal_regions_cleaned.py's own model_creation_vvapec() bundles every
temperature bin into a single mdefine'd model, so PyXspec can't plot the bins
individually. This script instead builds each temperature bin as its own
separately-numbered apec component, so Plot.addComp() can plot them one by one.

NOTE: figures are saved automatically via savefig() into paper_plots/, but the
final plot layout in the paper is a manually-assembled 2x2 grid of all four
regions (bkc_and_components.pdf, ar_and_components.pdf, core_and_components.pdf,
fl_avg_and_components.pdf) — this script produces each panel individually.
'''

import matplotlib.pyplot as plt
from xspec import *
from coronal_regions import *
from matplotlib.ticker import FuncFormatter


def format_func(x, pos):
    if x >= 1:
        s = f'{x:.0f}'
    elif x >= 0.1:
        s = f'{x:.1f}'
    else:
        s = f'{x:.2f}'
    try:
        v = float(s)
        if abs(v - 0.8) < 1e-9 or abs(v - 0.9) < 1e-9:
            return ''
    except ValueError:
        pass
    return s


def read_full_emd(region, distance=4.84814e-6):
    '''
    Reads the EMD CSV file already associated with `region` (a
    coronal_class_spectra object from coronal_regions_cleaned.py) and converts
    temperature/normalization exactly as model_creation_vvapec() does — but
    always on the FULL, unrestricted EMD (no iter-based trimming), matching
    Appendix A.1's use of the unrestricted models.
    '''
    distance_cm = distance * 3.086e18
    temp, norm = [], []
    with open(region.temp_and_norm_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)
        for row in csvreader:
            temp.append(float(row[0]) * 1e-3 / 11606)
            norm.append(float(row[1]) * 1e44 * 1e10 * (1e-14 / (4 * np.pi * distance_cm ** 2)))
    return np.array(temp), np.array(norm)


def plot_model_with_components(region, label, color, spec_file, savefig_name, AllData=AllData):
    '''
    Plots `region`'s composite model (colored) together with each of its
    individual isothermal apec components (grey), matching Appendix Fig. A.1.
    Abundance is fixed at 1.0 throughout, per the paper.
    '''
    temp, norm = read_full_emd(region)
    n = len(temp)

    AllData.clear()
    AllModels.clear()
    Xset.abund = 'feld'

    mo_str = "+".join(["apec"] * n)
    m_test = Model(mo_str)
    comps = m_test.componentNames

    for i, comp in enumerate(comps):
        comp_obj = getattr(m_test, comp)
        comp_obj.kT.values = f"{temp[i]},-1"
        comp_obj.Abundanc.values = "1.0,-1"
        comp_obj.Redshift.values = "0.0,-1"
        comp_obj.norm.values = f"{norm[i]},-1"

    # A real spectrum + response is only used to establish the energy grid —
    # no data points are shown in the final plot, only the model.
    AllData += spec_file
    AllData.ignore("**-0.7 8.0-**")
    AllData.ignore("bad")

    Plot.device = '/null'
    Plot.add = True
    Plot.xAxis = "keV"
    Plot('data')
    energy = Plot.x()
    folded = Plot.model()

    fig, ax = plt.subplots()
    for i, comp in enumerate(comps):
        ax.plot(energy, Plot.addComp(i + 1), color="grey", linewidth=0.7,
                 label="APEC components" if i == 0 else None)
    ax.plot(energy, folded, color=color, linewidth=2, label=label)

    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_ylabel('counts/sec/keV', fontsize=16)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xmin=0.7, xmax=8.0)
    ax.legend(fontsize=12)
    ax.set_title(f"{label} model with APEC components", fontsize=14)

    ax.xaxis.set_major_formatter(FuncFormatter(format_func))
    ax.xaxis.set_minor_formatter(FuncFormatter(format_func))

    plt.tight_layout()
    plt.savefig(savefig_name)
    plt.show()
    plt.close("all")


spec_file = "data/spectra/minxss_fm3_PHA_1hr_avg127.pha"  # only used to set the energy grid/response

# Maps each region to its coronal_regions_cleaned.py object (unrestricted,
# default iter=4), display label, plot color, and paper_plots/ output filename.
REGIONS = {
    "bkc": (bkc, "BKC", "red", "paper_plots/bkc_and_components.pdf"),
    "ar": (active_region, "AR", "purple", "paper_plots/ar_and_components.pdf"),
    "core": (core, "CO", "orange", "paper_plots/core_and_components.pdf"),
    "fl_avg": (fl_all_fl_wt_avg_red_full, "FL-AVG", "blue", "paper_plots/fl_avg_and_components.pdf"),
}

# Change ONLY this line to switch which region's panel gets produced.
ACTIVE_REGION = "fl_avg"  # options: "bkc", "ar", "core", "fl_avg"

region_obj, label, color, savefig_name = REGIONS[ACTIVE_REGION]
plot_model_with_components(region_obj, label, color, spec_file, savefig_name)