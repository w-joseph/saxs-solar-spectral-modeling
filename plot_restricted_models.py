'''
Produces Appendix Fig. D.1: full (unrestricted) vs. EMD-restricted models,
for each region, saved automatically as paper_plots/restricted_models_<region>.pdf

Must be run from the repo root, since data files below are referenced
relative to it.
'''
from xspec import *
from coronal_regions import *
import matplotlib.pyplot as plt
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


def spec_fit(spec_file, mo_str, AllData=AllData):
    '''Fits mo_str (a "<region>+<region>_em_more_restricted"-style model pair) to spec_file.'''
    Plot.device = '/null'
    AllData.clear()
    AllModels.clear()
    Xset.abund = 'feld'

    m_test = Model(mo_str)
    comps = m_test.componentNames
    m_test.show()

    AllData += spec_file
    s1 = AllData(1)
    Fit.renorm('auto')

    # ignore channels below 0.7 keV and above 8.0 keV
    AllData.ignore("**-0.7 8.0-**")
    AllData.ignore("bad")

    Fit.renorm('auto')

    return Plot, comps


def plot_model(Plot, comps):
    Plot.device = '/null'
    Plot.add = True
    Plot.xAxis = "KeV"
    Plot('data')

    energy = Plot.x()
    folded = Plot.model()

    fig, axs = plt.subplots()

    for i, comp in enumerate(comps):
        linestyle = "-"
        if len(comps) > 1:
            if "bkc" in comp:
                c = "red"
                comp_label = "BKC"
            elif "core" in comp:
                c = "orange"
                comp_label = "CO"
            elif "ar" in comp:
                c = "purple"
                comp_label = "AR"
            elif "fl_x_9_0" in comp:
                c = "green"
                comp_label = "FL X9.0"
            elif "fl_all_fl_wt_avg" in comp:
                c = "blue"
                comp_label = "FL-AVG"

            if "em_more_restricted" in comp:
                comp_label = "EMD 1 order restricted"
                linestyle = ":"
            elif "em_restricted" in comp:
                comp_label = "EMD 2 order restricted"
                linestyle = "-."
            elif "e_restricted" in comp:
                comp_label = "energy restricted"
                linestyle = "--"
            elif comp == "fl_all_fl_wt_avg_red":
                comp_label = "EMD 1 order restricted"
                linestyle = ":"
            else:
                linestyle = "-"

            axs.plot(energy, Plot.addComp(i + 1), label=comp_label, color=c, linestyle=linestyle)

    plt.xlabel('Energy (keV)', fontsize=16)
    axs.set_ylabel('counts/sec/KeV', fontsize=16)
    axs.set_xscale("log")
    axs.set_yscale("log")
    axs.set_xlim(xmin=0.7, xmax=8.0)
    axs.legend(fontsize=14)

    plt.title("EMD restricted model vs. original", fontsize=12)
    plt.tight_layout()

    plt.setp(axs.get_xticklabels(), visible=True)
    axs.xaxis.set_major_formatter(FuncFormatter(format_func))
    axs.xaxis.set_minor_formatter(FuncFormatter(format_func))
    axs.tick_params(axis='x', which='major', labelsize=14)
    axs.tick_params(axis='x', which='minor', labelsize=14)
    axs.tick_params(axis='y', which='major', labelsize=14)

    plt.savefig('paper_plots/restricted_models_' + comps[0] + '.pdf')
    plt.show()
    plt.close("all")


spec_file = "data/spectra/minxss_fm3_PHA_1hr_avg127.pha"

# Each entry: full model + its EMD-restricted counterpart, for one region.
CONFIGS = {
    "core": "core+core_em_more_restricted",
    "ar": "ar+ar_em_more_restricted",
    "bkc": "bkc+bkc_em_more_restricted",
    "fl_x_9_0": "fl_x_9_0+fl_x_9_0_em_more_restricted",
    "fl_avg": "fl_all_fl_wt_avg_red_full+fl_all_fl_wt_avg_red",
}

# Change ONLY this line to switch which region's panel gets produced.
ACTIVE_CONFIG = "bkc"  # options: "core", "ar", "bkc", "fl_x_9_0", "fl_avg"

mo_str = CONFIGS[ACTIVE_CONFIG]
Plot, comps = spec_fit(spec_file, mo_str, AllData=AllData)
plot_model(Plot, comps)

# NOTE: an earlier version of this script also had an "fl_low_med_wt_avg_red_em_more_restricted"
# option, built from an older, superseded flare dataset (low+medium flares only,
# not including X9.0). It was replaced by "fl_avg" above, which uses the final
# all-flares dataset instead.