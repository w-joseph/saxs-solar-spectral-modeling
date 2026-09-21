'''
Produces Fig. 4 of the paper: the quiescent- and flaring-Sun DAXSS spectra
with arbitrarily-scaled (not fit) SaXS coronal-region models overlaid, using
the FULL, unrestricted EMDs. This demonstrates the steeper observed
high-energy drop-off compared to the unrestricted models, motivating the
EMD restriction used everywhere else in the pipeline (see coronal_regions_cleaned.py).

NOTE: figure is saved manually, not via savefig(). Saved as
paper_plots/bkc_ar_core_fl_both_spectra_bad_demo.pdf.
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


def plot_spectra_with_unrestricted_models(spec_file_1, spec_file_2, mo_str, AllData=AllData):
    '''
    Plots two DAXSS spectra (flaring and quiescent Sun) together, with
    arbitrarily-scaled (unfit) coronal-region models overlaid, using full,
    unrestricted EMDs. Not a fit — the component normalizations below are
    manually chosen only to make each region's model visible on the plot.
    '''
    AllData.clear()
    Xset.abund = 'feld'
    m_test = Model(mo_str)
    comps = m_test.componentNames
    m_test.show()

    AllData += spec_file_1
    AllData += spec_file_2

    # ignore channels below 0.7 keV and above 8.0 keV
    AllData.ignore("**-0.7 8.0-**")
    AllData.ignore("bad")

    # Manually chosen normalizations (not fit) so each region's unrestricted
    # model is visible on the plot at a comparable scale.
    for comp in comps:
        comp_norm = getattr(m_test, comp).norm
        if "core" in comp:
            comp_norm.values = "1e11,,1e11,1e11,1.6e12,1.6e12"
        elif "ar" in comp:
            comp_norm.values = "1e10,,1e10,1e10,1.6e12,1.6e12"
        elif "bkc" in comp:
            comp_norm.values = "1e11,,1e11,1e11,1.6e12,1.6e12"
        elif "fl_all_fl_wt_avg" in comp:
            comp_norm.values = "1e9,,1e9,1e9,1.6e12,1.6e12"
        elif "fl_x_9_0" in comp:
            comp_norm.values = "1e9,,1e9,1e9,1.6e12,1.6e12"
        comp_norm.frozen = False

    m_test.show()

    Plot.device = '/null'
    Plot.add = True
    Plot.xAxis = "keV"
    Plot('data')

    # spectrum 1 (flaring Sun)
    energy_1 = Plot.x(1)
    energy_err_1 = Plot.xErr(1)
    rates_1 = Plot.y(1)
    rates_err_1 = Plot.yErr(1)

    # spectrum 2 (quiescent Sun)
    energy_2 = Plot.x(2)
    energy_err_2 = Plot.xErr(2)
    rates_2 = Plot.y(2)
    rates_err_2 = Plot.yErr(2)

    fig, ax = plt.subplots()

    ax.errorbar(energy_1, rates_1, xerr=energy_err_1, yerr=rates_err_1,
                c="black", label="Flaring Sun", alpha=0.7, fmt='o', zorder=3, ms=0)
    ax.errorbar(energy_2, rates_2, xerr=energy_err_2, yerr=rates_err_2,
                c="gray", label="Qui. Sun", alpha=0.7, fmt='o', zorder=3, ms=0)

    for i, comp in enumerate(comps):
        if "bkc" in comp:
            comp_label, c = "BKC", "red"
        elif "core" in comp:
            comp_label, c = "CO", "orange"
        elif "ar" in comp:
            comp_label, c = "AR", "purple"
        elif "fl_x_9_0" in comp:
            comp_label, c = "FL X9.0", "green"
        elif "fl_all_fl_wt_avg" in comp:
            comp_label, c = "FL-AVG", "blue"
        else:
            comp_label, c = comp, "black"

        # models should be the same for both spectra, so only plot once against spectrum 1's energy grid
        ax.plot(energy_1, Plot.addComp(i + 1, 1), label=comp_label, color=c, linewidth=2, zorder=4)

    ax.set_xlabel('Energy (keV)', fontsize=16)
    ax.set_ylabel('Counts/sec/keV', fontsize=16)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(ymin=1e-4)
    ax.set_xlim(xmin=0.7, xmax=8.0)
    ax.legend(fontsize=14, loc='best')

    plt.tight_layout()

    ax.xaxis.set_major_formatter(FuncFormatter(format_func))
    ax.xaxis.set_minor_formatter(FuncFormatter(format_func))
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='x', which='minor', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)

    plt.show()
    plt.close("all")


spec_file_1 = "data/spectra/minxss_fm3_PHA_1hr_avg127.pha"  # flaring Sun
spec_file_2 = "data/spectra/minxss_fm3_PHA_1hr_avg167.pha"  # quiescent Sun
mo_str = "bkc+ar+core+fl_all_fl_wt_avg_red_full"

plot_spectra_with_unrestricted_models(spec_file_1, spec_file_2, mo_str, AllData=AllData)