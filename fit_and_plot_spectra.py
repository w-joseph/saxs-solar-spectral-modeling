import matplotlib.pyplot as plt
from xspec import *
from coronal_regions import *
from matplotlib.ticker import FuncFormatter

# Custom formatter function
def format_func(x,pos):
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

def spec_fit_and_plot_var_comps(spec_file, spec_lc_path, mo_str, ignore_str, AllData=AllData):
    '''
    Performs spectral fitting for given model components
    Parameters:
    spec_file: filename of the spectrum to be fitted
    mo_str: model string defining the components to be used in the fit
    ignore_str: XSPEC ignore-range string (e.g. "**-0.7 3.0-**") — the energy
        range actually fit differs between quiescent (0.7-3.0 keV) and
        flaring (0.7-8.0 keV) spectra; see CONFIGS below.
    AllData: PyXspec AllData container to hold the spectrum data
    Returns:
    norm_values: list of normalization values for each component after fitting
    Plot: PyXspec Plot object containing the fit results
    comps: list of component names used in the model
    fit_stat: fit statistic value after fitting
    fit_dof: degrees of freedom in the fit
    '''
    AllData.clear()
    AllModels.clear()
    Xset.abund = 'feld'

    m_test = Model(mo_str)
    comps = m_test.componentNames
    m_test.show()

    AllData += spec_file
    s1 = AllData(1)

    AllData.ignore(ignore_str)
    AllData.ignore("bad")

    for comp in comps:
        comp_norm = getattr(m_test, comp).norm
        if "fl" in comp:
            comp_norm.values = "1e7,,1e7,1e7,1.6e12,1.6e12"
        elif "ar" in comp:
            comp_norm.values = "1e6,,1e6,1e6,1.6e12,1.6e12"
        elif "bkc" in comp:
            comp_norm.values = "1e6,,1e6,1e6,1.6e12,1.6e12"
        else:
            comp_norm.values = "1e6,,1e6,1e6,1.6e12,1.6e12"
        comp_norm.frozen = False

    for comp in comps:
        if "vnei" in comp:
            comp_tau = getattr(m_test, comp).tau
            comp_tau.values = "1e12,,1e8,1e8,1e13,1e13"
        if ("core" in comp or "fl" in comp or "ar" in comp) and "vabund" in comp and "vnei" not in comp and "norm" not in comp:
            comp_Mg = getattr(m_test, comp).Mg
            comp_Mg.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_Si = getattr(m_test, comp).Si
            comp_Si.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_S = getattr(m_test, comp).S
            comp_S.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_Ar = getattr(m_test, comp).Ar
            comp_Ar.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_Ca = getattr(m_test, comp).Ca
            comp_Ca.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_Fe =getattr(m_test, comp).Fe
            comp_Fe.values = "1.0,,0.8,0.8,1.2,1.2"
        if "bkc" in comp and "vabund" in comp and "vnei" not in comp and "norm" not in comp:
            comp_Mg = getattr(m_test, comp).Mg
            comp_Mg.values = "1.0,,0.8,0.8,1.2,1.2"
            comp_Si = getattr(m_test, comp).Si
            comp_Si.values = "1.0,,0.8,0.8,1.2,1.2"

    for elem in ["Mg", "Si", "S", "Ar", "Ca", "Fe"]:
        first_par = None
        for comp_name in comps:
            comp_obj = getattr(m_test, comp_name)
            if hasattr(comp_obj, elem):
                par = getattr(comp_obj, elem)
                if first_par is None:
                    first_par = par
                else:
                    par.link = first_par

    m_test.show()

    Fit.renorm('auto')
    Fit.nIterations = 500
    Fit.perform()

    norm_values = [getattr(m_test, comp).norm.values[0] for comp in comps]
    fit_stat = Fit.statistic
    fit_dof = Fit.dof

    m_test.show()
    print("filling factors:", comps,fit_stat/fit_dof)

    return norm_values, Plot, comps, fit_stat, fit_dof

def plot_spectra_res(Plot, comps, spec_file, xmax):
    '''
    Plots the spectra and residuals using Matplotlib
    Parameters:
    :param Plot: plot object from PyXspec containing the fit results
    :param comps: component names used in the model
    :param xmax: upper energy limit (keV) for the plot's x-axis — must match
        the upper limit actually used in the fit (see CONFIGS)
    '''
    Plot.device = '/null'
    Plot.add = True
    Plot.xAxis = "keV"
    Plot('data')

    energy = Plot.x()
    energy_err = Plot.xErr()
    rates = Plot.y()
    rates_err = Plot.yErr()
    folded = Plot.model()

    fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    axs[0].errorbar(energy, rates, xerr=energy_err, yerr=rates_err, c="black", fmt='o', zorder=3, ms=0)

    for i, comp in enumerate(comps):
        if "bkc" in comp:
            comp_label = "BKC"
            c = "red"
        elif "core" in comp:
            comp_label = "CO"
            c = "orange"
        elif "ar" in comp:
            comp_label = "AR"
            c = "purple"
        elif "fl_x_9_0_red" in comp:
            comp_label = "FL X9.0"
            c = "green"
        elif "fl_all_fl_wt_avg" in comp:
            comp_label = "FL-AVG"
            c = "blue"
        else:
            comp_label = comp
            c = "black"

        if len(comps) > 1:
            axs[0].plot(energy, Plot.addComp(i + 1), label=comp_label, color=c)
        else:
            axs[0].plot(energy, folded, label=comp_label, color=c)

    if len(comps) > 1:
        axs[0].plot(energy, folded, label="model")

    plt.xlabel('Energy (keV)', fontsize=16)
    axs[0].set_ylabel('Counts/sec/KeV', fontsize=16)
    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_ylim(ymin=1e-4)
    axs[0].set_xlim(xmin=0.7, xmax=xmax)
    axs[0].legend(fontsize=14)

    plt.tight_layout()

    Plot("delc")
    res = Plot.y()
    res_err = Plot.yErr()
    axs[1].errorbar(energy, res, xerr=energy_err, yerr=res_err, c="black", fmt='o', zorder=3, ms=0)
    axs[1].axhline(y=0)
    axs[1].set_ylabel("residuals (delc)", fontsize=14)

    plt.setp(axs[0].get_xticklabels(), visible=True)
    for ax in axs:
        ax.xaxis.set_major_formatter(FuncFormatter(format_func))
        ax.xaxis.set_minor_formatter(FuncFormatter(format_func))
        ax.tick_params(axis='x', which='major', labelsize=14)
        ax.tick_params(axis='x', which='minor', labelsize=14)
        ax.tick_params(axis='y', which='major', labelsize=14)

    plt.show()
    plt.close("all")

# Named model configurations for spectral fitting.
# Each entry: spectrum file, list of model strings to fit, and the energy
# range actually used for that fit. Quiescent-Sun fits use 0.7-3.0 keV;
# flaring-Sun fits use 0.7-8.0 keV (per the paper's Sect. 3.1/4.3 and the
# x-axis ranges shown in Figs. 5-8).
# Figures are saved manually (not via savefig()) — mapping to the paper's
# actual figures:
#
#   "quiescent"          -> Fig. 5  (apec-based quiescent Sun: AR alone, then AR+BKC)
#                            saved as ar_em_more_restricted_obs167.pdf (top)
#                            and bkc_ar_em_more_restricted_obs167.pdf (bottom)
#   "quiescent_vabund"    -> Fig. 6  (quiescent Sun, AR+BKC, variable abundances)
#                            saved as bkc_ar_em_more_restricted_obs167_vabund.pdf
#   "flare_apec"          -> Fig. 7  (apec-based flaring Sun: FL-AVG alone, FL-AVG+CO, FL-AVG+CO+AR)
#                            saved as fl_em_more_restricted_obs127.pdf (top),
#                            fl_core_em_more_restricted_obs127.pdf (middle),
#                            ar_em_more_restricted_core_em_more_restricted_fl_em_more_restricted_obs127.pdf (bottom)
#   "flare_vnei"          -> intermediate step only, no published figure
#                            (nonequilibrium ionization tested alone, before adding variable abundances)
#   "flare_vabund"         -> intermediate step only, no published figure
#                            (variable abundances tested alone, before adding vnei)
#   "flare_vabund_vnei"    -> Fig. 8  (final flaring-Sun result: AR+CO+FL-AVG,
#                            variable abundances for AR/CO, nonequilibrium ionization for FL-AVG)
#                            saved as ar_em_more_restricted_core_em_more_restricted_fl_em_more_restricted_obs127_vabund_vnei.pdf
CONFIGS = {
    "quiescent": (
        "data/spectra/minxss_fm3_PHA_1hr_avg167.pha",
        [
            "ar_em_more_restricted",
            "ar_em_more_restricted+bkc_em_more_restricted",
        ],
        "**-0.7 3.0-**",
        3.0,
    ),
    "quiescent_vabund": (
        "data/spectra/minxss_fm3_PHA_1hr_avg167.pha",
        [
            "ar_em_more_restricted_vabund+bkc_em_more_restricted_vabund",
        ],
        "**-0.7 3.0-**",
        3.0,
    ),
    "flare_apec": (
        "data/spectra/minxss_fm3_PHA_1hr_avg127.pha",
        [
            "fl_all_fl_wt_avg_red",
            "core_em_more_restricted+fl_all_fl_wt_avg_red",
            "ar_em_more_restricted+core_em_more_restricted+fl_all_fl_wt_avg_red",
        ],
        "**-0.7 8.0-**",
        8.0,
    ),
    "flare_vnei": (
        "data/spectra/minxss_fm3_PHA_1hr_avg127.pha",
        [
            "fl_all_fl_wt_avg_red_vnei",
            "core_em_more_restricted+fl_all_fl_wt_avg_red_vnei",
            "ar_em_more_restricted+core_em_more_restricted+fl_all_fl_wt_avg_red_vnei",
        ],
        "**-0.7 8.0-**",
        8.0,
    ),
    "flare_vabund": (
        "data/spectra/minxss_fm3_PHA_1hr_avg127.pha",
        [
            "fl_all_fl_wt_avg_red_vabund",
            "core_em_more_restricted_vabund+fl_all_fl_wt_avg_red_vabund",
            "ar_em_more_restricted_vabund+core_em_more_restricted_vabund+fl_all_fl_wt_avg_red_vabund",
        ],
        "**-0.7 8.0-**",
        8.0,
    ),
    "flare_vabund_vnei": (
        "data/spectra/minxss_fm3_PHA_1hr_avg127.pha",
        [
            "ar_em_more_restricted_vabund+core_em_more_restricted_vabund+fl_all_fl_wt_avg_red_vabund_vnei",
        ],
        "**-0.7 8.0-**",
        8.0,
    ),
}

# Change ONLY this line to switch which fit runs — no other edits needed.
ACTIVE_CONFIG = "flare_vabund_vnei"  # Options: "quiescent", "quiescent_vabund", "flare_apec", "flare_vnei", "flare_vabund", "flare_vabund_vnei"

spec_file, mo_str_array, ignore_str, xmax = CONFIGS[ACTIVE_CONFIG]
spec_lc_path = spec_file

for mo_str in mo_str_array:
    norm_values, Plot, comps, fit_stat, fit_dof = spec_fit_and_plot_var_comps(spec_file, spec_lc_path,
                                                                                       mo_str, ignore_str, AllData=AllData)
    plot_spectra_res(Plot, comps, spec_file, xmax)