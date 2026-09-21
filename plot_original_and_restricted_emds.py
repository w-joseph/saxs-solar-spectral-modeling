'''
Produces three of the paper's EMD figures, all from the original (raw,
unrestricted) per-region EMD CSVs:

  emd_old()       -> Fig. 1  (fig01_original_emds.pdf) - EMDs of BKC, AR, CO, FL-AVG
  fl_combined()    -> Fig. 2  (fig02_flare_emds_and_average.pdf) - individual flare
                     EMDs (C5.8 through X9.0) alongside their weighted average (FL-AVG)
  emd_old_vs_red() -> Appendix Fig. C.1 (appendix_c1_emd_restriction_comparison.pdf) - original EMDs (faded)
                     vs. restricted EMDs (opaque), for BKC/AR/CO/FL-AVG

NOTE: figures are saved manually, not via savefig() — set ACTIVE_PLOT below,
run, and save the resulting plot to paper_plots/ under the filename noted above.
'''

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def fl_combined():
    '''Fig. 2: individual flare EMDs plus their weighted average (FL-AVG).'''
    flc58 = pd.read_csv('data/emds/emperarea_C5.8_multi.csv')
    flm10 = pd.read_csv('data/emds/emperarea_M1.0_multi.csv')
    flm11 = pd.read_csv('data/emds/emperarea_M1.1_multi.csv')
    flm28 = pd.read_csv('data/emds/emperarea_M2.8_multi.csv')
    flm42 = pd.read_csv('data/emds/emperarea_M4.2_multi.csv')
    flm76 = pd.read_csv('data/emds/emperarea_M7.6_multi.csv')
    flx15 = pd.read_csv('data/emds/emperarea_X1.5_multi.csv')
    flx90 = pd.read_csv('data/emds/emperarea_X9.0_multi.csv')
    fl_all = pd.read_csv('data/emds/norm_FL_all_wt_avg_new.csv')

    fig, ax = plt.subplots()
    plt.step(flc58['temp (K)'], flc58['norm_log_avg (1e44 cm^-5)'], label='FL C5.8', color='violet', where='mid', alpha=0.6, linewidth=5)
    plt.step(flm10['temp (K)'], flm10['norm_log_avg (1e44 cm^-5)'], label='FL M1.0', color='indigo', where='mid', alpha=0.6, linewidth=5)
    plt.step(flm11['temp (K)'], flm11['norm_log_avg (1e44 cm^-5)'], label='FL M1.1', color='black', where='mid', alpha=0.6, linewidth=5)
    plt.step(flm28['temp (K)'], flm28['norm_log_avg (1e44 cm^-5)'], label='FL M2.8', color='green', where='mid', alpha=0.6, linewidth=5)
    plt.step(flm42['temp (K)'], flm42['norm_log_avg (1e44 cm^-5)'], label='FL M4.2', color='teal', where='mid', alpha=0.6, linewidth=5)
    plt.step(flm76['temp (K)'], flm76['norm_log_avg (1e44 cm^-5)'], label='FL M7.6', color='orange', where='mid', alpha=0.6, linewidth=5)
    plt.step(flx15['temp (K)'], flx15['norm_log_avg (1e44 cm^-5)'], label='FL X1.5', color='red', where='mid', alpha=0.6, linewidth=5)
    plt.step(flx90['temp (K)'], flx90['norm_log_avg (1e44 cm^-5)'], label='FL X9.0', color='brown', where='mid', alpha=0.6, linewidth=5)
    plt.step(fl_all['temp (K)'], fl_all['norm (1e44 cm^-5)'], label='FL-AVG', color='blue', where='mid', alpha=1, linewidth=5)

    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_xlabel('Temperature (K)', fontsize=16)
    ax.set_ylabel('EM/area ($10^{44}$ cm$^{-5}$)', fontsize=16)
    ax.set_title('Flare EMDs', fontsize=18)
    ax.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.show()
    plt.close()


def restricted_emd(df, iter=4, e_low_cut=0, e_up_cut=1e10):
    '''
    Restrict the EMD to only include EM bins within `iter` orders of magnitude
    of the max norm (all other bins zeroed out). Also supports an optional
    temperature-range cut (in K), off by default.

    :param df: input dataframe of a coronal region's EMD
    :param iter: orders of magnitude below the peak EM to keep (default 4 = effectively unrestricted)
    :param e_low_cut: lower temperature cut in K (default 0 = no cut)
    :param e_up_cut: upper temperature cut in K (default 1e10 = no cut)
    :return: dataframe with restricted EM bins
    '''
    df_used = df.copy()
    if 'norm_log_avg (1e44 cm^-5)' in df_used.columns:
        df_used = df_used.rename(columns={'norm_log_avg (1e44 cm^-5)': 'norm (1e44 cm^-5)'})

    max_em = df_used['norm (1e44 cm^-5)'].max()
    lower_bound = max_em / 10 ** iter
    df_used['norm (1e44 cm^-5)'] = df_used['norm (1e44 cm^-5)'].where(df_used['norm (1e44 cm^-5)'] >= lower_bound, 0.0)

    df_used['norm (1e44 cm^-5)'] = df_used['norm (1e44 cm^-5)'].where(
        (df_used['temp (K)'] >= e_low_cut) & (df_used['temp (K)'] <= e_up_cut), 0.0)
    return df_used


def emd_old_vs_red():
    '''Appendix Fig. C.1: original EMDs (faded) vs. restricted EMDs (opaque).'''
    fl_wt_avg_old = pd.read_csv('data/emds/norm_FL_all_wt_avg_new.csv')
    core_old = pd.read_csv('data/emds/norm_core_new.csv')
    bkc_old = pd.read_csv('data/emds/norm_BKC_new.csv')
    ar_old = pd.read_csv('data/emds/norm_AR_new.csv')

    fl_wt_avg = restricted_emd(fl_wt_avg_old, iter=1)
    core = restricted_emd(core_old, iter=1)
    bkc = restricted_emd(bkc_old, iter=1)
    ar = restricted_emd(ar_old, iter=1)

    fig, ax = plt.subplots()
    plt.step(core['temp (K)'], core['norm (1e44 cm^-5)'], label='CO', color='orange', where='mid', alpha=1, linewidth=5)
    plt.step(ar['temp (K)'], ar['norm (1e44 cm^-5)'], label='AR', color='purple', where='mid', alpha=1, linewidth=5)
    plt.step(bkc['temp (K)'], bkc['norm (1e44 cm^-5)'], label='BKC', color='red', where='mid', alpha=1, linewidth=5)
    plt.step(fl_wt_avg['temp (K)'], fl_wt_avg['norm (1e44 cm^-5)'], label='FL-AVG', color='blue', where='mid', alpha=1, linewidth=5)
    # original (unrestricted), faded
    plt.step(core_old['temp (K)'], core_old['norm (1e44 cm^-5)'], color='orange', where='mid', alpha=0.2, linewidth=5)
    plt.step(ar_old['temp (K)'], ar_old['norm (1e44 cm^-5)'], color='purple', where='mid', alpha=0.2, linewidth=5)
    plt.step(bkc_old['temp (K)'], bkc_old['norm (1e44 cm^-5)'], color='red', where='mid', alpha=0.2, linewidth=5)
    plt.step(fl_wt_avg_old['temp (K)'], fl_wt_avg_old['norm (1e44 cm^-5)'], color='blue', where='mid', alpha=0.2, linewidth=5)

    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_xlabel('Temperature (K)', fontsize=16)
    ax.set_ylabel('EM/area ($10^{44}$ cm$^{-5}$)', fontsize=16)
    ax.set_title('Restricted EMDs of coronal region types', fontsize=16)
    ax.legend(loc='upper left', fontsize=14)
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.show()
    plt.close()


def emd_old():
    '''Fig. 1: the original (unrestricted) EMDs of BKC, AR, CO, and FL-AVG.'''
    fl_wt_avg_old = pd.read_csv('data/emds/norm_FL_all_wt_avg_new.csv')
    core_old = pd.read_csv('data/emds/norm_core_new.csv')
    bkc_old = pd.read_csv('data/emds/norm_BKC_new.csv')
    ar_old = pd.read_csv('data/emds/norm_AR_new.csv')

    fig, ax = plt.subplots()
    plt.step(core_old['temp (K)'], core_old['norm (1e44 cm^-5)'], color='orange', where='mid', alpha=1, linewidth=5, label='CO')
    plt.step(ar_old['temp (K)'], ar_old['norm (1e44 cm^-5)'], color='purple', where='mid', alpha=1, linewidth=5, label='AR')
    plt.step(bkc_old['temp (K)'], bkc_old['norm (1e44 cm^-5)'], color='red', where='mid', alpha=1, linewidth=5, label='BKC')
    plt.step(fl_wt_avg_old['temp (K)'], fl_wt_avg_old['norm (1e44 cm^-5)'], color='blue', where='mid', alpha=1, linewidth=5, label='FL-AVG')

    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_xlabel('Temperature (K)', fontsize=16)
    ax.set_ylabel('EM/area ($10^{44}$ cm$^{-5}$)', fontsize=16)
    ax.set_title('EMDs of various coronal region types', fontsize=18)
    ax.legend(loc='upper left', fontsize=14)
    plt.ylim(3e-23, 1e-14)
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.show()
    plt.close()


# Change ONLY this line to switch which figure gets produced: "fig1", "fig2", or "appendix_c1"
ACTIVE_PLOT = "appendix_c1"

if ACTIVE_PLOT == "fig1":
    emd_old()
elif ACTIVE_PLOT == "fig2":
    fl_combined()
elif ACTIVE_PLOT == "appendix_c1":
    emd_old_vs_red()