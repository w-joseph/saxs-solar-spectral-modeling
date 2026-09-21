import csv
from xspec import *
import numpy as np
import pandas as pd

class coronal_class_spectra:
    '''
    Defines one spectral model object for each of the following coronal regions:
    coronal region    | object name     | model name
    --------------------------------------------
    1. Background     | bkc             | "bkc"
    2. Cores          | core            | "core"
    3. Active regions | active_region   | "ar"
    4. Flares         | fl_all_fl_wt_avg| "fl_all_fl_wt_avg_red"
    (wt.avg. flares)  |                 |
    5. C 5.8 flare    | fl_c_5_8        | "fl_c_5_8"
    6. M 1.0 flare    | fl_m_1_0        | "fl_m_1_0"
    7. M 1.1 flare    | fl_m_1_1        | "fl_m_1_1"
    8. M 2.8 flare    | fl_m_2_8        | "fl_m_2_8"
    9. M 4.2 flare    | fl_m_4_2        | "fl_m_4_2"
    10. M 7.6 flare   | fl_m_7_6        | "fl_m_7_6"
    11. X 1.5 flare   | fl_x_1_5        | "fl_x_1_5"
    12. X 9.0 flare   | fl_x_9_0        | "fl_x_9_0"

    Use the standard PyXspec method of creating a model from these components
    e.g., Model("ar")
          Model("core+ar+fl")
    Can also be used with standard XSPEC models:
    e.g., Model("apec+fl_all_fl_wt_avg_red")

    NOTE: this list covers the original, base set of regions only. Many more
    variants (EMD-restricted, variable-abundance, non-equilibrium-ionization
    versions of these same regions, plus the weighted-average flare model)
    are instantiated further down in this file, below the class definition —
    see the "NAMING CONVENTION" comment block there for what each variant's
    suffix means, and that section overall for the full, current list of
    every region object actually available for use in Model(...) calls.'''

    def __init__(self,temp_and_norm_file,mod_name):
        self.temp_and_norm_file = temp_and_norm_file
        self.mod_name=mod_name

    def _load_restricted_emd(self, distance, iter, e_low_cut, e_up_cut):
        '''
        Shared data-loading step used by all model_creation_* methods below.
        Reads this region's EMD CSV, converts temperature/normalization to
        the units XSPEC expects, then applies two restrictions:
          - EMD restriction: bins more than `iter` orders of magnitude below
            the peak emission measure are zeroed out (iter=4, the default
            used everywhere below, is effectively unrestricted; iter=1 is
            the "EMD 1 order restricted" models used throughout the paper).
          - Optional temperature-range cut (e_low_cut/e_up_cut, in eV),
            off by default (0 to 1e10 covers the full range).
        Returns (temp, norm) as numpy arrays, ready for building a model string.
        '''
        distance_cm = distance * 3.086e18

        # Open the CSV file for reading
        with open(self.temp_and_norm_file, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  # skip header row
            temp = []
            norm = []
            for row in csvreader:
                #convert from /cm2 to /km2 because values are too small by multiplying by 1e10
                temp.append(float(row[0]) * 1e-3 / 11606)
                norm.append(float(row[1]) * 1e44 * 1e10 * (1e-14 / (4 * 3.14 * distance_cm ** 2)))

        df = pd.DataFrame({'temp (K)': temp, 'norm (1e44 cm^-5)': norm})

        # EMD restriction: zero out bins more than `iter` orders of magnitude below the peak
        max_em = df['norm (1e44 cm^-5)'].max()
        lower_bound = max_em / 10 ** iter
        df['norm (1e44 cm^-5)'] = df['norm (1e44 cm^-5)'].where(df['norm (1e44 cm^-5)'] >= lower_bound, 0.0)

        # Optional temperature-range cut
        e_low_cut_kev = e_low_cut * 1e-3 / 11606
        e_up_cut_kev = e_up_cut * 1e-3 / 11606
        df['norm (1e44 cm^-5)'] = df['norm (1e44 cm^-5)'].where(
            (df['temp (K)'] >= e_low_cut_kev) & (df['temp (K)'] <= e_up_cut_kev), 0.0)

        return df['temp (K)'].to_numpy(), df['norm (1e44 cm^-5)'].to_numpy()

    def _finalize_model(self, mo):
        '''Shared finishing step used by all model_creation_* methods below: strips
        the trailing "+", stores and prints the assembled model string, and
        registers it with XSPEC via AllModels.mdefine().'''
        mo = mo[:-1]  # remove the extra "+" left over from the build loop
        print("mdef ", self.mod_name, mo)
        self.mo = mo
        AllModels.mdefine(self.mod_name + " " + self.mo)

    def model_creation_vvapec(self,abundance=1.0,redshift=0.0,distance=4.84814e-6,iter=4, e_low_cut=0, e_up_cut=1e10):
        '''
        Function to define the spectral models.
        :abundance: abundance parameter for the vvapec model
        :redshift: redshift parameter for the vvapec model
        :distance: distance to the object in parsec
        :iter: number of orders of magnitude below the maximum emission measure to include
        :e_low_cut: lower temperature cut in eV
        :e_up_cut: upper temperature cut in eV

        distance prox cen= 1.295 pc
        distance AD Leo= 4.965 pc
        distance eta eri =41.9 pc
        distance Sun= 4.84814e-6
        e.g., To set the abundance of the flare component to 1.0:
        flare_low.model_creation_vvapec(abundance=0.1)

        e.g., To set the abundance as a variable that can be fit, send any String value as a parameter say "a" for example:
        flare_low.model_creation_vvapec(abundance="a")

        NOTE: This function sets the abundance and redshift for one component at a time — call it individually on each region object to change multiple at once.'''
        temp, norm = self._load_restricted_emd(distance, iter, e_low_cut, e_up_cut)

        mo = ""
        for i in range(len(temp)):
            #using Schwab 2020 abund from table 4: Mg=1.02, Si=0.99,S=1,Fe=1.35
            #H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            mo += str(norm[i]) + "*vvapec(" + str(
                temp[i]) + ",1,1,1,1,1,1,1,1,1,1,1,1.02,1,0.99,1,1,1,1,1,1,1,1,1,1,1,1.35,1,1,1,1," + str(redshift) + ")+"

        self._finalize_model(mo)

    def model_creation_vvnei(self,abundance=1.0,redshift=0.0,distance=4.84814e-6,iter=4, e_low_cut=0, e_up_cut=1e10):
        '''
        Exactly like model_creation_vvapec but uses the vvnei model instead of vvapec,
        vvnei is vvapec in non-equilibrium ionization and has the extra parameter tau (ionization timescale)
        '''
        temp, norm = self._load_restricted_emd(distance, iter, e_low_cut, e_up_cut)

        mo = ""
        for i in range(len(temp)):
            #using Schwab 2020 abund from table 4: Mg=1.02, Si=0.99,S=1,Fe=1.35 (same as model_creation_vvapec, plus tau)
            #H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            mo += str(norm[i]) + "*vvnei(" + str(
                temp[i]) + ",1,1,1,1,1,1,1,1,1,1,1,1.02,1,0.99,1,1,1,1,1,1,1,1,1,1,1,1.35,1,1,1,1,tau," + str(redshift) + ")+"

        self._finalize_model(mo)

    def model_creation_vvapec_ar_core_fl_abundvar(self,abundance=1.0,redshift=0.0,distance=4.84814e-6,iter=4, e_low_cut=0, e_up_cut=1e10):
        '''
        Similar to model_creation_vvapec but uses variable Mg, Si, S, Ar, Ca, Fe abundances for ar, core and fl components
        on the recommendation of Woods2023, Table 1
        '''
        temp, norm = self._load_restricted_emd(distance, iter, e_low_cut, e_up_cut)

        mo = ""
        for i in range(len(temp)):
            #using Woods2023 abundance free recommendation from table 1: Mg, Si, S, Ar,Ca,Fe
            #H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            mo += str(norm[i]) + "*vvapec(" + str(
                temp[i]) + ",1,1,1,1,1,1,1,1,1,1,1,Mg,1,Si,1,S,1,Ar,1,Ca,1,1,1,1,1,Fe,1,1,1,1," + str(redshift) + ")+"

        self._finalize_model(mo)

    def model_creation_vvnei_abundvar(self,abundance=1.0,redshift=0.0,distance=4.84814e-6,iter=4, e_low_cut=0, e_up_cut=1e10):
        '''
        Similar to model_creation_vvapec but uses variable Mg, Si, S, Ar, Ca, Fe abundances for ar, core and fl components
        on the recommendation of Woods2023, Table 1.
        In addition, uses the vvnei model instead of vvapec,
        vvnei is vvapec in non-equilibrium ionization and has the extra parameter tau (ionization timescale)
        '''
        temp, norm = self._load_restricted_emd(distance, iter, e_low_cut, e_up_cut)

        mo = ""
        for i in range(len(temp)):
            #using Woods2023 abundance free recommendation from table 1: Mg, Si, S, Ar,Ca,Fe (same as model_creation_vvapec_ar_core_fl_abundvar, plus tau)
            #H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            mo += str(norm[i]) + "*vvnei(" + str(
                temp[i]) + ",1,1,1,1,1,1,1,1,1,1,1,Mg,1,Si,1,S,1,Ar,1,Ca,1,1,1,1,1,Fe,1,1,1,1,tau," + str(redshift) + ")+"

        self._finalize_model(mo)

    def model_creation_vvapec_bkc_abundvar(self,abundance=1.0,redshift=0.0,distance=4.84814e-6,iter=4, e_low_cut=0, e_up_cut=1e10):
        '''
        Similar to model_creation_vvapec but uses variable Mg and Si abundances for bkc components on the recommendation of Woods2023, Table 1.
        '''
        temp, norm = self._load_restricted_emd(distance, iter, e_low_cut, e_up_cut)

        mo = ""
        for i in range(len(temp)):
            #using Woods2023 abundance free recommendation from table 1: Mg=variable, Si=variable
            #H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            mo += str(norm[i]) + "*vvapec(" + str(
                temp[i]) + ",1,1,1,1,1,1,1,1,1,1,1,Mg,1,Si,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1," + str(redshift) + ")+"

        self._finalize_model(mo)


#creating all the objects corresponding to the coronal regions
#
# NAMING CONVENTION — every object below is one of the four base regions
# (core, active_region/ar, bkc, or a flare), with a suffix indicating which
# variant it is. Suffixes combine additively in the order listed:
#
#   (no suffix)          full, unrestricted EMD (iter=4, the default) —
#                        equivalent to a plain apec model (all abundances
#                        frozen at 1.0)
#   _e_restricted        temperature-range cut only (e_low_cut/e_up_cut),
#                        EMD itself not restricted
#   _em_restricted        EMD restricted to iter=2 (2 orders of magnitude
#                        below the peak emission measure)
#   _em_more_restricted   EMD restricted to iter=1 — this is the restriction
#                        level used throughout most of the paper's fits and
#                        figures (see Sect. 3.2, "EMD 1 order restricted")
#   _vabund               variable (free-fitting) elemental abundances,
#                        following Woods et al. (2023) Table 1: all six of
#                        Mg/Si/S/Ar/Ca/Fe for core/ar/flares, but only Mg/Si
#                        for bkc (see Sect. 3.3)
#   _vnei                 non-equilibrium ionization (vvnei model, with an
#                        extra free ionization-timescale parameter, tau)
#                        instead of equilibrium (vvapec)
#   _vabund_vnei           both of the above combined — this is the model
#                        used for the paper's best-fit flaring-Sun result
#
# e.g. ar_em_more_restricted_vabund = active regions, EMD restricted to
# iter=1, with variable abundances.
#
# Flares are the one exception to the "iter=1 needs the suffix" rule:
# fl_all_fl_wt_avg_red (no "_em_more_restricted" suffix) is ALREADY built
# with iter=1 directly — there's a separate _full object
# (fl_all_fl_wt_avg_red_full) for the unrestricted version instead.
core=coronal_class_spectra('data/emds/norm_core_new.csv',"core")
core_e_restricted=coronal_class_spectra('data/emds/norm_core_new.csv',"core_e_restricted")
core_em_restricted=coronal_class_spectra('data/emds/norm_core_new.csv',"core_em_restricted")
core_em_more_restricted=coronal_class_spectra('data/emds/norm_core_new.csv',"core_em_more_restricted")
core_em_more_restricted_vabund=coronal_class_spectra('data/emds/norm_core_new.csv',"core_em_more_restricted_vabund")
active_region=coronal_class_spectra('data/emds/norm_AR_new.csv',"ar")
active_region_e_restricted=coronal_class_spectra('data/emds/norm_AR_new.csv',"ar_e_restricted")
active_region_em_restricted=coronal_class_spectra('data/emds/norm_AR_new.csv',"ar_em_restricted")
active_region_em_more_restricted=coronal_class_spectra('data/emds/norm_AR_new.csv',"ar_em_more_restricted")
active_region_em_more_restricted_vabund=coronal_class_spectra('data/emds/norm_AR_new.csv',"ar_em_more_restricted_vabund")
fl_c_5_8=coronal_class_spectra('data/emds/norm_FLC5.8_new.csv',"fl_c_5_8")
fl_m_1_0=coronal_class_spectra('data/emds/norm_FLM1.0_new.csv',"fl_m_1_0")
fl_m_1_1=coronal_class_spectra('data/emds/norm_FLM1.1_new.csv',"fl_m_1_1")
fl_m_2_8=coronal_class_spectra('data/emds/norm_FLM2.8_new.csv',"fl_m_2_8")
fl_m_4_2=coronal_class_spectra('data/emds/norm_FLM4.2_new.csv',"fl_m_4_2")
fl_m_7_6=coronal_class_spectra('data/emds/norm_FLM7.6_new.csv',"fl_m_7_6")
fl_x_1_5=coronal_class_spectra('data/emds/norm_FLX1.5_new.csv',"fl_x_1_5")
fl_x_9_0=coronal_class_spectra('data/emds/norm_FLX9.0_new.csv',"fl_x_9_0")
fl_x_9_0_em_more_restricted=coronal_class_spectra('data/emds/norm_FLX9.0_new.csv',"fl_x_9_0_em_more_restricted")
fl_x_9_0_red=coronal_class_spectra('data/emds/norm_FLX9.0_red.csv',"fl_x_9_0_red")
bkc=coronal_class_spectra('data/emds/norm_BKC_new.csv',"bkc")
bkc_e_restricted=coronal_class_spectra('data/emds/norm_BKC_new.csv',"bkc_e_restricted")
bkc_em_restricted=coronal_class_spectra('data/emds/norm_BKC_new.csv',"bkc_em_restricted")
bkc_em_more_restricted=coronal_class_spectra('data/emds/norm_BKC_new.csv',"bkc_em_more_restricted")
bkc_em_more_restricted_vabund=coronal_class_spectra('data/emds/norm_BKC_new.csv',"bkc_em_more_restricted_vabund")

#Defining all the spectral models of the coronal regions.
core.model_creation_vvapec()
core_e_restricted.model_creation_vvapec(e_up_cut=1e7)
core_em_restricted.model_creation_vvapec(iter=2)
core_em_more_restricted.model_creation_vvapec(iter=1)
core_em_more_restricted_vabund.model_creation_vvapec_ar_core_fl_abundvar(iter=1)
active_region.model_creation_vvapec()
active_region_e_restricted.model_creation_vvapec(e_up_cut=1e7)
active_region_em_restricted.model_creation_vvapec(iter=2)
active_region_em_more_restricted.model_creation_vvapec(iter=1)
active_region_em_more_restricted_vabund.model_creation_vvapec_ar_core_fl_abundvar(iter=1)
fl_c_5_8.model_creation_vvapec()
fl_m_1_0.model_creation_vvapec()
fl_m_1_1.model_creation_vvapec()
fl_m_2_8.model_creation_vvapec()
fl_m_4_2.model_creation_vvapec()
fl_m_7_6.model_creation_vvapec()
fl_x_1_5.model_creation_vvapec()
fl_x_9_0.model_creation_vvapec()
fl_x_9_0_red.model_creation_vvapec()
fl_x_9_0_em_more_restricted.model_creation_vvapec(iter=1)
bkc.model_creation_vvapec()
bkc_e_restricted.model_creation_vvapec(e_low_cut=1e6, e_up_cut=1e7)
bkc_em_restricted.model_creation_vvapec(iter=2)
bkc_em_more_restricted.model_creation_vvapec(iter=1)
bkc_em_more_restricted_vabund.model_creation_vvapec_bkc_abundvar(iter=1)

#wt avg of all flare unrestricted
fl_all_fl_wt_avg_red_full=coronal_class_spectra('data/emds/norm_FL_all_wt_avg_new.csv',"fl_all_fl_wt_avg_red_full")
#wt avg of all flares
fl_all_fl_wt_avg_red=coronal_class_spectra('data/emds/norm_FL_all_wt_avg_new.csv',"fl_all_fl_wt_avg_red")
#wt avg of all flares but with vnei model defined later
fl_all_fl_wt_avg_red_vnei=coronal_class_spectra('data/emds/norm_FL_all_wt_avg_new.csv',"fl_all_fl_wt_avg_red_vnei")
#var abund and vnei
fl_all_fl_wt_avg_red_vabund_vnei=coronal_class_spectra('data/emds/norm_FL_all_wt_avg_new.csv',"fl_all_fl_wt_avg_red_vabund_vnei")
#wt avg of all flares but with var abund
fl_all_fl_wt_avg_red_vabund=coronal_class_spectra('data/emds/norm_FL_all_wt_avg_new.csv',"fl_all_fl_wt_avg_red_vabund")

fl_all_fl_wt_avg_red.model_creation_vvapec(iter=1)
fl_all_fl_wt_avg_red_full.model_creation_vvapec()
#vnei for flares
fl_all_fl_wt_avg_red_vnei.model_creation_vvnei(iter=1)
# var abund for fl on woods 2023 recommendation
fl_all_fl_wt_avg_red_vabund.model_creation_vvapec_ar_core_fl_abundvar(iter=1)
#var abund and vnei
fl_all_fl_wt_avg_red_vabund_vnei.model_creation_vvnei_abundvar(iter=1)