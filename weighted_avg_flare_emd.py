'''
Creates the weighted-average flare emission measure distribution,
norm_FL_all_wt_avg_new.csv — the "FL-AVG" model used throughout
the paper's spectral fits. This is the cross-class weighting step: it combines
the eight already within-class-averaged per-flare-class EMDs (C5.8 through
X9.0, each already averaged in log-space across multiple observed instances
of that flare class — see the paper's Sect. 1, "time-averaged EMDs") into one
weighted average, using the slope of the flare-frequency distribution (FFD)
as the weighting.

Input files (provided directly as data, already within-class-averaged, not
derived elsewhere in this repo):
emperarea_C5.8_multi.csv, emperarea_M1.0_multi.csv, emperarea_M1.1_multi.csv,
emperarea_M2.8_multi.csv, emperarea_M4.2_multi.csv, emperarea_M7.6_multi.csv,
emperarea_X1.5_multi.csv, emperarea_X9.0_multi.csv
'''

import pandas as pd
import numpy as np

def weighted_avg_flare_em_dist_all(csv_list, output_file, alpha=-1.54):
    '''weighted average emission measure distribution of flares.
    weights based on flare frequency distribution given by formula between flare energy and frequency:
    weight_flare=N_flare/N_ref_flare=(E_flare/E_ref_flare)^alpha, where E_ref_flare is the reference flare energy,
    here that would be lowest energy flare, C5.8., and the rest of the flares are M1.0,M1.1,M2.8,M4.2,M7.6,X1.5 and X9.0
    flare energy from x-ray flux, flare energy empirical relation in Stelzer 2022:
    log(E_flare)= a + b · log(goes_flux_flare) , b = 1.150 ± 0.005 and a = 34.711 ± 0.007
    Here I use the Namekata2017 values instead, a=33.67 and b=0.87

    alpha=-1.54 (Aschwanden & Parnell 2002) is the value used in the paper.
    Other candidate slopes considered at the time, kept here for reference:
    alpha ranges from -0.44 to -1.45 (Ram 2025), Ram 2025 uses -1.69, Orlando 2017 uses -1.53.
    '''

    #flare classes and corresponding flare values
    flare_classes = ["A","B","C","M","X"]
    flare_values = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]  # Corresponding to A, B, C, M, X flares
    goes_flux_flare_ref = 5.8* 1e-6  # C5.8 flare as reference, in W/m^2
    a = 33.67
    b = 0.87 #Beate had these Namekata fit pars in her code
    flare_energy_ref = 10 ** (a + b * np.log10(goes_flux_flare_ref))  # Reference flare energy in ergs

    csv_list=csv_list
    df = pd.read_csv(csv_list[0])
    df = df.rename(columns={'norm_log_avg (1e44 cm^-5)': 'norm C5.8'})
    # keep only the temp and norm_log_avg columns
    df = df[['temp (K)', 'norm C5.8']]
    df['weight_C5.8'] = 1  # Weight for the reference flare is always 1
    for i in range(len(csv_list)-1):
        #     # Read the CSV file
        # Merge the DataFrames on column 'temp (K)'
        df2 = pd.read_csv(csv_list[i+1])
        # keep only the temp and norm_log_avg columns
        df2 = df2[['temp (K)', 'norm_log_avg (1e44 cm^-5)']]

        # get flare type from file name, e.g, M1.0 from emperarea_M1.0_multi.csv
        # NOTE: relies on the flare-type segment being exactly 4 characters
        # (e.g. "M1.0", "X9.0") when the filename is split on '_' — matches all
        # current input filenames, but would need adjusting for a differently
        # named file.
        flare_type = [segment for segment in csv_list[i+1].split('_') if len(segment)==4][0][0:]  # Extract the part after 'FL'
        flare_type_class = flare_type[0]  # First character is the flare class (A, B, C, M, X)
        flare_type_value = float(flare_type[1:])  # The rest is the subtype (e.g., 1.0, 2.8, etc.)

        # convert flare type to goes_flux_flare
        # get goes_flux_flare from flare type: [A,B,C,M,X]=[1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
        # subtype comes after the letter: e.g. m1.1= M*1.1= 1e-5*1.1
        goes_flux_flare = flare_values[flare_classes.index(flare_type_class)] * flare_type_value

        #convert goes_flux_flare to flare energy using log(E_flare)= a + b · log(goes_flux_flare) , b = 1.150 ± 0.005 and a = 34.711 ± 0.007
        flare_energy = 10 ** (a + b * np.log10(goes_flux_flare))

        # calculate weight_flare=N_flare/N_ref_flare=(E_flare/E_ref_flare)^alpha, where E_ref_flare is the reference flare energy
        weight_flare= (flare_energy / flare_energy_ref) ** alpha

        #store the weight in the dataframe as weight_<flare_type>
        df2['weight_' + flare_type] = weight_flare

        # replace nans in norm_log_avg with 0
        df2['norm_log_avg (1e44 cm^-5)'] = df2['norm_log_avg (1e44 cm^-5)'].fillna(0)
        df2 = df2.rename(columns={'norm_log_avg (1e44 cm^-5)': 'norm '+flare_type})
        df = pd.merge(df, df2, on='temp (K)')

    # find the weighted mean norm
    # for each row, multiply each norm by its corresponding weight and then sum the results, divide by the sum of the weights
    norm_columns = [col for col in df.columns if col.startswith('norm ')]
    weight_columns = [col for col in df.columns if col.startswith('weight_')]
    df['weighted_avg_norm'] = df[norm_columns].multiply(df[weight_columns].values, axis=0).sum(axis=1) / df[weight_columns].sum(axis=1)

    # keep the columns temp (K), avg_norm
    df = df[['temp (K)', 'weighted_avg_norm']]
    df = df.rename(columns={'weighted_avg_norm': 'norm (1e44 cm^-5)'})
    # replace all NaN values with 0
    df = df.fillna(0)
    # save to csv
    df.to_csv(output_file, index=False)

# Produces the FL-AVG dataset used throughout the paper, with alpha=-1.54
# (Aschwanden & Parnell 2002) — the value cited in the paper's methodology.
weighted_avg_flare_em_dist_all(
    ["data/emds/emperarea_C5.8_multi.csv","data/emds/emperarea_M1.0_multi.csv","data/emds/emperarea_M1.1_multi.csv",
     "data/emds/emperarea_M2.8_multi.csv","data/emds/emperarea_M4.2_multi.csv","data/emds/emperarea_M7.6_multi.csv",
     "data/emds/emperarea_X1.5_multi.csv","data/emds/emperarea_X9.0_multi.csv"],
    "data/emds/norm_FL_all_wt_avg_new.csv",
    alpha=-1.54,
)