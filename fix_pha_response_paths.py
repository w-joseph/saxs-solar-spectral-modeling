'''
One-off fix: updates the RESPFILE/ANCRFILE header keywords in both PHA files
so they point into data/spectra/ instead of a bare filename — needed once the
PHA and response files were moved into that subfolder, since PyXspec resolves
these header paths relative to the working directory the script is run from
(repo root), not relative to the PHA file's own location.

Run this once, from the repo root, after moving the files into data/spectra/.
'''
from astropy.io import fits

pha_files = [
    "data/spectra/minxss_fm3_PHA_1hr_avg127.pha",
    "data/spectra/minxss_fm3_PHA_1hr_avg167.pha",
]

for pha_file in pha_files:
    with fits.open(pha_file, mode='update') as hdul:
        hdul['SPECTRUM'].header['RESPFILE'] = 'data/spectra/minxss_fm3_RMF.fits'
        hdul['SPECTRUM'].header['ANCRFILE'] = 'data/spectra/minxss_fm3_ARF.fits'
        hdul.flush()
    print(f"Updated {pha_file}")
