# daxss_pha_spectrum_extraction.py
# Extracts one DAXSS spectrum (selected via SPECTRUM_INDEX/ACTIVE) into a
# standard OGIP Type-I PHA FITS file for use with XSPEC/PyXspec.

import netCDF4 as nc
import numpy as np
from astropy.io import fits

SPECTRUM_INDEX = {"flare": 127, "quiescent": 167}
ACTIVE = "quiescent"  # change this to extract a different spectrum

DAXSS_LEVEL2_FILE = "data/spectra/daxss_solarSXR_level2_1hour_average_2022-02-14-mission_v2.1.0.ncdf"

daxss = nc.Dataset(DAXSS_LEVEL2_FILE)
spectrum_index = SPECTRUM_INDEX[ACTIVE]

hdr_dummy = fits.Header()
hdr_data = fits.Header()
for hdr in (hdr_dummy, hdr_data):
    hdr["MISSION"] = "InspireSat-1"
    hdr["TELESCOP"] = "InspireSat-1"
    hdr["INSTRUME"] = "DAXSS"
    hdr["ORIGIN"] = "LASP"
    hdr["CREATOR"] = "DAXSSPlotterUtility_v1"
hdr_dummy["CONTENT"] = "Type-I PHA file"

hdr_data["CONTENT"] = "SPECTRUM"
hdr_data["HDUCLASS"] = "OGIP"
hdr_data["LONGSTRN"] = "OGIP 1.0"
hdr_data["HDUCLAS1"] = "SPECTRUM"
hdr_data["HDUVERS1"] = "1.2.1"
hdr_data["HDUVERS"] = "1.2.1"
hdr_data["AREASCAL"] = "1"
hdr_data["BACKSCAL"] = "1"
hdr_data["CORRSCAL"] = "1"
hdr_data["BACKFILE"] = "none"
hdr_data['RESPFILE'] = "data/spectra/minxss_fm3_RMF.fits"
hdr_data['ANCRFILE'] = "data/spectra/minxss_fm3_ARF.fits"
hdr_data["CHANTYPE"] = "PHA"
hdr_data["POISSERR"] = "FALSE"  # TODO: confirm this is correct — flagged as "broken" in an earlier version
hdr_data["CORRFILE"] = "none"
hdr_data["EXTNAME"] = "SPECTRUM"
hdr_data["FILTER"] = "Be/Kapton"
hdr_data["EXPOSURE"] = "9"
hdr_data["DETCHANS"] = "1000"
hdr_data["GROUPING"] = "0"

channel_number_array = []
systematic_error_array = []
for i in range(1, 1001):
    channel_number_array.append(np.int32(i))
    systematic_error_array.append(np.float32(
        daxss["SPECTRUM_CPS_ACCURACY"][spectrum_index, i + 5] /
        daxss["SPECTRUM_CPS"][spectrum_index, i + 5]
    ))

c1 = channel_number_array
c2 = daxss["SPECTRUM_CPS"][spectrum_index, 6:1006]
c3 = daxss["SPECTRUM_CPS_PRECISION"][spectrum_index, 6:1006]
c4 = systematic_error_array

time_iso_string = "".join(
    daxss["TIME_ISO"][int(spectrum_index)][i].decode("utf-8") for i in range(20)
)
hdr_dummy["FILENAME"] = f"minxss_fm3_PHA_{time_iso_string.replace(':', '-')}.pha"
hdr_dummy["DATE"] = time_iso_string.replace(":", "-")
hdr_data["FILENAME"] = hdr_dummy["FILENAME"]
hdr_data["DATE"] = hdr_dummy["DATE"]

hdu_data = fits.BinTableHDU.from_columns(
    [fits.Column(name="CHANNEL", format="J", array=c1),
     fits.Column(name="RATE", format="E", array=c2),
     fits.Column(name="STAT_ERR", format="E", array=c3),
     fits.Column(name="SYS_ERR", format="E", array=c4)],
    header=hdr_data,
)
hdul = fits.HDUList([fits.PrimaryHDU(header=hdr_dummy), hdu_data])

output_filename = f"data/spectra/minxss_fm3_PHA_1hr_avg{spectrum_index}.pha"
hdul.writeto(output_filename, overwrite=True)
print(f"Saved: {output_filename}")