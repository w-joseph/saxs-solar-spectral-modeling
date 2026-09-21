'''
This code plots the Hinode solar image corresponding to the Flaring Sun and Quiescent Sun DAXSS spectra.
It also compares the filling factors of coronal regions we get from fitting DAXSS spectra with SaXS spectral models,
with the filling factors obtained independently by Hinode image segmentation done by Adithya2021 on Hinode/XRT maps for AR,
BP, and limb regions, and the HEK (multi-wavelength solar database) for AR and FL regions and plots these region comparisons.
The Hinode solar image and Hinode segmentation maps need to be downloaded separately and stored in the same directory as this script.
The HEK data is queried live from the HEK database using the date and time of the Hinode map.
'''

import sunpy.map
from sunpy.net import Fido, attrs as a
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from xrtpy.image_correction.remove_lightleak import remove_lightleak
from matplotlib.colors import LogNorm
from sunpy.net import hek
from shapely import wkt
from astropy.coordinates import SkyCoord
from matplotlib.path import Path

# XRT segmentation contour levels shared by both the flare and quiescent maps
XRT_SEG_CONTOURS = [
    {"level": 4, "color": "teal", "label": "Limb (XRT seg.)"},
    {"level": 16, "color": "red", "label": "BP (XRT seg.)"},
    {"level": 32, "color": "black", "label": "AR (XRT seg.)"},
]

# Approx. total pixels on the solar disk for the Hinode maps used here,
# per V1_XRT_Irradiance_Data_20080201T181156.9_20230617T181910.2.csv
TOTAL_DISK_PIXELS = 616473.0


def add_hek_overlays_and_finalize(ax, fig, hinode_map_wo_light_leak, additional_map_filename, ars, fls, title):
    '''
    Shared finishing steps for both the flare and quiescent Hinode plots:
    - draws XRT segmentation contours (limb, bright points, active regions)
    - overlays HEK-catalog active-region and flare bounding boxes
    - prints total pixel counts and bbox-derived filling factors for AR/FL
    - assembles the shared legend, axis labels, tick sizes, and title
    - shows the finished figure

    Region-specific plotting (the colored filling-factor overlays themselves)
    happens in the caller before this is invoked.
    '''
    additional_map = sunpy.map.Map(additional_map_filename)

    for c in XRT_SEG_CONTOURS:
        mask = (additional_map.data == c["level"])
        plt.contour(mask.astype(float), levels=[0.5],
                    colors=[c["color"]], linewidths=1)
        plt.plot([], [], color=c["color"], label=c["label"])

    flag_ar = False
    flag_fl = False
    tot_pixels_ar = 0
    tot_pixels_fl = 0

    for ar in ars:
        print(ar['obs_observatory'])
        if ar["hpc_bbox"] != '':
            poly = wkt.loads(ar["hpc_bbox"])
            x, y = poly.exterior.xy

            # Create coordinates with HEK's observer (usually Earth/SDO), then
            # transform to Hinode's observer location
            coords_hek = SkyCoord(
                x * u.arcsec, y * u.arcsec,
                frame="helioprojective",
                obstime=hinode_map_wo_light_leak.date,
                observer="earth",
            )
            coords_hinode = coords_hek.transform_to(hinode_map_wo_light_leak.coordinate_frame)

            x_pix, y_pix = hinode_map_wo_light_leak.world_to_pixel(coords_hinode)
            verts = np.vstack((x_pix.value, y_pix.value)).T
            poly_path = Path(verts)

            ny, nx = hinode_map_wo_light_leak.data.shape
            Y, X = np.mgrid[:ny, :nx]
            points = np.vstack((X.ravel(), Y.ravel())).T
            mask = poly_path.contains_points(points).reshape(hinode_map_wo_light_leak.data.shape)

            region_data = np.ma.array(hinode_map_wo_light_leak.data, mask=~mask)
            mean_dn = region_data.mean()
            n_pix = mask.sum()
            tot_pixels_ar += n_pix
            print(f"Region mean DN/s = {mean_dn:.2f}, area (pixels) = {n_pix}")

            if not flag_ar:
                ax.plot_coord(coords_hinode, color="black", linewidth=1, label="AR (HEK)")
                flag_ar = True
            else:
                ax.plot_coord(coords_hinode, color="black", linewidth=1)

    for fl in fls:
        print(fl['obs_observatory'])
        if fl["hpc_bbox"] != '':
            poly = wkt.loads(fl["hpc_bbox"])
            x, y = poly.exterior.xy

            coords_hek = SkyCoord(
                x * u.arcsec, y * u.arcsec,
                frame="helioprojective",
                obstime=hinode_map_wo_light_leak.date,
                observer="earth",
            )
            coords_hinode = coords_hek.transform_to(hinode_map_wo_light_leak.coordinate_frame)

            x_pix, y_pix = hinode_map_wo_light_leak.world_to_pixel(coords_hinode)
            verts = np.vstack((x_pix.value, y_pix.value)).T
            poly_path = Path(verts)

            ny, nx = hinode_map_wo_light_leak.data.shape
            Y, X = np.mgrid[:ny, :nx]
            points = np.vstack((X.ravel(), Y.ravel())).T
            mask = poly_path.contains_points(points).reshape(hinode_map_wo_light_leak.data.shape)

            region_data = np.ma.array(hinode_map_wo_light_leak.data, mask=~mask)
            mean_dn = region_data.mean()
            n_pix = mask.sum()
            tot_pixels_fl += n_pix
            print(f"Region mean DN/s = {mean_dn:.2f}, area (pixels) = {n_pix}")

            if n_pix > 0:
                if not flag_fl:
                    ax.plot_coord(coords_hinode, color="blue", linewidth=1, label="FL (HEK)")
                    flag_fl = True
                else:
                    ax.plot_coord(coords_hinode, color="blue", linewidth=1)

    print("Total AR pixels from bbox:", tot_pixels_ar)
    print("Total FL pixels from bbox:", tot_pixels_fl)
    print("filling factors AR from bbox:", tot_pixels_ar / TOTAL_DISK_PIXELS * 100)
    print("filling factors FL from bbox:", tot_pixels_fl / TOTAL_DISK_PIXELS * 100)

    plt.tight_layout()

    handles, labels = ax.get_legend_handles_labels()
    fig.subplots_adjust(bottom=0.1)
    fig.legend(handles, labels, loc='lower center', ncol=3,
               bbox_to_anchor=(0.5, -0.01), fancybox=True, shadow=True, fontsize=14)

    ax.set_xlabel("Helioprojective Longitude (Solar-X)", fontsize=14)
    ax.set_ylabel("Helioprojective Latitude (Solar-Y)", fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)

    for cb_ax in fig.findobj(plt.Axes):
        if cb_ax != ax:
            cb_ax.tick_params(axis='both', which='major', labelsize=12)
            if cb_ax.yaxis.label.get_text() != '':
                cb_ax.yaxis.label.set_size(14)
            if cb_ax.xaxis.label.get_text() != '':
                cb_ax.xaxis.label.set_size(14)

    plt.title(title, fontsize=16)
    plt.show()


def plot_base_map(hinode_map_wo_light_leak, colorbar_ticks):
    '''Plots the raw Hinode/XRT map on its own, before any region overlays.'''
    hinode_map_wo_light_leak.plot(title="Hinode/XRT (DN/s)")
    plt.title("Hinode/XRT (DN/s)", fontsize=16)
    plt.tick_params(labelsize=12)
    plt.ylabel("Helioprojective Latitude (Solar-Y)", fontsize=14)
    plt.xlabel("Helioprojective Longitude (Solar-X)", fontsize=14)
    cb = plt.colorbar()
    cb.set_ticks(colorbar_ticks)
    cb.ax.tick_params(labelsize=12)


def comp_map_regs_w_contours(hinode_map_wo_light_leak, ff_1, ff_2, ff_3, ars, fls):
    '''Flare observation: three filling-factor regions (AR, CO, FL-AVG).'''
    p_ff_1 = 100 - ff_1 if ff_1 is not None else 100
    p_ff_2 = 100 - ff_2 if ff_2 is not None else 100
    p_ff_3 = 100 - ff_3 if ff_3 is not None else 100

    p_x = np.percentile(hinode_map_wo_light_leak.data, p_ff_1)
    p_y = np.percentile(hinode_map_wo_light_leak.data, p_ff_2)
    p_z = np.percentile(hinode_map_wo_light_leak.data, p_ff_3)
    print(p_x, p_y, p_z, hinode_map_wo_light_leak.data.max(), hinode_map_wo_light_leak.data.min())

    bright_mask_1 = (hinode_map_wo_light_leak.data >= p_x)
    bright_mask_2 = (hinode_map_wo_light_leak.data >= p_y) & (hinode_map_wo_light_leak.data < p_x)
    bright_mask_3 = (hinode_map_wo_light_leak.data >= p_z) & (hinode_map_wo_light_leak.data < p_y)

    masked_data_1 = np.ma.array(hinode_map_wo_light_leak.data, mask=~bright_mask_1)
    masked_data_2 = np.ma.array(hinode_map_wo_light_leak.data, mask=~bright_mask_2)
    masked_data_3 = np.ma.array(hinode_map_wo_light_leak.data, mask=~bright_mask_3)

    plot_base_map(hinode_map_wo_light_leak, colorbar_ticks=[0, 1000, 10000, 100000])

    fig = plt.figure(figsize=(8, 5.2))
    ax = plt.subplot(projection=hinode_map_wo_light_leak)

    plt.imshow(masked_data_3, origin="lower", norm=LogNorm(vmin=p_z, vmax=p_y), cmap="summer", interpolation='none')
    cb = plt.colorbar(shrink=0.75)
    cb.set_label('DN/s per pixel, AR')

    plt.imshow(masked_data_2, origin="lower", norm=LogNorm(vmin=p_y, vmax=p_x), cmap='autumn', interpolation='none')
    cb_1 = plt.colorbar(shrink=0.75, pad=0.1)
    cb_1.set_label('DN/s per pixel, CO')

    plt.imshow(masked_data_1, origin="lower", norm=LogNorm(vmin=p_x), cmap='winter', interpolation='none')
    cb_2 = plt.colorbar(shrink=0.75)
    cb_2.set_label('DN/s per pixel, FL AVG')

    add_hek_overlays_and_finalize(
        ax, fig, hinode_map_wo_light_leak,
        additional_map_filename='data/hinode/comp_XRT20220425_054836.5_seg.fits',
        ars=ars, fls=fls, title="AR, CO and FL AVG",
    )


def comp_map_regs_qui_w_contours(hinode_map_wo_light_leak, ff_1, ars, fls):
    '''Quiescent observation: one filling-factor region (AR).'''
    p_ff_1 = 100 - ff_1 if ff_1 is not None else 100
    p_x = np.percentile(hinode_map_wo_light_leak.data, p_ff_1)
    bright_mask_1 = (hinode_map_wo_light_leak.data >= p_x)
    masked_data_1 = np.ma.array(hinode_map_wo_light_leak.data, mask=~bright_mask_1)

    plot_base_map(hinode_map_wo_light_leak, colorbar_ticks=[0, 1e2, 1e3, 1e4, 1e5])
    plt.clim(vmin=0)

    fig = plt.figure(figsize=(7, 6.5))
    ax = plt.subplot(projection=hinode_map_wo_light_leak)

    plt.imshow(masked_data_1, origin="lower", norm=LogNorm(vmin=p_x), cmap='summer', interpolation='none')
    cb_1 = plt.colorbar(shrink=0.85)
    cb_1.set_label('DN/s per pixel, AR')

    add_hek_overlays_and_finalize(
        ax, fig, hinode_map_wo_light_leak,
        additional_map_filename='data/hinode/comp_XRT20220629_180339.8_seg.fits',
        ars=ars, fls=fls, title="AR",
    )


# NOTE: figures from this script were fine-tuned and saved manually (not via
# savefig()). Set obs to "flare" or "qui" below, run, and save the resulting
# plot. Based on paper_plots/, this generally corresponds to:
#   obs="flare" -> hinode_flare_ar_co_fl_hek_contours_vabund_vnei.pdf (and the _zoom variant)
#   obs="qui"   -> low_cts_ar_hek_contours_vabund.pdf
obs = "qui"

if obs == "flare":
    # highest counts — Be_thin filter, no light-leak correction needed (see paper Sect. 3.4)
    hinode_map = sunpy.map.Map('data/hinode/comp_XRT20220425_020655.1.fits')
    hinode_map_wo_light_leak = hinode_map
    date_start = "2022-04-25 02:01"
    date_end = "2022-04-25 02:11"
elif obs == "qui":
    # lowest counts — Al_mesh filter, needs light-leak correction
    hinode_map = sunpy.map.Map('data/hinode/comp_XRT20220629_180339.8.fits')
    hinode_map_wo_light_leak = remove_lightleak(hinode_map)
    date_start = "2022-06-29 17:58"
    date_end = "2022-06-29 18:08"

# Query ARs and flares from HEK at the same time as the Hinode map
client = hek.HEKClient()
ars = client.search(a.Time(date_start, date_end), hek.attrs.EventType('AR'))
fls = client.search(a.Time(date_start, date_end), hek.attrs.EventType('FL'))

if obs == "flare":
    comp_map_regs_w_contours(hinode_map_wo_light_leak, ff_1=0.07, ff_2=3.37, ff_3=14.72, ars=ars, fls=fls)
elif obs == "qui":
    comp_map_regs_qui_w_contours(hinode_map_wo_light_leak, ff_1=21.31, ars=ars, fls=fls)