
import numpy as np
import scipy.ndimage as nd

import astropy.units as u
from astropy.convolution import Gaussian1DKernel

from spectral_cube import SpectralCube
from radio_beam import Beam

from regions import Regions
import seaborn as sns

import matplotlib.pyplot as plt


from uvcombine.scale_factor import find_scale_factor

from cube_analysis.feather_cubes import feather_compare_cube
from cube_analysis.register_cubes import cube_registration


from pathlib import Path

sd_data_path = Path("/reduction10/erickoch/LGLBS/hi_feathering/")
vla_data_path = Path("/reduction10/erickoch/LGLBS/C+D_HI_2023/")

galaxy_dict = { 'ngc6822': ["NGC6822-center_cube.fits",
                           "ngc6822_C+D_hilores.fits",
                           [-75, -25] * u.km / u.s],
               'wlm': ['WLM_GBT.FITS',
                       'wlm_C+D_hilores.fits',
                       [-170, -70] * u.km / u.s],
               'ic10': ['IC10_GBT_Jy.fits',
                        'ic10_C+D_hilores.fits',
                        [-380, -280] * u.km / u.s],
                'ic1613': ['IC1613_GBT_vegas_K_noresample_lsrk.fits',
                        'ic1613_C+D_hilores.fits',
                        [-265, -210] * u.km / u.s],
               }

def taper_weights(pb_plane,
                erosion_interations=5):

    # mask = np.logical_or(np.isfinite(pb_plane), pb_plane > 0.)
    mask = pb_plane > 0.

    if erosion_interations > 0:
        mask = nd.binary_erosion(mask, iterations=erosion_interations)

    smoothed_weights = nd.gaussian_filter(mask.astype(float), sigma=erosion_interations/2.)

    weight_arr = pb_plane * smoothed_weights

    return weight_arr


for this_gal in galaxy_dict:

    this_gbt_filename = galaxy_dict[this_gal][0]
    this_vla_filename = galaxy_dict[this_gal][1]

    this_specslice_low, this_specslice_high = galaxy_dict[this_gal][2]

    vla_cube = SpectralCube.read(vla_data_path / this_vla_filename)
    vla_cube.allow_huge_operations = True

    # Nick's reprocessed and gridded GBT cube.
    gbt_cube = SpectralCube.read(sd_data_path / this_gbt_filename)
    gbt_cube.allow_huge_operations = True

    # Use the proper beam model size, not the one in the header!
    gbt_beam_model = Beam(area=3.69e5 *u.arcsec**2)
    gbt_beam_model.major.to(u.arcmin)

    gbt_cube = gbt_cube.with_beam(gbt_beam_model, raise_error_jybm=False)
    gbt_cube = gbt_cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

    fwhm_factor = np.sqrt(8*np.log(2))
    current_resolution = np.abs(np.diff(gbt_cube.spectral_axis)[0]).to(u.km / u.s)
    target_resolution = np.abs(np.diff(vla_cube.spectral_axis)[0]).to(u.km / u.s)

    if current_resolution < target_resolution:
        gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                        current_resolution / fwhm_factor)
        print(gaussian_width)
        kernel = Gaussian1DKernel(gaussian_width.value)
        gbt_cube_specsmooth = gbt_cube.spectral_smooth(kernel)
    else:
        gbt_cube_specsmooth = gbt_cube


    gbt_cube_specinterp = gbt_cube_specsmooth.spectral_interpolate(vla_cube.spectral_axis)

    # gbt_cube_specinterp.mean(axis=(1, 2)).quicklook()
    # vla_cube.mean(axis=(1, 2)).quicklook()

    # To speed up the computation, select a velocity range corresponding to the target of interest.

    chan_min = vla_cube.closest_spectral_channel(this_specslice_low)
    chan_max = vla_cube.closest_spectral_channel(this_specslice_high)

    if chan_min > chan_max:
        chan_min, chan_max = chan_max, chan_min

    print(f"Channels to slice: {chan_min}, {chan_max}")

    spec_slicer = slice(chan_min, chan_max)

    # This function will smooth to the lowsest resolution and reproject per channel
    # https://github.com/e-koch/CubeAnalysis/blob/master/cube_analysis/register_cubes.py#L58

    spatial_offsets = cube_registration(gbt_cube_specinterp[spec_slicer],
                                        vla_cube[spec_slicer],
                                        verbose=True, num_cores=1,
                                        restfreq=1.42040575177 * u.GHz,
                                        check_specaxis=True,)


    cmap = sns.color_palette("icefire", as_cmap=True)

    fig, ax = plt.subplots()
    points = ax.scatter(spatial_offsets[:, 1],
                        spatial_offsets[:, 0],
                        c=vla_cube[spec_slicer].spectral_axis.to(u.km / u.s).value,
                        s=50,
                        cmap=cmap)
    fig.colorbar(points, label="Velocity (km/s)")

    ax.set_ylabel("Dec pixel offset")
    ax.set_xlabel("RA pixel offset")

    fig.savefig(sd_data_path / f"{this_gal}_spatialoffsets.pdf")
    plt.close(fig)

    if (sd_data_path / f"{this_gal}_spatialoffsets.npy").exists():
        (sd_data_path / f"{this_gal}_spatialoffsets.npy").unlink()

    np.save(sd_data_path / f"{this_gal}_spatialoffsets.npy",
            spatial_offsets)


    vla_reproj_plane = vla_cube[0].reproject(gbt_cube_specinterp[0].header)
    gbt_cube_specinterp_masked = gbt_cube_specinterp.with_mask(np.isfinite(vla_reproj_plane)).minimal_subcube(spatial_only=True)

    # gbt_cube_specinterp_masked[spec_slicer].moment0().quicklook()

    hi_freq = 1.42040575 * u.GHz
    las = (hi_freq.to(u.cm, u.spectral()) / (40 * u.m)).to(u.arcsec, u.dimensionless_angles())

    weights = taper_weights(np.isfinite(vla_cube[0]), erosion_interations=20)

    radii, ratios, high_pts, low_pts, chan_out = \
                feather_compare_cube(vla_cube.to(u.K),
                                    gbt_cube_specinterp_masked.to(u.K),
                                    las,
                                    lowresfwhm=None,
                                    num_cores=1,
                                    chunk=250,
                                    verbose=True,
                                    weights=weights,
                                    relax_spectral_check=False,
                                    spec_check_kwargs={'rtol': 0.01})

    if (sd_data_path / f"{this_gal}_overlap_vla_samples.npy").exists():
        (sd_data_path / f"{this_gal}_overlap_vla_samples.npy").unlink()

    np.save(sd_data_path / f"{this_gal}_overlap_vla_samples.npy",
<<<<<<< HEAD
            spatial_offsets)
=======
            high_pts)
>>>>>>> 24ab4b9 (Long overdue update, including foreground separation, feathering scripts, rms estimation, and cluster imaging)

    if (sd_data_path / f"{this_gal}_overlap_gbt_samples.npy").exists():
        (sd_data_path / f"{this_gal}_overlap_gbt_samples.npy").unlink()

    np.save(sd_data_path / f"{this_gal}_overlap_gbt_samples.npy",
<<<<<<< HEAD
            spatial_offsets)
=======
            low_pts)
>>>>>>> 24ab4b9 (Long overdue update, including foreground separation, feathering scripts, rms estimation, and cluster imaging)

    sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[:-1]),
                                        np.hstack(high_pts[:-1]),
                                        method='distrib',
                                        verbose=True)
    print(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")
    plt.title(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")

    plt.savefig(sd_data_path / f"{this_gal}_sdfactor_allchans.pdf")
    plt.close()


    sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[spec_slicer]),
                                        np.hstack(high_pts[spec_slicer]),
                                        method='distrib',
                                        verbose=True)
    print(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
    plt.title(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
    plt.xlim([-1, 1])

    plt.savefig(sd_data_path / f"{this_gal}_sdfactor_galchans.pdf")
    plt.close()


    # Check for a dependence on velocity:

    fig, axs = plt.subplots(2, 1, sharex=False)

    cmap = sns.color_palette("icefire", n_colors=chan_max - chan_min, as_cmap=False)


    for ii, this_ratio in enumerate(ratios[spec_slicer]):

        _ = axs[0].hist(this_ratio, label=f"{ii}", color=cmap[ii])

        axs[1].plot(2 * [np.median(this_ratio)], [0, 1], color=cmap[ii])

    fig.savefig(sd_data_path / f"{this_gal}_sdfactor_perchannel.pdf")
    plt.close(fig)

