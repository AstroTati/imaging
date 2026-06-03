"""
Created on Fri May 12 13:26:56 2023
Tatiana M. Rodriguez

==================================================================
plot_radio_continuum.py
==================================================================
Plot a VLA radio continuum image with two datasets:
  - data_1: rendered as a color map (imshow)
  - data_2: overlaid as contours at multiples of the rms noise

Requirements: 
    2-axis FITS files (see drop2axis.py if needed).

Usage example:
    python plot_radio_continuum.py --file1 continuum_color.fits \
                                   --file2 continuum_contours.fits \
                                   [options]

Run 'python plot_radio_continuum.py --help' for all options.

Note:
    I created this script and used genAI to generalize, comment, 
    & optimize it. I tested it on my data but let me know if you
    encounter any issues. :)
    
"""

import argparse
import warnings
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib.patches import Ellipse

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_fits(path: str):
    """Return (data, WCS, primary HDU) for a FITS file."""
    hdu = fits.open(path)[0]
    return hdu.data, WCS(hdu.header), hdu


def cutout(data, wcs, hdu, center_pix: tuple, box_pix: tuple):
    """
    Trim data to a rectangular region.

    Parameters
    ----------
    center_pix : (x, y) pixel coordinates of the cutout centre.
    box_pix    : (nx, ny) cutout size in pixels.
    """
    size = u.Quantity(box_pix, u.pix)
    cut = Cutout2D(data, position=center_pix, size=size, wcs=wcs)
    hdu.header.update(cut.wcs.to_header())
    hdu.data = cut.data
    return cut.data, WCS(hdu.header), hdu


# ---------------------------------------------------------------------------
# Header parsers
# ---------------------------------------------------------------------------

def get_beam(hdr) -> tuple:
    """Return (major, minor, pa) beam parameters in degrees."""
    return hdr["BMAJ"], hdr["BMIN"], hdr["BPA"]


def get_pixel_scale(hdr) -> tuple:
    """Return (CDELT1, CDELT2) pixel scales in deg/pixel."""
    return hdr["CDELT1"], hdr["CDELT2"]


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def add_beam(ax, wcs_ref, beam_coord_deg: tuple, major: float, minor: float,
             pa: float, **ellipse_kwargs):
    """
    Overlay a beam ellipse on *ax*.

    Parameters
    ----------
    beam_coord_deg : (ra_deg, dec_deg) position of the ellipse centre.
    major, minor   : beam axes in degrees.
    pa             : beam position angle in degrees (BPA convention).
    ellipse_kwargs : forwarded to matplotlib Ellipse.
    """
    ra, dec = beam_coord_deg
    patch = Ellipse(
        (ra, dec), major, minor,
        angle=90.0 - pa,
        transform=ax.get_transform("fk5"),
        **ellipse_kwargs,
    )
    ax.add_patch(patch)


def add_scale_bar(ax, distance_kpc: float, bar_au: float, pix_scale_arcsec: float,
                  color: str = "w", fontsize: int = 22):
    """
    Draw a physical scale bar in the upper-right corner.

    Parameters
    ----------
    distance_kpc      : source distance in kpc.
    bar_au            : desired bar length in au.
    pix_scale_arcsec  : image pixel scale in arcsec/pixel.
    """
    bar_arcsec = bar_au / (distance_kpc * 1e3)          # angular size in arcsec
    bar_pix = bar_arcsec / pix_scale_arcsec

    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    xr, yr = x2 - 2, y2 - 1
    ax.plot([xr - bar_pix, xr], [yr, yr], lw=3, color=color)
    ax.text(
        (xr - bar_pix + xr) / 2.0, yr - 1.5,
        f"{bar_au:.0f} au",
        fontsize=fontsize, ha="center", color=color,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot a VLA radio continuum image (colour + contours)."
    )
    p.add_argument("--file1", required=True,
                   help="FITS file rendered as a colour image.")
    p.add_argument("--file2", required=True,
                   help="FITS file overlaid as contours.")

    # Cutout centres and sizes
    p.add_argument("--center1", nargs=2, type=int, default=[1048, 1224],
                   metavar=("X", "Y"), help="Pixel centre for file1 cutout.")
    p.add_argument("--box1", nargs=2, type=int, default=[55, 55],
                   metavar=("NX", "NY"), help="Cutout box size for file1 (pixels).")
    p.add_argument("--center2", nargs=2, type=int, default=[935, 1020],
                   metavar=("X", "Y"), help="Pixel centre for file2 cutout.")
    p.add_argument("--box2", nargs=2, type=int, default=[300, 300],
                   metavar=("NX", "NY"), help="Cutout box size for file2 (pixels).")

    # Plot window (pixels within the cutout)
    p.add_argument("--xlim", nargs=2, type=float, default=[7, 47],
                   metavar=("XMIN", "XMAX"))
    p.add_argument("--ylim", nargs=2, type=float, default=[5, 45],
                   metavar=("YMIN", "YMAX"))

    # Colour scale
    p.add_argument("--vmin", type=float, default=-1e-7,
                   help="imshow vmin (Jy/beam).")

    # Contour levels as multiples of rms
    p.add_argument("--sigma-levels", nargs="+", type=float,
                   default=[-3, 3, 5, 10, 15],
                   help="Contour levels in units of image rms.")

    # Scale bar
    p.add_argument("--distance", type=float, default=2.3,
                   help="Source distance in kpc.")
    p.add_argument("--bar-au", type=float, default=500,
                   help="Physical length of the scale bar in au.")
    p.add_argument("--pix-scale", type=float, default=0.05,
                   help="Pixel scale in arcsec/pixel.")

    # Beam position (FK5, sexagesimal)
    p.add_argument("--beam-ra", default="19h06m01.663s",
                   help="RA of beam ellipse centre.")
    p.add_argument("--beam-dec", default="6d46m35.35s",
                   help="Dec of beam ellipse centre.")

    # Labels
    p.add_argument("--title", default="Source name", help="Plot title.")
    p.add_argument("--cbar-unit", default=r"$\mu$Jy/beam",
                   help="Colour-bar axis label.")

    # Output
    p.add_argument("--output", default=None,
                   help="Save path (e.g. plot.pdf). Omit to display interactively.")

    return p.parse_args()


def main():
    args = parse_args()

    # ------------------------------------------------------------------
    # Load and trim data
    # ------------------------------------------------------------------
    data1, wcs1, hdu1 = load_fits(args.file1)
    data1 = np.ma.masked_invalid(data1)
    data1, wcs1, hdu1 = cutout(data1, wcs1, hdu1, args.center1, args.box1)

    data2, wcs2, hdu2 = load_fits(args.file2)
    data2 = np.ma.masked_invalid(data2)
    rms = float(data2.std())
    data2, wcs2, hdu2 = cutout(data2, wcs2, hdu2, args.center2, args.box2)

    # ------------------------------------------------------------------
    # Figure setup  — WCS of data1 (colour image) drives the projection
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs1)

    ra_ax = ax.coords[0]
    dec_ax = ax.coords[1]
    ra_ax.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
    dec_ax.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
    ra_ax.set_major_formatter("hh:mm:ss.ss")
    ra_ax.display_minor_ticks(True)
    dec_ax.display_minor_ticks(True)
    ax.tick_params(which="both", direction="in", color="white",
                   length=10, width=2, labelsize=20)
    ax.tick_params(which="minor", length=5)
    ax.set_title(args.title, fontsize=25)
    ax.set_xlim(*args.xlim)
    ax.set_ylim(*args.ylim)

    # ------------------------------------------------------------------
    # Colour image
    # ------------------------------------------------------------------
    im = ax.imshow(
        data1, origin="lower", cmap="plasma",
        vmin=args.vmin, vmax=float(data1.max()),
        transform=ax.get_transform(wcs1),
    )

    cbar = fig.colorbar(im, shrink=0.8)
    cbar.ax.set_ylabel(args.cbar_unit, fontsize=20)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(18)

    # ------------------------------------------------------------------
    # Contours
    # ------------------------------------------------------------------
    levels = [rms * s for s in args.sigma_levels]
    plt.rcParams["lines.linewidth"] = 2
    ax.contour(
        data2, levels=levels, colors="k",
        transform=ax.get_transform(wcs2),
    )

    # ------------------------------------------------------------------
    # Scale bar
    # ------------------------------------------------------------------
    add_scale_bar(
        ax,
        distance_kpc=args.distance,
        bar_au=args.bar_au,
        pix_scale_arcsec=args.pix_scale,
    )

    # ------------------------------------------------------------------
    # Beam ellipses
    # ------------------------------------------------------------------
    beam_ra = Angle(args.beam_ra).deg
    beam_dec = Angle(args.beam_dec).deg
    beam_coord = (beam_ra, beam_dec)

    maj1, min1, pa1 = get_beam(hdu1.header)
    add_beam(ax, wcs1, beam_coord, maj1, min1, pa1,
             edgecolor="w", facecolor="none", lw=3)

    maj2, min2, pa2 = get_beam(hdu2.header)
    plt.rcParams["hatch.linewidth"] = 3
    add_beam(ax, wcs1, beam_coord, maj2, min2, pa2,
             edgecolor="w", facecolor="w")

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    if args.output:
        plt.savefig(args.output, bbox_inches="tight")
        print(f"Saved → {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
