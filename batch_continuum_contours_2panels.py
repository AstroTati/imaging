"""
Created on Thu Dec  8 15:34:53 2022
Tatiana M. Rodriguez

=================================================================================
plot_continuum_batch.py
=================================================================================
Batch-generate two-panel continuum plots for a list of VLA sources.

Panel 1 (wide field): two continuum fits files plotted in contours, with
                      a zoom rectangle indicating the panel-2 region.
                      
Panel 2 (zoom):       contours of one data set only, with source labels and masers.

Both panels share the same beam, scale-bar, maser, and label logic via
helper functions to avoid code duplication.

Required inputs:
    - Per-source FITS files:  <fits_dir>/<src>_data1.fits
                              <fits_dir>/<src>_data2.fits
    - continuum_info.csv      tab-separated catalogue (see column list below)
    - CIIMMs/all.csv          catalogue of CH3OH/H2O masers

continuum_info.csv columns:
    Source, center, box, d, data1_cont, data2_cont, labels, lab_RA, lab_Dec, WM_if, 
    WM_RA, WM_Dec

Usage example:
    python plot_continuum_batch.py [options]

    python plot_continuum_batch.py --fits-dir /data/fits \
                                   --catalogue continuum_info.csv \
                                   --masers CIIMMs/all.csv \
                                   --output-dir plots/ \
                                   --sources 18345 20293

Note: 
    I created this script and used genAI to generalize, comment and 
    optimize it. I tested it on my own data but let me know if you 
    encounter any issues. Good luck! :)
    
"""

import argparse
import ast
import warnings
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib.patches import Ellipse, Rectangle

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_fits(path: str):
    """Return (data, WCS, primary HDU). WCS is forced to 2 axes."""
    hdu = fits.open(path)[0]
    return hdu.data, WCS(hdu.header, naxis=2), hdu


def read_maser_catalogue(path: str):
    """Return (ra_array, dec_array) from a tab-separated maser catalogue."""
    df = pd.read_csv(path, sep="\t")
    ra  = df["RA"].apply(ast.literal_eval).explode().astype(float).to_numpy()
    dec = df["DEC"].apply(ast.literal_eval).explode().astype(float).to_numpy()
    return ra, dec


def safe_eval(value):
    """eval() wrapper that also handles plain scalars stored as strings."""
    if isinstance(value, str):
        return ast.literal_eval(value)
    return value


# ---------------------------------------------------------------------------
# FITS helpers
# ---------------------------------------------------------------------------

def apply_cutout(data, wcs, hdu, center_pix: tuple, box_pix: int):
    """
    Trim *data* to a square cutout of side *box_pix* centred on *center_pix*.
    Returns updated (data, wcs, hdu).
    """
    size = u.Quantity((box_pix, box_pix), u.pix)
    cut = Cutout2D(data, position=center_pix, size=size, wcs=wcs)
    hdu.header.update(cut.wcs.to_header())
    hdu.data = cut.data
    return cut.data, WCS(hdu.header), hdu


def get_beam(hdr) -> tuple:
    """Return (major_deg, minor_deg, pa_deg) from a FITS header."""
    return hdr["BMAJ"], hdr["BMIN"], hdr["BPA"]


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def setup_axes(ax, wcs, title: str, fontsize: int = 24):
    """Apply standard WCS axis formatting to *ax*."""
    ra_ax  = ax.coords[0]
    dec_ax = ax.coords[1]
    ra_ax.set_axislabel("RA (J2000)",  minpad=0.8,  fontsize=fontsize + 1)
    dec_ax.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=fontsize + 1)
    ra_ax.set_major_formatter("hh:mm:ss.ss")
    ra_ax.display_minor_ticks(True)
    dec_ax.display_minor_ticks(True)
    ax.tick_params(which="both", direction="in", color="k",
                   length=15, width=2, labelsize=fontsize)
    ax.tick_params(which="minor", length=5)
    ax.set_title(title, fontsize=fontsize)


def add_beam_ellipse(ax, wcs, beam_pos_deg: tuple,
                     major: float, minor: float, pa: float,
                     **ellipse_kwargs):
    """Overlay a beam ellipse at *beam_pos_deg* = (ra_deg, dec_deg)."""
    patch = Ellipse(
        beam_pos_deg, major, minor,
        angle=90.0 - pa,
        transform=ax.get_transform("icrs"),
        **ellipse_kwargs,
    )
    ax.add_patch(patch)


def add_scale_bar(ax, wcs, box_pix: int,
                  distance_kpc: float, bar_au: float,
                  pix_scale_arcsec: float,
                  color: str = "k", fontsize: int = 24):
    """
    Draw a physical scale bar anchored to the upper-right of the cutout.

    Parameters
    ----------
    box_pix          : cutout side length in pixels (sets anchor position).
    distance_kpc     : source distance in kpc.
    bar_au           : desired bar length in au.
    pix_scale_arcsec : pixel scale in arcsec/pixel.
    """
    bar_pix = bar_au / (distance_kpc * 1e3) / pix_scale_arcsec  # arcsec → pixels
    margin_x = box_pix / 10
    margin_y = box_pix / 25

    xr = box_pix - margin_x
    yr = box_pix - margin_y
    ax.plot([xr - bar_pix, xr], [yr, yr], color=color, lw=8)
    ax.text(
        xr - bar_pix / 2, yr - margin_y,
        f"{bar_au:.0f} au",
        color=color, ha="center", va="center", fontsize=fontsize,
    )


def plot_masers(ax, maser_ra, maser_dec, wcs, center_pix, box_pix: int,
                marker: str = "r+", markersize: int = 30, mew: int = 3):
    """
    Plot masers that fall within ~20 arcsec of the cutout centre.

    Uses a simple RA offset filter to avoid projecting every catalogue entry.
    """
    centre_sky = wcs.pixel_to_world(center_pix[0], center_pix[1])
    for ra, dec in zip(maser_ra, maser_dec):
        pos = SkyCoord(ra, dec, frame="fk5", unit="deg")
        delta_ra = (pos.ra - centre_sky.ra).to(u.arcsec)
        if abs(delta_ra.value) <= 20:
            ax.plot(ra, dec, marker, fillstyle="full",
                    markersize=markersize, mew=mew, alpha=1,
                    markerfacecolor="r",
                    transform=ax.get_transform("fk5"))


def plot_labels(ax, wcs, lab_ra, lab_dec, labels,
                xlim=None, ylim=None, fontsize: int = 35):
    """
    Annotate continuum components.

    If *xlim*/*ylim* are provided (pixel ranges), only labels inside the
    current view are drawn — prevents crowding on the zoom panel.
    """
    if xlim is not None and ylim is not None:
        # Convert pixel limits to sky coordinates for the boundary check
        blc = wcs.pixel_to_world(xlim[0], ylim[0])
        trc = wcs.pixel_to_world(xlim[1], ylim[1])
        ra_min  = min(blc.ra.value,  trc.ra.value)
        ra_max  = max(blc.ra.value,  trc.ra.value)
        dec_min = min(blc.dec.value, trc.dec.value)
        dec_max = max(blc.dec.value, trc.dec.value)
        visible = lambda ra, dec: ra_min <= ra <= ra_max and dec_min <= dec <= dec_max
    else:
        visible = lambda ra, dec: True

    for ra, dec, lbl in zip(lab_ra, lab_dec, labels):
        if visible(ra, dec):
            ax.text(ra, dec, lbl, size=fontsize, color="k",
                    transform=ax.get_transform("fk5"))


def format_contour_annotation(levels_jy: list, rms_jy: float) -> str:
    """
    Return a LaTeX-ready string like '[−3.0, 3.0, 5.0] × 12.34 μJy beam⁻¹'.
    """
    sigma_levels = [f"{v / rms_jy:.1f}" for v in levels_jy]
    rms_uJy = f"{rms_jy * 1e6:.2f}"
    return (
        "[" + ", ".join(sigma_levels) + "]"
        + rf" $\times$ {rms_uJy} $\mu$Jy$\,$beam$^{{-1}}$"
    )


# ---------------------------------------------------------------------------
# Per-source plotting
# ---------------------------------------------------------------------------

def plot_source(row, fits_dir: Path, output_dir: Path,
                maser_ra, maser_dec,
                pix_scale: float = 0.013):
    """
    Generate and save a two-panel figure for one source.

    Parameters
    ----------
    row       : pandas Series from continuum_info.csv.
    fits_dir  : directory containing the FITS files.
    output_dir: directory where the PDF will be saved.
    maser_ra/dec : full maser catalogue arrays.
    pix_scale : pixel scale in arcsec/pixel (same for all data).
    """
    src  = row["Source"]
    d_kpc = float(row["d"])

    # Parse catalogue columns that store Python literals as strings
    center          = safe_eval(row["center"])          # [ra_deg, dec_deg]
    box             = safe_eval(row["box"])             # [wide_pix, zoom_pix]
    data1_levels    = safe_eval(row["data1_cont"])
    data2_levels    = safe_eval(row["data2_cont"])
    labels          = safe_eval(row["labels"])
    lab_ra          = safe_eval(row["lab_RA"])
    lab_dec         = safe_eval(row["lab_Dec"])

    # ------------------------------------------------------------------
    # Load FITS data
    # ------------------------------------------------------------------
    data1_fits  = str(fits_dir / f"{src}_data1.fits")
    data2_fits = str(fits_dir / f"{src}_data2.fits")

    data1_data,  data1_wcs,  data1_hdu  = load_fits(data1_fits)
    data1_data = np.ma.masked_invalid(data1_data)
    rms = float(data1_data.std())

    data2_data, data2_wcs, data2_hdu = load_fits(data2_fits)
    data2_data = np.ma.masked_invalid(data2_data)
    # Source-specific rms override (e.g. known artefact in 18345)
    rmsK = 3e-5 if src == "18345" else float(data2_data.std())

    # Common cutout centre in sky coords
    center_sky = SkyCoord(center[0], center[1], frame="fk5", unit="deg")

    # ------------------------------------------------------------------
    # Panel 1: wide-field cutout (box[0] pixels)
    # ------------------------------------------------------------------
    wide_box = box[0]
    d1, w1, h1 = load_fits(data1_fits)   # fresh copy — cutout modifies hdu in place
    d1 = np.ma.masked_invalid(d1)
    c1 = w1.world_to_pixel(center_sky)
    d1, w1, h1 = apply_cutout(d1, w1, h1, c1, wide_box)

    # ------------------------------------------------------------------
    # Panel 2: zoom cutout (box[1] pixels)
    # ------------------------------------------------------------------
    zoom_box = box[1]
    d2, w2, h2 = load_fits(data1_fits)
    d2 = np.ma.masked_invalid(d2)
    c2 = w2.world_to_pixel(center_sky)
    d2, w2, h2 = apply_cutout(d2, w2, h2, c2, zoom_box)

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(36, 15))
    plt.subplots_adjust(hspace=0.2)

    # --- Panel 1 ---
    ax1 = fig.add_subplot(1, 2, 1, projection=w1)
    setup_axes(ax1, w1, src)
    ax1.set_xlim(-0.5, wide_box - 0.5)
    ax1.set_ylim(-0.5, wide_box - 0.5)

    # Data1 contours (black)
    ax1.contour(d1, levels=data1_levels,  colors="k",
                transform=ax1.get_transform(w1))
    # Data2 contours (green)
    ax1.contour(data2_data, levels=data2_levels, colors="g",
                transform=ax1.get_transform(data2_wcs))

    # Beams
    beam_pos1 = w1.pixel_to_world(wide_box / 12, wide_box / 12)
    bp1 = (beam_pos1.ra.deg, beam_pos1.dec.deg)
    add_beam_ellipse(ax1, w1, bp1, *get_beam(data2_hdu.header),
                     edgecolor="g", facecolor="none", lw=3)
    add_beam_ellipse(ax1, w1, bp1, *get_beam(h1.header),
                     edgecolor="k", facecolor="k", lw=3)

    add_scale_bar(ax1, w1, wide_box, d_kpc, bar_au=1000,
                  pix_scale_arcsec=pix_scale)

    plot_labels(ax1, w1, lab_ra, lab_dec, labels)

    if row["WM_if"]:
        wm_ra  = safe_eval(row["WM_RA"])
        wm_dec = safe_eval(row["WM_Dec"])
        plot_masers(ax1, wm_ra, wm_dec, w1,
                    (wide_box / 2, wide_box / 2), wide_box)

    plot_masers(ax1, maser_ra, maser_dec, w1,
                (wide_box / 2, wide_box / 2), wide_box,
                marker="ko", markersize=10, mew=1)

    # Zoom-region rectangle
    offset = (wide_box - zoom_box) / 2
    ax1.add_patch(Rectangle(
        (offset, offset), zoom_box, zoom_box,
        edgecolor="k", facecolor="none", lw=1.5,
    ))

    # Contour-level annotation below the panel
    x_lim = ax1.get_xlim()
    ax1.text(
        min(x_lim), -(wide_box / 6),
        format_contour_annotation(data2_levels, rmsK),
        ha="left", fontsize=22, color="g",
    )

    # --- Panel 2 ---
    ax2 = fig.add_subplot(1, 2, 2, projection=w2)
    setup_axes(ax2, w2, src)
    ax2.set_xlim(-0.5, zoom_box - 0.5)
    ax2.set_ylim(-0.5, zoom_box - 0.5)

    ax2.contour(d2, levels=data1_levels, colors="k",
                transform=ax2.get_transform(w2))

    beam_pos2 = w2.pixel_to_world(zoom_box / 14, zoom_box / 15)
    bp2 = (beam_pos2.ra.deg, beam_pos2.dec.deg)
    add_beam_ellipse(ax2, w2, bp2, *get_beam(h2.header),
                     edgecolor="k", facecolor="k", lw=3)

    add_scale_bar(ax2, w2, zoom_box, d_kpc, bar_au=500,
                  pix_scale_arcsec=pix_scale)

    plot_labels(ax2, w2, lab_ra, lab_dec, labels,
                xlim=ax2.get_xlim(), ylim=ax2.get_ylim())

    if row["WM_if"]:
        plot_masers(ax2, wm_ra, wm_dec, w2,
                    (zoom_box / 2, zoom_box / 2), zoom_box,
                    markersize=30, mew=4)

    plot_masers(ax2, maser_ra, maser_dec, w2,
                (zoom_box / 2, zoom_box / 2), zoom_box,
                marker="ko", markersize=12, mew=1)

    x_lim2 = ax2.get_xlim()
    ax2.text(
        min(x_lim2), -(zoom_box / 6),
        format_contour_annotation(data1_levels, rms),
        ha="left", fontsize=22,
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    out_path = output_dir / f"{src}_continuum.pdf"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved → {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Batch-generate two-panel VLA continuum plots."
    )
    p.add_argument("--fits-dir", default=".",
                   help="Directory containing the per-source FITS files.")
    p.add_argument("--catalogue", default="continuum_info.csv",
                   help="Tab-separated source catalogue (continuum_info.csv).")
    p.add_argument("--masers", default="CIIMMs/all.csv",
                   help="Tab-separated maser catalogue.")
    p.add_argument("--output-dir", default=".",
                   help="Directory for output PDFs.")
    p.add_argument("--sources", nargs="*", default=None,
                   help="Process only these source names. Omit to run all.")
    p.add_argument("--pix-scale", type=float, default=0.013,
                   help="Pixel scale in arcsec/pixel (default: 0.013).")
    return p.parse_args()


def main():
    args = parse_args()

    fits_dir   = Path(args.fits_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.catalogue, sep="\t")
    maser_ra, maser_dec = read_maser_catalogue(args.masers)

    # Filter to requested sources (preserving catalogue order)
    if args.sources:
        requested = set(args.sources)
        df = df[df["Source"].isin(requested)]
        if df.empty:
            raise SystemExit(f"None of {args.sources} found in {args.catalogue}.")

    for _, row in df.iterrows():
        try:
            plot_source(row, fits_dir, output_dir, maser_ra, maser_dec,
                        pix_scale=args.pix_scale)
        except Exception as exc:
            print(f"[WARN] Skipping {row['Source']}: {exc}")


if __name__ == "__main__":
    main()
