"""
Created on Tue Dec 12 11:10:49 2023
Tatiana M. Rodriguez

============================================================================
fits_cutout.py
============================================================================
Extract a square (or rectangular) cutout from a FITS file and write it to
a new file with a corrected WCS header.

Usage example:
    python fits_cutout.py --input data.fits --output cutout.fits \
                          --center 1218 1284

    # Custom box size
    python fits_cutout.py --input data.fits --output cutout.fits \
                          --center 1218 1284 --size 500

    # Non-square box
    python fits_cutout.py --input data.fits --output cutout.fits \
                          --center 1218 1284 --size 600 400

Note:
    I used copilot to optimize this script I wrote a few years back.
    I tested it on my data, but let me know if you encounter any issues.

"""

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def fits_cutout(input_path: str, output_path: str,
                center_x: int, center_y: int,
                size_x: int = 350, size_y: int = None) -> None:
    """
    Write a cutout of a FITS image to *output_path*.

    The cutout is clipped to the image boundary if the requested box
    extends beyond the edge, and the WCS reference pixel is updated
    accordingly so the astrometry stays correct.

    Parameters
    ----------
    input_path  : Path to the source FITS file.
    output_path : Path for the output cutout FITS file.
    center_x    : Cutout centre in pixel X (NAXIS1 direction).
    center_y    : Cutout centre in pixel Y (NAXIS2 direction).
    size_x      : Width of the cutout in pixels (default 350).
    size_y      : Height of the cutout in pixels. Defaults to size_x (square).
    """
    if size_y is None:
        size_y = size_x

    with fits.open(input_path) as hdul:
        data   = hdul[0].data
        header = hdul[0].header.copy()

    if data.ndim not in (2, 3, 4):
        raise ValueError(f"Unexpected data shape {data.shape}. "
                         "Expected 2-, 3-, or 4-dimensional FITS data.")

    # Work on the last two axes (handles cubes well)
    ny = data.shape[-2]
    nx = data.shape[-1]

    # Clip the cutout window to image boundaries
    x0 = max(center_x - size_x // 2, 0)
    x1 = min(x0 + size_x, nx)
    y0 = max(center_y - size_y // 2, 0)
    y1 = min(y0 + size_y, ny)

    cutout = data[..., y0:y1, x0:x1]

    # Update the WCS reference pixel to match the new array origin
    header["NAXIS1"] = x1 - x0
    header["NAXIS2"] = y1 - y0
    header["CRPIX1"] = header.get("CRPIX1", 0) - x0
    header["CRPIX2"] = header.get("CRPIX2", 0) - y0

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fits.writeto(output_path, cutout, header, overwrite=True)
    print(f"Cutout [{x0}:{x1}, {y0}:{y1}] → {output_path}")


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract a cutout from a FITS file with corrected WCS."
    )
    p.add_argument("--input",  required=True, help="Input FITS file path.")
    p.add_argument("--output", required=True, help="Output FITS file path.")
    p.add_argument("--center", nargs=2, type=int, required=True,
                   metavar=("X", "Y"), help="Cutout centre in pixels.")
    p.add_argument("--size", nargs="+", type=int, default=[350],
                   metavar="N",
                   help="Box size in pixels. One value → square. "
                        "Two values → width height.")
    return p.parse_args()


def main():
    args = parse_args()

    if len(args.size) == 1:
        size_x = size_y = args.size[0]
    elif len(args.size) == 2:
        size_x, size_y = args.size
    else:
        raise SystemExit("--size accepts one or two values.")

    fits_cutout(
        input_path=args.input,
        output_path=args.output,
        center_x=args.center[0],
        center_y=args.center[1],
        size_x=size_x,
        size_y=size_y,
    )


if __name__ == "__main__":
    main()