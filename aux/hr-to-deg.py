"""
Created on Thu Mar 17 13:05:50 2022
Tatiana M. Rodriguez

Going between degrees and hour angles is very annoying but also a very 
common occurrance. I made this short script to make quick transformations.
It looks more complicated than it actually is ;)

Usage example:
    # Sexagesimal -> degrees (default)
    python convert_coords.py --ra "18h21m08.9907s" --dec "-14d31m46.685s"
 
    # Degrees -> sexagesimal
    python convert_coords.py --ra 275.287 --dec -14.529 --to-sex
 
    # From a tab-separated CSV (columns: RA, Dec)
    python convert_coords.py --csv coords.csv
    python convert_coords.py --csv coords.csv --to-sexag
    
"""
 
import argparse
import numpy as np
import pandas as pd
from astropy.coordinates import Angle
import astropy.units as u
 
 
def to_degrees(ra_str: str, dec_str: str) -> tuple:
    """Convert sexagesimal RA/Dec strings to decimal degrees."""
    return Angle(ra_str).degree, Angle(dec_str).degree
 
 
def to_sexagesimal(ra_deg: float, dec_deg: float) -> tuple:
    """Convert decimal degrees to sexagesimal strings (hms / dms)."""
    ra_sex  = Angle(ra_deg,  unit=u.deg).to_string(unit=u.hour,       sep="hms", precision=7)
    dec_sex = Angle(dec_deg, unit=u.deg).to_string(unit=u.degree,     sep="dms", precision=7)
    return ra_sex, dec_sex
 
 
def convert_list(ra_values, dec_values, to_sexag: bool) -> tuple:
    """
    Convert a list of RA and Dec values in either direction.
 
    Parameters
    ----------
    ra_values, dec_values : sequences of strings (sexagesimal) or floats (degrees).
    to_sexag              : if True, degrees -> sexagesimal; otherwise the reverse.
    """
    ra_out, dec_out = [], []
    for ra, dec in zip(ra_values, dec_values):
        if to_sexag:
            r, d = to_sexagesimal(float(ra), float(dec))
        else:
            r, d = to_degrees(str(ra), str(dec))
        ra_out.append(r)
        dec_out.append(d)
    return ra_out, dec_out
 
 
def print_results(ra_out, dec_out):
    print("RA:  ", ", ".join(str(v) for v in ra_out))
    print("Dec: ", ", ".join(str(v) for v in dec_out))
 
 
def parse_args():
    p = argparse.ArgumentParser(
        description="Convert RA/Dec between sexagesimal and decimal degrees."
    )
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--ra",  nargs="+", help="RA value(s).")
    source.add_argument("--csv", help="Tab-separated CSV with RA and Dec columns.")
 
    p.add_argument("--dec", nargs="+",
                   help="Dec value(s): required when using --ra.")
    p.add_argument("--to-sexag", action="store_true",
                   help="Convert degrees -> sexagesimal (default: sexagesimal -> degrees).")
    return p.parse_args()
 
 
def main():
    args = parse_args()
 
    if args.csv:
        df = pd.read_csv(args.csv, sep="\t")
        ra_in  = df["RA"].tolist()
        dec_in = df["Dec"].tolist()
    else:
        if not args.dec:
            raise SystemExit("--dec is required when using --ra.")
        ra_in  = args.ra
        dec_in = args.dec
 
    ra_out, dec_out = convert_list(ra_in, dec_in, to_sexag=args.to_sexag)
    print_results(ra_out, dec_out)
 
 
if __name__ == "__main__":
    main()



