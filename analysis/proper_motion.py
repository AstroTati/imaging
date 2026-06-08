"""
Created on November, 2025
Author: Tatiana M. Rodriguez

================================================================================
proper_motion.py
================================================================================
Compute angular separation, proper motion, projected and deprojected velocity
between two epochs of VLA astrometry, with full error propagation.

Error budget per position:
        sigma_total = sqrt(measurement_error^2 + phase_calibrator_error^2)
RA errors are multiplied by cos(dec) before combining with Dec errors.

Usage example:
    # Run the built-in source catalogue
    python proper_motion.py

    # Custom inclination or output file
    python proper_motion.py --inclination 45 --output results.csv

    # Single source from the command line
    python proper_motion.py --ra1 "07h05m10.94082s" --dec1 "-12d19m00.47433s" \
                            --date1 24/11/2013 \
                            --ra2 "07h05m10.93861s" --dec2 "-12d19m00.4672s" \
                            --date2 14/05/2022 \
                            --distance 1000 --pos-acc 2
"""

import argparse
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PMResult:
    source:        str
    ang_sep_mas:   float
    ang_sep_err:   float
    lin_sep_au:    float
    lin_sep_err:   float
    v_proj_kms:    float
    v_proj_err:    float
    v_deproj_kms:  float
    v_deproj_err:  float
    pm_mas_yr:     float
    pm_err_mas_yr: float
    date1:         str
    date2:         str


# ---------------------------------------------------------------------------
# Core calculation
# ---------------------------------------------------------------------------

def proper_motion_calc(
    ra1: str, dec1: str, date1: str,
    ra2: str, dec2: str, date2: str,
    distance_pc: float, pos_acc_mas: float,
    inclination_deg: float = 57.0,
    err_ra1: float = 0.0, err_dec1: float = 0.0,
    err_ra2: float = 0.0, err_dec2: float = 0.0,
) -> dict:
    """
    Compute proper motion and velocity between two astrometric epochs.

    Parameters:
        ra1, dec1       : Epoch-1 position strings (hms / dms).
        date1           : Epoch-1 date string 'dd/mm/yyyy'.
        ra2, dec2       : Epoch-2 position strings.
        date2           : Epoch-2 date string 'dd/mm/yyyy'.
        distance_pc     : Source distance in parsecs.
        pos_acc_mas     : Phase-calibrator positional accuracy in mas.
        inclination_deg : Disk/outflow inclination in degrees.
        err_ra1/2       : RA  measurement errors in arcsec (default 0).
        err_dec1/2      : Dec measurement errors in arcsec (default 0).

    Returns:
        dict with keys: ang_sep_mas, ang_sep_err_mas, lin_sep_au, lin_sep_err_au,
                        v_proj_kms, v_proj_err_kms, v_deproj_kms, v_deproj_err_kms,
                        pm_mas_yr, pm_err_mas_yr
    """
    t1 = datetime.strptime(date1, "%d/%m/%Y")
    t2 = datetime.strptime(date2, "%d/%m/%Y")
    delta_t = abs((t2 - t1).total_seconds()) * u.s
    if delta_t.value == 0:
        raise ValueError("date1 and date2 must differ.")

    coord1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg), distance=distance_pc * u.pc)
    coord2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg), distance=distance_pc * u.pc)

    ang_sep  = coord1.separation(coord2)
    lin_sep  = coord1.separation_3d(coord2).to(u.AU)

    # --- Error propagation ---
    pos_acc_arcsec = (pos_acc_mas * u.mas).to(u.arcsec).value

    # Quadrature sum of measurement error and phase-calibrator error per axis
    def total_err(meas_err):
        return np.sqrt(meas_err**2 + pos_acc_arcsec**2)  # arcsec

    sig_ra1  = total_err(err_ra1)
    sig_dec1 = total_err(err_dec1)
    sig_ra2  = total_err(err_ra2)
    sig_dec2 = total_err(err_dec2)

    # RA errors scaled to true angle on sky
    cos_dec = np.cos((coord1.dec.radian + coord2.dec.radian) / 2)
    sig_ra1 *= cos_dec
    sig_ra2 *= cos_dec

    # Positional uncertainty per epoch, then combined for the separation
    sigma_pos1  = np.sqrt(sig_ra1**2  + sig_dec1**2)
    sigma_pos2  = np.sqrt(sig_ra2**2  + sig_dec2**2)
    sigma_theta = np.sqrt(sigma_pos1**2 + sigma_pos2**2) * u.arcsec

    # Linear separation uncertainty
    distance_au = (distance_pc * u.pc).to(u.AU)
    sigma_lin   = sigma_theta.to(u.rad).value * distance_au

    # Velocity and proper motion
    v_proj       = (lin_sep   / delta_t).to(u.km / u.s)
    sigma_v_proj = (sigma_lin / delta_t).to(u.km / u.s)
    pm           = (ang_sep   / delta_t).to(u.arcsec / u.yr)
    sigma_pm     = (sigma_theta / delta_t).to(u.arcsec / u.yr)

    cos_inc      = np.cos(np.deg2rad(inclination_deg))
    v_deproj     = v_proj       / cos_inc
    sigma_deproj = sigma_v_proj / cos_inc

    return dict(
        ang_sep_mas      = ang_sep.to(u.mas).value,
        ang_sep_err_mas  = sigma_theta.to(u.mas).value,
        lin_sep_au       = lin_sep.value,
        lin_sep_err_au   = sigma_lin.value,
        v_proj_kms       = v_proj.value,
        v_proj_err_kms   = sigma_v_proj.value,
        v_deproj_kms     = v_deproj.value,
        v_deproj_err_kms = sigma_deproj.value,
        pm_mas_yr        = pm.to(u.mas / u.yr).value,
        pm_err_mas_yr    = sigma_pm.to(u.mas / u.yr).value,
    )


def print_result(name: str, r: dict) -> None:
    print(f"{name}:")
    print(f"  Angular separation : {r['ang_sep_mas']:.3f} ± {r['ang_sep_err_mas']:.3f} mas")
    print(f"  Linear separation  : {r['lin_sep_au']:.3f} ± {r['lin_sep_err_au']:.3f} AU")
    print(f"  Projected velocity : {r['v_proj_kms']:.3f} ± {r['v_proj_err_kms']:.3f} km/s")
    print(f"  Deprojected veloc. : {r['v_deproj_kms']:.3f} ± {r['v_deproj_err_kms']:.3f} km/s")
    print(f"  Proper motion      : {r['pm_mas_yr']:.3f} ± {r['pm_err_mas_yr']:.3f} mas/yr")


# ---------------------------------------------------------------------------
# Source catalogue
# ---------------------------------------------------------------------------

SOURCES = dict(
    src_name = [
        'UYSO1 A', 'UYSO1 B', 'G11.11 A', 'G11.11 B'
    ],
    ra1 = [
        '07h05m10.94082s','07h05m10.81042s','18h10m28.39686s','18h10m28.33157s'
    ],
    err_ra1 = [
        0.00019,0.00027,0.00078,0.00094,0.000099
    ],
    dec1 = [
        '-12d19m00.47433s','-12d18m56.78213s','-19d22m29.97343s','-19d22m30.63570s'
    ],
    err_dec1 = [
        0.00527,0.01315,0.00980,0.01882
    ],
    date1 = [
        '24/11/2013','24/11/2013','20/03/2011','20/03/2011'
    ],
    ra2 = [
        '07h05m10.93861s','07h05m10.80961s','18h10m28.39523s','18h10m28.32875s'
    ],
    err_ra2 = [
        0.0000848,0.0002542,0.0002104,0.0005032
    ],
    dec2 = [
        '-12d19m00.4672s','-12d18m56.7640s','-19d22m29.9117s','-19d22m30.5491s'
    ],
    err_dec2 = [
        0.001604,0.004494,0.003883,0.009939
    ],
    date2 = [
        '14/05/2022','14/05/2022','25/04/2022','25/04/2022'
    ],
    # Distance in pc
    d = [
        1000,1000,2870,2870
    ],
    # Phase calibrator accuracy: A=2 mas, B=5 mas, C=10 mas (NRAO convention)
    pos_acc = [
        2,2,2,2
    ],
)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute proper motions and velocities from two-epoch VLA astrometry."
    )
    p.add_argument("--inclination", type=float, default=57.0,
                   help="Inclination angle in degrees (default: 57).")
    p.add_argument("--output", default=None,
                   help="Save results to this CSV file.")

    # Optional single-source mode
    single = p.add_argument_group("single-source mode (overrides built-in catalogue)")
    single.add_argument("--ra1");   single.add_argument("--dec1")
    single.add_argument("--date1"); single.add_argument("--ra2")
    single.add_argument("--dec2");  single.add_argument("--date2")
    single.add_argument("--distance", type=float)
    single.add_argument("--pos-acc", type=float, dest="pos_acc")
    single.add_argument("--err-ra1",  type=float, default=0.0, dest="err_ra1")
    single.add_argument("--err-dec1", type=float, default=0.0, dest="err_dec1")
    single.add_argument("--err-ra2",  type=float, default=0.0, dest="err_ra2")
    single.add_argument("--err-dec2", type=float, default=0.0, dest="err_dec2")
    single.add_argument("--name", default="custom source")
    return p.parse_args()


def run_catalogue(inclination_deg: float) -> list[dict]:
    cat = SOURCES
    results = []
    n = len(cat["src_name"])
    for i in range(n):
        name = cat["src_name"][i]
        try:
            r = proper_motion_calc(
                cat["ra1"][i],  cat["dec1"][i],  cat["date1"][i],
                cat["ra2"][i],  cat["dec2"][i],  cat["date2"][i],
                distance_pc     = cat["d"][i],
                pos_acc_mas     = cat["pos_acc"][i],
                inclination_deg = inclination_deg,
                err_ra1=cat["err_ra1"][i], err_dec1=cat["err_dec1"][i],
                err_ra2=cat["err_ra2"][i], err_dec2=cat["err_dec2"][i],
            )
            print_result(name, r)
            r["source"] = name
            r["date1"]  = cat["date1"][i]
            r["date2"]  = cat["date2"][i]
            results.append(r)
        except Exception as exc:
            print(f"[WARN] Skipping {name}: {exc}")
        print()
    return results


def save_csv(results: list[dict], path: str) -> None:
    df = pd.DataFrame(results)
    # Reorder for readability
    cols = ["source", "date1", "date2",
            "ang_sep_mas", "ang_sep_err_mas",
            "lin_sep_au",  "lin_sep_err_au",
            "v_proj_kms",  "v_proj_err_kms",
            "v_deproj_kms","v_deproj_err_kms",
            "pm_mas_yr",   "pm_err_mas_yr"]
    df[cols].to_csv(path, index=False, float_format="%.4f")
    print(f"Results saved → {path}")


def main():
    args = parse_args()

    single_mode = all(
        v is not None for v in
        [args.ra1, args.dec1, args.date1, args.ra2, args.dec2, args.date2,
         args.distance, args.pos_acc]
    )

    if single_mode:
        r = proper_motion_calc(
            args.ra1, args.dec1, args.date1,
            args.ra2, args.dec2, args.date2,
            distance_pc     = args.distance,
            pos_acc_mas     = args.pos_acc,
            inclination_deg = args.inclination,
            err_ra1=args.err_ra1, err_dec1=args.err_dec1,
            err_ra2=args.err_ra2, err_dec2=args.err_dec2,
        )
        print_result(args.name, r)
        if args.output:
            r["source"] = args.name
            save_csv([r], args.output)
    else:
        results = run_catalogue(args.inclination)
        if args.output:
            save_csv(results, args.output)


if __name__ == "__main__":
    main()