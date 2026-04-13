'''
This script calculates the distance between two positions.
You can get the velocity from going to point A to point B as well for a given inclination.
Error propagation is done taking into account both positional accuracy and individual measurement errors.
'''

from astropy.coordinates import SkyCoord
from astropy import units as u
from datetime import datetime
import numpy as np


def proper_motion_calc(
    ra1, dec1, date1, ra2, dec2, date2, distance_pc, pos_acc, inclination_deg,
    err_ra1=0, err_dec1=0, err_ra2=0, err_dec2=0):
    """
    Compute angular separation, linear separation, and deprojected velocity with error estimates.

    Parameters:
    - ra1, dec1, ra2, dec2: RA/Dec (hh:mm:ss.sss dd:mm:ss.ss)
    - date1, date2: Observing dates in 'dd/mm/yyyy'
    - distance_pc: Source distance (pc)
    - inclination_deg: Inclination angle (deg)
    - pos_acc: Phase calibrator positional uncertainty (mas)
    - err_ra1, err_ra2: RA measurement errors (arcsec)
    - err_dec1, err_dec2: Dec measurement errors (arcsec)

    Returns:
    - Angular and linear separation in the sky (proj)
    - Projected velocity
        Optional:
        - Time between observations in days
        - Deprojected velocity
    """
    # Parse time with units
    t1 = datetime.strptime(date1, '%d/%m/%Y')
    t2 = datetime.strptime(date2, '%d/%m/%Y')
    delta_t = abs((t2 - t1).total_seconds()) * u.s
    days = delta_t.to(u.day).value

    if delta_t == 0 * u.s:
        raise ValueError("Dates must not be the same.")

    # SkyCoords
    coord1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg), distance=distance_pc * u.pc)
    coord2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg), distance=distance_pc * u.pc)

    # Angular and linear separation
    angular_sep = coord1.separation(coord2)
    linear_sep = coord1.separation_3d(coord2).to(u.AU)

    # Convert individual measurement errors to common units (arcsec)
    # RA errors are already in arcsec, Dec errors are in arcsec
    err_ra1_arcsec = err_ra1 * u.arcsec
    err_dec1_arcsec = err_dec1 * u.arcsec
    err_ra2_arcsec = err_ra2 * u.arcsec
    err_dec2_arcsec = err_dec2 * u.arcsec
    
    # Phase calibrator uncertainty in arcsec
    pos_acc_arcsec = pos_acc * u.mas
    pos_acc_arcsec = pos_acc_arcsec.to(u.arcsec)
    
    # Combine errors quadratically for each epoch
    # Total positional error = sqrt(individual_measurement_error^2 + phase_calibrator_error^2)
    sigma_ra1_total = np.sqrt(err_ra1_arcsec.value**2 + pos_acc_arcsec.value**2) * u.arcsec
    sigma_dec1_total = np.sqrt(err_dec1_arcsec.value**2 + pos_acc_arcsec.value**2) * u.arcsec
    sigma_ra2_total = np.sqrt(err_ra2_arcsec.value**2 + pos_acc_arcsec.value**2) * u.arcsec
    sigma_dec2_total = np.sqrt(err_dec2_arcsec.value**2 + pos_acc_arcsec.value**2) * u.arcsec
    
    # Convert RA errors to true angular separation (multiply by cos(dec))
    # Use average declination for the conversion
    dec1_rad = coord1.dec.radian
    dec2_rad = coord2.dec.radian
    avg_dec = (dec1_rad + dec2_rad) / 2
    cos_dec = np.cos(avg_dec)
    
    sigma_ra1_angular = sigma_ra1_total * cos_dec
    sigma_ra2_angular = sigma_ra2_total * cos_dec
    
    # Total angular uncertainty for separation measurement
    # Propagate errors: sigma_separation = sqrt(sigma_pos1^2 + sigma_pos2^2)
    # For each position: sigma_pos = sqrt(sigma_ra^2 + sigma_dec^2)
    sigma_pos1 = np.sqrt(sigma_ra1_angular.value**2 + sigma_dec1_total.value**2) * u.arcsec
    sigma_pos2 = np.sqrt(sigma_ra2_angular.value**2 + sigma_dec2_total.value**2) * u.arcsec
    
    # Total uncertainty in angular separation
    sigma_theta = np.sqrt(sigma_pos1.value**2 + sigma_pos2.value**2) * u.arcsec
    sigma_theta_rad = sigma_theta.to(u.rad)

    # Convert angular uncertainty to linear (AU)
    distance_au = (distance_pc * u.pc).to(u.AU)
    sigma_s = sigma_theta_rad.value * distance_au

    # Velocity and error (projected)
    v_proj = (linear_sep / delta_t).to(u.km/u.s)
    sigma_v_proj = (sigma_s / delta_t).to(u.km/u.s)

    # Proper motion
    pm = (angular_sep / delta_t).to(u.arcsec/u.yr)
    sigma_pm = (sigma_theta / delta_t).to(u.arcsec/u.yr)
    
    # Deprojected velocity
    inclination_rad = np.deg2rad(inclination_deg)
    v_deproj = v_proj / np.cos(inclination_rad)
    sigma_v_deproj = sigma_v_proj / np.cos(inclination_rad)

    # Output with enhanced error information
    print(f"Angular separation     : {angular_sep.to(u.mas).value:.3f} mas ± {sigma_theta.to(u.mas).value:.3f} mas")
    print(f"Linear separation      : {linear_sep:.3f} ± {sigma_s:.3f}")
    print(f"Projected velocity     : {v_proj:.3f} ± {sigma_v_proj:.3f}")
    print(f"Proper motion          : {pm.to(u.mas/u.yr).value:.3f} ± {sigma_pm.to(u.mas/u.yr).value:.3f} mas/yr")
    # print(f"Error breakdown:")
    # print(f"  Phase cal. error     : {pos_acc:.1f} mas")
    # print(f"  Epoch 1 total error  : {sigma_pos1.to(u.mas).value:.3f} mas")
    # print(f"  Epoch 2 total error  : {sigma_pos2.to(u.mas).value:.3f} mas")
    # # print(f"Time span              : {days:.2f} days")
    # # print(f"Inclination angle      : {inclination_deg:.1f}°")
    # # print(f"Deprojected velocity   : {v_deproj:.3f} ± {sigma_v_deproj:.3f}")

    return(angular_sep.to(u.arcsec).value, linear_sep.value, v_proj.value, v_deproj.value, pm.value)

# === INPUT ARRAYS ===

src_name = ['UYSO1 A', 'UYSO1 B', 'G11.11 A', 'G11.11 B', '18151 B', 'G23.01 A','18470 B',\
    '18517 A','18521 A','18440 A','19012 A','G35.39 A','G53.25 mm2', 'G53.25mm4',\
    '20343 B','20293 E',\
    'G23.01 A - C band','G53.25 mm2 - C band','G53.25mm4 - C band','20343 B - C band'] # this line here is to compare Rosero's C and K band

## Rosero+16,19 data
ra1 = ['07h05m10.94082s','07h05m10.81042s','18h10m28.39686s','18h10m28.33157s','18h17m58.123256s','18h34m40.28520s','18h49m37.74987s',\
    '18h54m14.24390s','18h54m40.74066s','18h46m36.60419s','19h03m45.26661s','18h56m59.06148s','19h29m33.52824s','19h29m34.20071s',
    '20h36m07.52489s','20h31m12.89296',\
    '18h34m40.28520s','19h29m33.52824s','19h29m34.20071s','20h36m07.52489s']
err_ra1 = [0.00019,0.00027,0.00078,0.00094,0.000099,0.00018,0.00019,\
    0.00032,0.00021,0.00015,0.00023,0.00038,0.00060,0.00044,\
    0.00080,0.00081,\
    0.00018,0.00060,0.00044,0.00080]# this line here is to compare Rosero's C and K band
dec1 = ['-12d19m00.47433s','-12d18m56.78213s','-19d22m29.97343s','-19d22m30.63570s','-12d07m24.724407s','-09d00m38.29499s','-00d41m00.19731s',\
    '04d41m40.87830s','01d38m06.54255s','-01d45m22.63078s','05d40m42.75203s','02d04m54.56150s','18d00m54.21179s','18d01m39.114055s',\
    '41d40m09.04707','40d03m22.71089s',\
    '-09d00m38.29499s','18d00m54.21179s','18d01m39.114055s','41d40m09.04707']# this line here is to compare Rosero's C and K band
err_dec1 = [0.00527,0.01315,0.00980,0.01882,0.002459,0.00423,0.00380,\
    0.00765,0.00381,0.00321,0.00479,0.00704,0.00663,0.00565,\
    0.00296,0.00443,\
    0.00423,0.00663,0.00565,0.00296]# this line here is to compare Rosero's C and K band
date1 = ['24/11/2013','24/11/2013','20/03/2011','20/03/2011','20/03/2011','18/01/2014','02/05/2011',\
    '14/04/2011','14/04/2011','02/05/2011','15/04/2011','14/04/2011','17/01/2014','17/01/2014',
    '18/01/2014','18/01/2014',\
    '18/01/2014','17/01/2014','17/01/2014','18/01/2014']# this line here is to compare Rosero's C and K band

## VLA/22A-092
ra2 = ['07h05m10.93861s','07h05m10.80961s','18h10m28.39523s','18h10m28.32875s','18h17m58.12329s','18h34m40.28342s','18h49m37.74770s',\
    '18h54m14.24301s','18h54m40.73879s','18h46m36.60282s','19h03m45.26517s','18h56m59.05763s','19h29m33.52532s','19h29m34.20008s',\
    '20h36m07.52084s','20h31m12.88999s',\
    '18h34m40.28649s','19h29m33.52619s','19h29m34.20447s','20h36m07.52334s']# this line is Rosero's C band; 'G23.01 A','G53.25 mm2','G53.25mm4','20343 B','20293 E'
err_ra2 = [0.0000848,0.0002542,0.0002104,0.0005032,0.0000224,0.0000902,0.0001003,\
0.0001020,0.0000913,0.0001887,0.0000478,0.0001158,0.0001765,0.0003187,\
0.0006730,0.0005295,\
0.00013,0.00034,0.00029,0.00066]# this line is Rosero's C band
dec2 = ['-12d19m00.4672s','-12d18m56.7640s','-19d22m29.9117s','-19d22m30.5491s','-12d07m24.7754s','-09d00m38.3314s','-00d41m00.2466s',\
    '04d41m40.8391s','01d38m06.4614s','-01d45m22.6980s','05d40m42.6856s','02d04m54.4730s','18d00m54.1546s','18d01m39.0826s',\
    '41d40m09.0160s','40d03m22.6881s',\
    '-9d00m38.28183s','18d00m54.21661s','18d01m39.13154s','41d40m09.07713s']# this line is Rosero's C band
err_dec2 = [0.001604,0.004494,0.003883,0.009939,0.000421,0.001395,0.001658,\
0.001950,0.001545,0.003008,0.000743,0.001883,0.002413,0.004354,\
0.007122,0.004689,0.00916,\
0.00346,0.00916,0.00396,0.00761]# this line is Rosero's C band
date2 = ['14/05/2022','14/05/2022','25/04/2022','25/04/2022','18/05/2022','06/05/2022','09/04/2022',\
    '09/04/2022','09/04/2022','07/03/2022','30/04/2022','27/03/2022','05/05/2022','05/05/2022',\
    '06/06/2022','06/06/2022',\
    '11/08/2011','28/08/2011','28/08/2011','07/08/2011']# this line is Rosero's C band

d = [1000,1000,2870,2870,1990,4590,6530,1900,9100,5200,4200,2300,2000,2000,1400,2000,\
4590,2000,2000,1400]# this line here is to compare Rosero's C and K band

# Phase calibrator astrometry precision. I took: A = 2 mas, B = 5 mas, C = 10 mas.
# In the NRAO manual this is: A < 2 mas, B = 2-10 mas, C = 10-150 mas.
pos_acc = [2,2,2,2,10,10,10,10,10,10,10,10,2,2,5,5,\
10,2,2,5]# this line here is to compare Rosero's C and K band

ang_sep = []
lin_sep = []
vel = []
vel_deproj = []
proper_motion = []

for i in range(len(src_name)):
    print(f"{src_name[i]}:")
    aux_ang_sep, aux_lin_sep, aux_vel, aux_vel_deproj, aux_proper_motion = proper_motion_calc(
        ra1[i], dec1[i], date1[i], ra2[i], dec2[i], date2[i], d[i], pos_acc[i], 
        inclination_deg=57, err_ra1=err_ra1[i], err_dec1=err_dec1[i], 
        err_ra2=err_ra2[i], err_dec2=err_dec2[i])
    
    ang_sep = np.append(ang_sep, aux_ang_sep)
    lin_sep = np.append(lin_sep, aux_lin_sep)
    vel = np.append(vel, aux_vel)
    vel_deproj = np.append(vel_deproj, aux_vel_deproj)
    proper_motion = np.append(proper_motion, aux_proper_motion*1000)

    print('\n')

# Optional summary table (uncomment if needed)
# print('SOURCE \t\t ANG SEP \t\t LIN SEP \t\t Obs date 1 \t Obs date 2 \t\t V_proj \t\t PM (mas/yr)')
# print('\t\t (mas) \t\t\t (au)')
# for j in range(len(src_name)):
#     print(f"{src_name[j]}\t{ang_sep[j]*1000:.3f}\t{lin_sep[j]:.3f}\t{date1[j]}\t{date2[j]}\t{vel[j]:.3f}\t{proper_motion[j]:.3f}")