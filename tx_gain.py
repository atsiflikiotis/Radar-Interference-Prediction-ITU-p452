import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline


def tx_gain(atrdeg, theta_t, antennas_db, antenna, azim, eltilt, mechtilt):

    e_pt = theta_t / 1000
    # in P.452-16 eq 69a doesnt seem right, and since theta_t is defined in function
    # "path_parameters" fo both los and trans-horizon cases, theta_t (mrad) / 1000 is used instead

    # elevation transmitter angle
    e_pt_deg = -np.rad2deg(e_pt)  # - sign is used because in most vertical antenna patterns, dowtilt is positive

    # antenna_gain = # get antenna gain from pandas database, from antenna name and el tilt

    anglehor = atrdeg - azim
    anglever = e_pt_deg - mechtilt

    if anglehor < 0:
        anglehor += 360

    if anglever < 0:
        anglever += 360

    # get antenna pattern
    row = antennas_db.loc[(antenna, eltilt)]
    maxgain = row['Gain'].item()
    horizontal = row['Pattern'].item()[:, 0]
    vertical = row['Pattern'].item()[:, 1]

    x = np.linspace(0, 360, 361)
    fhor = InterpolatedUnivariateSpline(x, horizontal)
    fver = InterpolatedUnivariateSpline(x, vertical)

    gain = maxgain - fhor(anglehor) - fver(anglever)

    return gain

