# Input parameters:
# f:        Frequency (GHz)
# p:        Required time percentage for which the calculated basic transmission loss is not exceeded
# d:        array of distances di of the i:th profile point (km)
# h:        array of heights hi of the i-th profile point (meters) above mean sea level.
#           Both d and h vectors contain n+1 profile points
# zone:     Zone type: Coastal land (A1), Inland (A2) or Sea (B) ((n+1) size numpy array)
# htg:      Tx Antenna mid height above ground level (m)
# hrg:      Rx Antenna mid height above ground level (m)
# phi_t:    Latitude of Tx station (decimal degrees)
# phi_r:    Latitude of Rx station (decimal degrees)
# Gt:       Tx Antenna gain in the direction of the horizon along the greatcircle interference path (dBi)
# Gr:       Rx Antenna gain in the direction of the horizon along the greatcircle interference path (dBi)
# pol:      polarization (1) horizontal, (2) vertical
# dct, dcr: Distance over land from the transmit and receive antennas to the coast along the great-circle
#           interference path (km).
# DN:       Average radio:refractive index
# N0:       sea-level surface refractivity
# pressure: Dry air pressure (hPa)
# temp:     Air temperature (°C)

# kwargs:
# ha_t:     Clutter nominal height (m) at the Tx side
# ha_r:     Clutter nominal height (m) at the Rx side
# dk_t:     Clutter nominal distance (km) at the Tx side
# dk_r:     Clutter nominal distance (km) at the Rx side
#
# examples:
#  p452_loss(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp)
#
#  p452_loss(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp,
#            ha_t, ha_r, dk_r, dk_t)
#
# Output:
# Lb:       Basic  transmission loss according to P.452-16

import itertools
import operator

import numpy as np

import srtm_path
from p452_modules import p676_12, diff_loss, clutter_corr, path_parameters
from p452_modules.ducting_layer_loss import duct_loss as duct_layer
from p452_modules.inv_cumulative import inv_cumulative as inv_cumul
from tx_gain import tx_gain


def p452_loss(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gr, pol, dct, dcr, DN, N0, pressure, temp,
              Gt=None, psi_t=None, psi_r=None, antennasdb=None, antennaname=None, azim=None, eltilt=None, mechtilt=None,
              ha_t=0, dk_t=0, ha_r=0, dk_r=0):

    phi_center = (phi_t + phi_r) / 2
    dtot = d[-1] - d[0]
    T = 273.15 + temp

    # longest continuous land (inland + coastal) section of the great-circle path (km)

    A2 = (zone == 'A2')
    A1A2 = np.logical_or(zone == 'A1', zone == 'A2')
    B = (zone == 'B')

    idx_A1A2 = [[i for i, value in it] for key, it in itertools.groupby(enumerate(A1A2), key=operator.itemgetter(1)) if
                key]
    idx_A2 = [[i for i, value in it] for key, it in itertools.groupby(enumerate(A2), key=operator.itemgetter(1)) if key]
    idx_B = [[i for i, value in it] for key, it in itertools.groupby(enumerate(B), key=operator.itemgetter(1)) if key]

    nA1A2 = len(idx_A1A2)
    nA2 = len(idx_A2)
    nB = len(idx_B)

    dtm = dlm = dsea = 0
    for i in range(0, nA1A2):
        startidx = idx_A1A2[i][0]
        stopidx = idx_A1A2[i][-1]
        dd = 0
        if d[stopidx] < d[-1]:
            dd = dd + (d[stopidx + 1] - d[stopidx]) / 2
        if d[startidx] > 0:
            dd = dd + (d[stopidx] - d[stopidx - 1]) / 2

        dtm = max(d[stopidx] - d[startidx] + dd, dtm)

    for i in range(0, nA2):
        startidx = idx_A2[i][0]
        stopidx = idx_A2[i][-1]
        dd = 0
        if d[stopidx] < d[-1]:
            dd = dd + (d[stopidx + 1] - d[stopidx]) / 2
        if d[startidx] > 0:
            dd = dd + (d[stopidx] - d[stopidx - 1]) / 2

        dlm = max(d[stopidx] - d[startidx] + dd, dlm)

    # b0 calculation
    tau = 1 - np.exp(-(4.12 * 1e-4 * dlm ** 2.41))
    mu1 = (10 ** (-dtm / (16 - 6.6 * tau)) + 10 ** (-5 * (0.496 + 0.354 * tau))) ** 0.2

    if mu1 > 1:
        mu1 = 1

    if abs(phi_center) < 70:
        mu4 = 10 ** ((-0.935 + 0.0176 * abs(phi_center)) * np.log10(mu1))
        b0 = 10 ** (-0.015 * abs(phi_center) + 1.67) * mu1 * mu4
    else:
        mu4 = 10 ** (0.3 * np.log10(mu1))
        b0 = 4.17 * mu1 * mu4

    # median value of the effective earth radius,
    # and the effective Earth radius exceeded for beta0% of time (eq. 5-6)
    ae = 6371 * 157 / (157 - DN)
    kb = 3
    ab = 6371 * kb
    Ce = 1 / ae

    # path fraction over sea (zone == A3)
    if any(B):
        for i in range(nB):
            startidx = idx_B[i][0]
            stopidx = idx_B[i][-1]
            dd = 0
            if d[stopidx] < d[-1]:
                dd = dd + (d[stopidx + 1] - d[stopidx]) / 2
            if d[startidx] > 0:
                dd = dd + (d[stopidx] - d[stopidx - 1]) / 2

            dsea = dsea + d[stopidx] - d[startidx] + dd

    omega = dsea / dtot

    if (ha_t != 0) or (ha_t != 0):
        dcorr, hcorr, htgcorr, hrgcorr, Aht, Ahr = clutter_corr.clutterloss(f, d, h, htg, hrg, ha_t, ha_r, dk_t, dk_r)

        # replace values with new values after clutter corrections
        d = dcorr
        h = hcorr
        htg = htgcorr
        hrg = hrgcorr
    else:
        Aht = 0
        Ahr = 0

    # Total Tx and Rx antenna heights (above mean sea level)
    hts = h[0] + htg
    hrs = h[-1] + hrg

    # path line parameters
    hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype = \
        path_parameters.path_parameters(d, h, htg, hrg, ae, f)

    # Get Tx antenna gain, from antennas db, after calculating bearing angle phi(tx,rx) receiver,
    # if None is passed as Gt parameter:
    if Gt is None:
        atrdeg = srtm_path.get_path_geometry(phi_t, psi_t, phi_r, psi_r, 0)[1]
        Gt = tx_gain(atrdeg, theta_t, antennasdb, antennaname, azim, eltilt, mechtilt)

    #
    # 4.1 line of sight transmission loss including short-term effects
    # water vapor density eq 9a
    rho = 7.5 + 2.5 * omega

    # Specific attenuation due to dry air and water vapour, ITU-R P.676-12
    g_0, g_w = p676_12.specatten(f, pressure, rho, T)

    # total gaseous absorption (eq. 9)
    Ag = (g_0 + g_w) * dtot

    # Basic transmission loss due to free-space propagation and
    # attenuation by atmospheric gases (eq 8)
    Lbfsg = 92.5 + 20 * np.log10(f) + 20 * np.log10(dtot) + Ag

    # Corrections for multipath and focusing effects at p and b0
    # percentage  times (eq 10)

    Esp = 2.6 * (1 - np.exp(-0.1 * (dlt + dlr))) * np.log10(p / 50)
    Esb = 2.6 * (1 - np.exp(-0.1 * (dlt + dlr))) * np.log10(b0 / 50)

    # Basic transmission loss (dB) not exceeded for time percentage p%, due to
    # LoS propagation eq 11
    Lb0p = Lbfsg + Esp

    # Basic transmission loss (dB) not exceeded for time percentage b0%, due to
    # LoS propagation eq 12
    Lb0b = Lbfsg + Esb

    # 4.2 Diffraction loss

    # intermediate profile point with the highest slope of the line from the transmitter to the
    # point (eq 14)
    Stim = max((h[1:-1] + 500 * Ce * d[1:-1] * (dtot - d[1:-1]) - hts) / d[1:-1])

    # slope of the line from transmitter to receiver assuming a LoS path
    Str = (hrs - hts) / dtot

    # Ldp: The diffraction loss not exceeded for p% time (4.2.4)
    Ldp, Ld50 = diff_loss.diffloss(d, h, hts, hrs, hstd, hsrd, f, omega, p, b0, ae, ab, pol)

    # The median basic transmission loss associated with diffraction is guben by eq 43:
    Lbd50 = Lbfsg + Ld50

    # The basic tranmission loss associated with diffraction not exceeded for p% time is givenn by eq 44:
    Lbd = Lb0p + Ldp

    # 4.3 Tropospheric scatter
    # Frequency dependent loss, eq 45a
    Lf = 25 * np.log10(f) - 2.5 * (np.log10(f / 2)) ** 2

    # aperture to medium coupling loss (dB), eq 45b
    Lc = 0.051 * np.exp(0.055 * (Gt + Gr))

    # Ag, gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
    # whole path length, as given in eq (9):
    g_0, g_w = p676_12.specatten(f, pressure, 3, T)
    Ag_3 = (g_0 + g_w) * dtot

    # The basic transmission loss due to troposcatter, Lbs (dB) not exceeded for any time percentage, p,
    # below 50%, is given by eq (45):
    Lbs = 190 + Lf + 20 * np.log10(dtot) + 0.573 * theta_tot - 0.15 * N0 + Lc + Ag_3 - 10.1 * (-np.log10(p / 50)) ** 0.7

    #
    # 4.4 Ducting/layer reflection
    Lba = duct_layer(f, p, dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre,
                     hm, theta_t, theta_r, omega, ae, b0, Ag)

    #
    #
    # 4.6 The overall prediction section
    #
    # Calculate an interpolation factor Fj to take account of the path angular distance (eq 58)
    Fj = 1 - 0.5 * (1 + np.tanh(3 * 0.8 * (Stim - Str) / 0.3))

    # interpolation factor Fk to take account of the great circle path distance (eq 59)
    Fk = 1 - 0.5 * (1.0 + np.tanh(3.0 * 0.5 * (dtot - 20) / 20))

    # Calculate a notional minimum basic transmission loss, Lminb0p (dB) associated with LoS propagation
    # and over-sea sub-path diffraction (eq 60):

    if p < b0:
        Lminb0p = Lb0p + (1 - omega) * Ldp
    else:
        # apply eq 41a and eq 60
        Fi = inv_cumul(p / 100) / inv_cumul(b0 / 100)
        Lminb0p = Lbd50 + (Lb0b + (1 - omega) * Ldp - Lbd50) * Fi

    # Calculate a notional minimum basic transmission loss, Lminbap (dB), associated with LoS and
    # transhorizon signal enhancements, eq (61)
    # eta = 2.5
    Lminbap = 2.5 * np.log(np.exp(Lba / 2.5) + np.exp(Lb0p / 2.5))

    # Calculate a notional basic transmission loss, Lbda (dB), associated with diffraction and LoS or
    # ducting/layer-reflection enhancements, eq (62):
    if Lminbap <= Lbd:
        Lbda = Lminbap + (Lbd - Lminbap) * Fk
    else:
        Lbda = Lbd

    # Calculate a modified basic transmission loss, Lbam (dB), which takes diffraction and LoS or
    # ducting/layer-reflection enhancements into account, eq (63):
    Lbam = Lbda + (Lminb0p - Lbda) * Fj

    #
    # Calculate the final basic transmission loss not exceed for p% time, Lb (dB), as given by eq (64):
    Lb = -5 * np.log10(10 ** (-0.2 * Lbs) + 10 ** (-0.2 * Lbam)) + Aht + Ahr

    # 4.7 Calculation of transmission loss
    # The following procedure provides a method for the calculation of transmission loss between two
    # terrestrial stations. As intermediate steps in the method, it also provides formulae for the calculation
    # of the great-circle path length and angular distance based on the stations’ geographic coordinates

    return Lb, Lbfsg, Lb0p, Lb0b, Ld50, Ldp, Lbs, Lba, theta_t, Gt, pathtype
