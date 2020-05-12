# This function gives the prediction of the basic transmission loss, Lba (dB) occurring during periods of anomalous
# propagation (ducting and layer reflection) as given in Section 4.4 eq 46-56
import numpy as np


def duct_loss(freq, p, dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre,
      hm, theta_t, theta_r, omega, ae, b0, Ag):

    # site-shielding diffraction losses for the interfering and interfered-with stations
    # respectively eq (48a) and (48):
    theta_t1 = theta_t - 0.1 * dlt
    theta_r1 = theta_r - 0.1 * dlr

    if theta_t1 > 0:
        Ast = 20 * np.log10(1 + 0.361 * theta_t1 * np.sqrt(freq * dlt)) + 0.264 * theta_t1 * freq ** (1/3)
    else:
        Ast = 0

    if theta_r1 > 0:
        Asr = 20 * np.log10(1 + 0.361 * theta_r1 * np.sqrt(freq * dlr)) + 0.264 * theta_r1 * freq ** (1 / 3)
    else:
        Asr = 0

    # empirical correction to account for the increasing attenuation with wavelength
    # inducted propagation, eq 47a
    if freq < 0.5:
        Alf = 45.375 - 137.0 * freq + 92.5 * freq**2
    else:
        Alf = 0

    # over-sea surface duct coupling corrections for the interfering and interferedwith stations respectively, eq 49:
    if (dct <= 5) and (dct <= dlt) and (omega >= 0.75):
        Act = -3 * np.exp(-0.25 * dct * dct) * (1 + np.tanh(0.07 * (50 - hts)))
    else:
        Act = 0

    if (dcr <= 5) and (dcr <= dlr) and (omega >= 0.75):
        Acr = -3 * np.exp(-0.25 * dcr * dcr) * (1 + np.tanh(0.07 * (50 - hrs)))
    else:
        Acr = 0

    # total of fixed coupling losses (except for local clutter losses) between the
    # antennas and the anomalous propagation structure within the atmosphere, Af, eq (47):
    Af = 102.45 + 20 * np.log10(freq) + 20 * np.log10(dlt + dlr) + Alf + Ast + Asr + Act + Acr

    # specific attenuation eq (51)
    g_d = 5e-5 * ae * freq ** (1/3)

    # angular distance (corrected where appropriate (via equation (52a)) to allow for
    # the application of the site shielding model in equation (48)):
    theta_t1 = min(theta_t, 0.1 * dlt)
    theta_r1 = min(theta_r, 0.1 * dlr)

    # eq 52:
    theta1 = 1e3 * dtot / ae + theta_t1 + theta_r1

    # eq 56a
    dI = min(dtot - dlt - dlr, 40)

    if hm > 10:
        mu3 = np.exp(-4.6e-5 * (hm - 10) * (43 + 6 * dI))   # eq (56)
    else:
        mu3 = 1

    # tau as defined in eq 2:
    tau = 1 - np.exp(-(4.12e-4 * dlm ** 2.41))

    # alpha parameter as in eq 55a
    alpha = -0.6 - 3.5 * 1e-9 * (dtot ** 3.1) * tau
    if alpha < -3.4:
        alpha = -3.4

    # μ2: correction for path geometry, eq (55)
    mu2 = (500 / ae * dtot ** 2 / (np.sqrt(hte) + np.sqrt(hre)) ** 2) ** alpha
    if mu2 > 1:
        mu2 = 1

    # eq (54):
    beta = max(b0 * mu2 * mu3, 1e-10)

    # eq 53a:
    G = 1.076 / (2.0058 - np.log10(beta)) ** 1.012 * np.exp(-(9.51 - 4.8 * np.log10(beta) +
                                                              0.198 * (np.log10(beta)) ** 2) * 1e-6 * dtot ** 1.13)

    # time percentage variability (cumulative distribution), eq 53:
    Ap = -12 + (1.2 + 3.7e-3 * dtot) * np.log10(p / beta) + 12 * (p / beta) ** G

    # time percentage and angular-distance dependent losses within the anomalous
    # propagation mechanism, Adp, eq (50):
    Adp = g_d * theta1 + Ap


    # basic transmission loss, Lba (dB) occurring during periods of anomalous
    # propagation (ducting and layer reflection) is based on the following function, eq (46):
    Lba = Af + Adp + Ag

    return Lba
