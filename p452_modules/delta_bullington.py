# 4.2.3 Complete 'delta-Bullington' diffraction model

from p452_modules import diff_spherical_earth
import numpy as np


def delta_bull(d, h, hts, hrs, hstd, hsrd, a, f, omega, pol):
    
    # Use the method in 4.2.1 for the actual terrain profile and antenna
    # heights. Set the resulting Bullington diffraction loss for the actual
    # path to Lbulla
    
    Lbulla = bull_loss(d, h, hts, hrs, a, f)

    # Use the method in § 4.2.1 for a second time, with all profile heights,
    # h i , set to zero, and modified antenna heights given by eqs 38:
    
    h0 = np.zeros(np.shape(h))
    hts0 = hts - hstd
    hrs0 = hrs - hsrd
    
    Lbulls = bull_loss(d, h0, hts0, hrs0, a, f)

    # Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
    # for the actual path length dtot (eq 39):
    dtot = d[-1] - d[0]
    hte = hts0           
    hre = hrs0
    
    Lbull_sph = diff_spherical_earth.sph_earth_loss(dtot, hte, hre, a, f, omega, pol)

    # Diffraction loss for the general path is now given by equation 40:
    Ld = Lbulla + max(Lbull_sph - Lbulls, 0)

    return Ld


def bull_loss(d, h, hts, hrs, a, f):
    # Bullington diffraction loss (section 4.2.1)

    dtot = d[-1] - d[0]
    C = 1 / a
    lam = 0.3 / f

    di = d[1:-1]
    hi = h[1:-1]

    Stim = max((h[1:-1] + 500 * C * d[1:-1] * (dtot - d[1:-1]) - hts) / d[1:-1])
    Str = (hrs - hts) / dtot

    if Stim < Str:
        # Los Path

        # Find the intermediate profile point with the highest diffraction
        # parameter v:

        vmax = np.max((hi + 500 * C * di * (dtot - di) - (hts * (dtot - di) + hrs * di) / dtot) *
                      np.sqrt(0.002 * dtot / (lam * di * (dtot - di))))

        # apply eq 17:
        if vmax > -0.78:
            Luc = 6.9 + 20 * np.log10(np.sqrt((vmax - 0.1) ** 2 + 1) + vmax - 0.1)
        else:
            Luc = 0

    else:

        # Trans-horizon path

        Srim = np.max((hi + 500 * C * di * (dtot - di) - hrs) / (dtot - di))

        # distance of the Bullington point from the transmitter # eq 19:

        dbp = (hrs - hts + Srim * dtot) / (Stim + Srim)

        # diffraction parameter vb (eq 20)

        vb = (hts + Stim * dbp - (hts * (dtot - dbp) + hrs * dbp) / dtot) * np.sqrt(
            0.002 * dtot / (lam * dbp * (dtot - dbp)))

        # The knife-edge loss for the Bullington point is ( eq 21)

        # apply eq 22
        if vb > -0.78:
            Luc = 6.9 + 20 * np.log10(np.sqrt((vb - 0.1) ** 2 + 1) + vb - 0.1)
        else:
            Luc = 0

    # eq 22:
    Lbull = Luc + (1 - np.exp(-Luc / 6.0)) * (10 + 0.02 * dtot)

    return Lbull
