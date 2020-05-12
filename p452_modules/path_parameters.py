import numpy as np

def path_parameters(d, h, htg, hrg, ae, f):
    # derive smooth-earth antenna geights according to section 4: Path classification and Section 5 of ITU P.452-16

    dtot = d[-1] - d[0]
    N = len(d)
    
    hts = h[0] + htg
    hrs = h[-1] + hrg

    di = d[1:N-1]
    hi = h[1:N-1]
    # theta = np.zeros(len(d) - 2)

    # Path classification
    # theta in mrads
    theta = 1000 * np.arctan((hi - hts) / (1000 * di) - di / (2 * ae))

    theta_t = max(theta)
    theta_td = 1000 * np.arctan((hrs - hts) / (1000 * dtot) - dtot / (2 * ae))
    theta_rd = 1000 * np.arctan((hts - hrs) / (1000 * dtot) - dtot / (2 * ae))  # not included in ITU-R P.452-16

    if theta_t > theta_td:
        pathtype = 'trans-horizon'
    else:
        pathtype = 'los'

    # Interfering antenna horizon distance
    ltindex = np.argwhere(theta == theta_t).item() + 1
    dlt = d[ltindex]

    # Interfered-with antenna horizon elevation angle eq 157
    theta = 1000 * np.arctan((hi - hrs) / (1000 * (dtot - di)) - (dtot - di) / (2 * ae))
    theta_r = max(theta)

    # Interfered-with antenna horizon distance eq 158
    # if more than 1 maxmimum indices found, must use the largest index as it is closer to receiver
    lrindex = np.argwhere(theta == theta_r)[-1].item() + 1
    dlr = dtot - d[lrindex]

    if pathtype == 'los':
        theta_t = theta_td
        theta_r = theta_rd
        # # # #
        lam = 0.3 / f
        Ce = 1 / ae

        nu = (hi + 500 * Ce * di * (dtot-di) - (hts * (dtot - di) + hrs * di) / dtot) * \
            np.sqrt(0.002 * dtot / (lam * di * (dtot-di)))

        numax = max(nu)
    
        ltindex = np.argmax(nu == numax) + 1
        dlt = d[ltindex]
        dlr = dtot - dlt
        lrindex = np.argwhere(dlr <= dtot - di)[-1].item() + 1
        # # #

    # angular distance (eq 159)
    theta_tot = 1e3 * dtot / ae + theta_t + theta_r

    v1 = 0
    for k in range(1, N):
        v1 = v1 + (d[k]-d[k-1]) * (h[k] + h[k-1])  
        
    v2 = 0
    for k in range(1, N):
        v2 = v2 + (d[k]-d[k-1])*(h[k] * (2*d[k] + d[k-1]) + h[k-1] * (d[k] + 2*d[k-1]))

    # height of the smooth-Earth surface at the interfering station (eq 163)
    hst = (2 * v1 * dtot - v2) / dtot**2

    # height of the smooth-Earth surface at the interfered-with station (eq 164)
    hsr = (v2 - v1 * dtot) / dtot**2          

    H = h - (hts * (dtot - d) + hrs * d) / dtot  
    
    h_obs = max(H[1:N-1])

    # horizon elevation angles:
    a_obt = max(H[1:N-1] / d[1:N-1])
    a_obr = max(H[1:N-1] / (dtot - d[1:N-1]))

    # provisional values for the heights of the smooth surface at the transmitter and receiver
    # ends of the path eq 166
    
    gt = a_obt/(a_obt + a_obr)
    gr = a_obr/(a_obt + a_obr)

    # eq 166:
    if h_obs <= 0:
        hstp = hst
        hsrp = hsr
    else:
        hstp = hst - h_obs*gt
        hsrp = hsr - h_obs*gr

    # final values for diffraction model (eqs 167)

    if hstp >= h[0]:
        hstd = h[0]
    else:
        hstd = hstp
    
    if hsrp > h[-1]:
        hsrd = h[-1]
    else:
        hsrd = hsrp

    #  Parameters for the ducting/layer-reflection model

    hst = min(hst, h[0])
    hsr = min(hsr, h[-1])
    
    #  Slope of the smooth-Earth surface eq 169
    m = (hsr - hst) / dtot
    
    #  The terminal effective heigts for the ducting/layer-reflection model eq170
    
    hte = htg + h[0] - hst
    hre = hrg + h[-1] - hsr

    if (h[ltindex:lrindex+1] - (hst + m * d[ltindex:lrindex+1])).any():
        hm = np.max(h[ltindex:lrindex+1] - (hst + m * d[ltindex:lrindex+1]))
    else:
        hm = 0


    return hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype
