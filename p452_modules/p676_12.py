# Attenuation by atmospheric gases, ITU-R P.676-12 recommendation

import numpy as np

def specatten(f, pressure, rho, T):
    oxygen = np.genfromtxt('p676-12_table_a.csv', delimiter=";")
    vapour = np.genfromtxt('p676-12_table_b.csv', delimiter=";") 
    
    f0_oxygen = oxygen[:, 0]
    a1 = oxygen[:, 1]
    a2 = oxygen[:, 2]
    a3 = oxygen[:, 3]
    a4 = oxygen[:, 4]
    a5 = oxygen[:, 5]
    a6 = oxygen[:, 6]
    
    f0_vapour = vapour[:, 0]
    b1 = vapour[:, 1]
    b2 = vapour[:, 2]
    b3 = vapour[:, 3]
    b4 = vapour[:, 4]
    b5 = vapour[:, 5]
    b6 = vapour[:, 6]
    
    theta = 300/T
    e = rho * T / 216.7

    ## Oxygen calculations (eqs 5-8 in P.676-12 recommendation)
    
    Si_o = a1 * 1e-7 * pressure * theta**3 * np.exp(a2 * (1 - theta))
    
    df_o = a3 * 1e-4 * (pressure * theta ** (0.8 - a4) + 1.1 * e * theta)
    
    # Zeeman splitting of oxygen lines and Doppler broadening
    
    df_o = np.sqrt(df_o**2 + 2.25e-6)
    
    delta = (a5 + a6 * theta) * 1e-4 * (pressure + e) * theta**0.8

    # equation 5:
    Fi_oxygen = (f / f0_oxygen) * ((df_o - delta * (f0_oxygen - f)) /
                                   ((f0_oxygen - f)**2 + df_o**2) + (df_o - delta * (f0_oxygen + f))
                                   / ((f0_oxygen + f)**2 + df_o**2))
    
    d0 = 5.6e-4 * (pressure + e) * theta**0.8
    
    Ndf = f * pressure * theta**2 * (6.14e-5/(d0 * (1 + (f/d0)**2)) +
                                     1.4e-12 * pressure * theta**1.5 / (1 + 1.9e-5 * f**1.5))

    # vapour:
    Si_v = b1 * 1e-1 * e * theta**3.5 * np.exp(b2 * (1 - theta))
    df = b3 * 1e-4 * (pressure * theta ** b4 + b5 * e * theta**b6)
    df = 0.535 * df + np.sqrt(0.217 * df**2 + 2.1316e-12 * f0_vapour * f0_vapour/theta)
    delta = 0

    Fi_vapour = (f / f0_vapour) * ((df - delta * (f0_vapour - f)) /
                                   ((f0_vapour - f)**2 + df**2) + (df - delta * (f0_vapour + f))
                                   / ((f0_vapour + f)**2 + df**2))
                     
    # specific attenuations (dB/km) due to dry air (oxygen, pressure-induced 
    # nitrogen and non-resonant Debye attenuation) and water vapour (eq 1-2)
    g_0 = 0.182 * f * (sum(Si_o * Fi_oxygen) + Ndf)
    g_w = 0.182 * f * sum(Si_v * Fi_vapour)
    
    return g_0, g_w
