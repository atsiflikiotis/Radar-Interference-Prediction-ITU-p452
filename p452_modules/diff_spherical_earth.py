import numpy as np


def sph_earth_loss(dtot, hte, hre, ap, freq, omega, pol):
    # The spherical-Earth diffraction loss not exceeded for p% time for antenna
    # heights hte and hre, equations 23-28

    lam = 0.3/freq
    
    # marginal LoS distance for a smooth path
    dlos = np.sqrt(2*ap) * (np.sqrt(0.001*hte) + np.sqrt(0.001*hre))
    
    if dtot >= dlos:
        # calculate diffraction loss Ldft using the method in Sec. 4.2.2.1 for 
        # adft = ap to obtain Ldft

        # Calculate Ldftland using eqs 30 to 37
        er = 22
        sigma = 0.003
        Ldftland = se_first_term(dtot, hte, hre, ap, freq, er, sigma, pol)

        # Calculate Ldftland using eqs 30 to 37
        er = 80
        sigma = 5
        Ldftsea = se_first_term(dtot, hte, hre, ap, freq, er, sigma, pol)

        # Calculation of ldft eq (29):
        Ldsph = omega * Ldftsea + (1 - omega) * Ldftland
    else:
        # calculate the smallest clearance between the curved-Earth path and
        # the ray between the antennas, hse
        
        # eq 25a-25e:
        c = (hte - hre)/(hte + hre)  
        m = 250*dtot**2 / (ap * (hte + hre))
        b = 2*np.sqrt((m+1)/(3*m)) * np.cos(np.pi/3 + 1/3 * np.arccos(3*c/2 * np.sqrt(3*m/(m+1)**3)))
        dse1 = dtot / 2*(1+b)
        dse2 = dtot - dse1
        
        # eq 24:
        hse = ((hte - 500*dse1*dse1 / ap) * dse2 + (hre - 500*dse2*dse2 / ap) * dse1)/dtot

        # required clearance for zero diffraction loss, eq 26
        hreq = 17.456*np.sqrt(dse1 * dse2 * lam/dtot)
        
        if hse > hreq:
            Ldsph = 0
        else:
            # calculate the modified effective Earth radius aem, which gives
            # marginal LoS at distance dtot (eq 27):
            
            aem = 500*(dtot/(np.sqrt(hte) + np.sqrt(hre)))**2

            ## Use the method in Sec. 4.2.2.1 for adft = aem to obtain Ldft
            
            # Calculate Ldftland using eqs 30 to 37
            er = 22
            sigma = 0.003
            Ldftland = se_first_term(dtot, hte, hre, aem, freq, er, sigma, pol)
            
            # Calculate Ldftland using eqs 30 to 37
            er = 80
            sigma = 5
            Ldftsea = se_first_term(dtot, hte, hre, aem, freq, er, sigma, pol)
            
            # Calculation of ldft eq (29):
            Ldft = omega * Ldftsea + (1-omega)*Ldftland
            
            if Ldft < 0:
                Ldsph = 0
            else:
                # eq 28
                Ldsph = (1 - hse/hreq)*Ldft
    return Ldsph


def se_first_term(dtot, hte, hre, aem, freq, er, sigma, pol):
    # 4.2.2.1 First-term part of spherical-Earth diffraction loss
    # This sub-section gives the method for calculating spherical-Earth diffraction using only the first
    # term of the residue series. It forms part of the overall diffraction method described in § 4.2.2 above
    # to give the first-term diffraction loss, L dft , for a given value of effective Earth radius a dft . The value
    # of a dft to use is given in § 4.2.2.

    # Normalized factorfor surface admittance, eq 30
    # horizontal polarization
    Kh = 0.036 * (aem * freq)**(-1 / 3) * ((er - 1)**2 + (18 * sigma / freq)**2) ** (-1 / 4)

    if pol == 1:
        K = Kh
    elif pol == 2:
        # horizontal polarization
        K = Kh * (er**2 + (18 * sigma / freq)**2) ** 0.5
    else:
        # slant polatization (+/- 45)
        Kv = Kh * (er**2 + (18 * sigma / freq)**2) ** 0.5
        K = np.sqrt(Kv**2 + Kh**2)

    # Earth ground / polarization parameter, eq 31
    beta_dft = (1 + 1.6 * K**2 + 0.67 * K**4) / (1 + 4.5 * K**2 + 1.53 * K**4)

    # Normalized distance eq 32
    X = 21.88 * beta_dft * (freq/aem**2)**(1 / 3) * dtot

    # Normalized transmitter and receiver heights eq 33a,33b
    Yt = 0.9575 * beta_dft * (freq**2 / aem) ** (1 / 3) * hte
    Yr = 0.9575 * beta_dft * (freq**2 / aem) ** (1 / 3) * hre

    # Calculate the distance term given by eq 34:

    if X >= 1.6:
        Fx = 11 + 10 * np.log10(X) - 17.6 * X
    else:
        Fx = -20 * np.log10(X) - 5.6488 * (X ** 1.425)

    # eq 36:
    Bt = beta_dft * Yt
    Br = beta_dft * Yr

    # Define a function of normalized height given by eq 35
    if Bt > 2:
        GYt = 17.6 * (Bt - 1.1)**0.5 - 5 * np.log10(Bt - 1.1) - 8
    else:
        GYt = 20 * np.log10(Bt + 0.1 * Bt ** 3)

    if Br > 2:
        GYr = 17.6 * (Br - 1.1) ** 0.5 - 5 * np.log10(Br - 1.1) - 8
    else:
        GYr = 20 * np.log10(Br + 0.1 * Br ** 3)

    if GYr < 2 + 20 * np.log10(K):
        GYr = 2 + 20 * np.log10(K)

    if GYt < 2 + 20 * np.log10(K):
        GYt = 2 + 20 * np.log10(K)

    # The first-term spherical-Earth diffraction loss is now given by eq 37:
    Ldft = -Fx - GYt - GYr

    return Ldft





