# 4.2.4 The diffraction loss not exceeded for p% time

from p452_modules import inv_cumulative, delta_bullington


def diffloss(d, h, hts, hrs, hstd, hsrd, f, omega, p, b0, ae, ab, pol):
    # I give pol as parameter to extract only 1 value depending on polarization
    # a = ae
    Ld50 = delta_bullington.delta_bull(d, h, hts, hrs, hstd, hsrd, ae, f, omega, pol)


    if p == 50:
        Ldp = Ld50
    else:
        # 1st step Use the method in § 4.2.3 to calculate diffraction loss Ld 
        # for effective Earth radius a as given in
        # equation (6b). Set diffraction loss not exceeded for β0% time Ldβ=Ld.
        # a = ab
        Ldb = delta_bullington.delta_bull(d, h, hts, hrs, hstd, hsrd, ab, f, omega, pol)
        
        # interpolation factor Fi (eq 41):
        if p > b0:
            Fi = inv_cumulative.inv_cumulative(p / 100) / inv_cumulative.inv_cumulative(b0 / 100)
        else:
            Fi = 1
        
        # The diffraction loss not exceeded for p% time eq 42
        Ldp = Ld50 + Fi*(Ldb - Ld50)

    return Ldp, Ld50