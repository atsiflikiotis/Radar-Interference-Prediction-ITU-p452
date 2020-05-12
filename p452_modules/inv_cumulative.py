import numpy as np


def inv_cumulative(x):
    if x < 0.000001:
        x = 0.000001
            
    # if x > 0.5: #warning 
    
    T = np.sqrt(-2*np.log(x))
    
    C0 = 2.515516698
    C1 = 0.802853
    C2 = 0.010328
    D1 = 1.432788
    D2 = 0.189269
    D3 = 0.001308
    
    k = ( (C2 * T + C1) * T + C0 )/ ( ((D3*T + D2)*T + D1)*T + 1)
    
    I = k - T

    return I
