import numpy as np

def clutterloss(f, d, h, htg, hrg, ha_t, ha_r, dk_t, dk_r):
    
    idx1 = 0
    idx2 = len(d)
    
    htgcorr = htg
    hrgcorr = hrg
    
    Aht = 0
    Ahr = 0
        
    if ha_t > htg:
        
        Ffct = 0.25 + 0.375 * (1 + np.tanh(7.5 * (f - 0.5)))  # eq.57a
        
        Aht = 10.25 * Ffct * np.exp(-dk_t) * (1 - np.tanh(6 * (htg / ha_t - 0.625))) - 0.33 # eq. 57
        
        idxs = np.argwhere(d >= dk_t)
        
        if idxs.size > 0:
            idx1 = idxs[0].item()
        else:
            idx1 = len(d)
        
        htgcorr = ha_t
    
    
    if ha_r > hrg:
        
        Ffcr = 0.25 + 0.375 * (1 + np.tanh(7.5 * (f - 0.5)))
        
        Ahr = 10.25 * Ffcr * np.exp(-dk_r) * (1 - np.tanh(6 * (hrg / ha_r - 0.625))) - 0.33  
        
        idxs = np.argwhere(d <= d[-1] - dk_r)
       
        if idxs.size > 0:
            idx2 = idxs[-1].item()
        else:
            idx2 = 0
        
        hrgcorr = ha_r

        
    # Adjust path length
 
    dcorr = d[idx1:idx2+1] - d[idx1]
    
    hcorr = h[idx1:idx2+1]
    
    
    return dcorr, hcorr, htgcorr, hrgcorr, Aht, Ahr

