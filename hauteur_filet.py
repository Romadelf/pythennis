import numpy as np

def hauteur_filet(dist_centre):
    
    bord = 8.23 / 2 + 0.914 # = 5.029
    h_nominale = 1.07
    delta_h = h_nominale - 0.914
    
    if dist_centre >= bord or dist_centre <= -bord :
        return h_nominale
    else :
        return h_nominale - np.cos(dist_centre * (np.pi / 2) / bord) * delta_h