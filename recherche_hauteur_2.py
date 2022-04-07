from scipy.integrate import solve_ivp

from RechercheRacine import bissection
from Trajectoire import oderhs

marge_ratio = 10
tol = 0.001 # 1 mm
ligne_fond = 11.89
coef_restitution = 0.7
t_max = 5  #TODO essayer sys.maxsize ?

"""
#FIXME
Défauts connus :
    - Ne prends pas en compte le filet
    - Fonction non réplicable (du genre rate la premiere fois puis réussi à la suivante)
    - Plusieurs incohérences dans le code marqués par des TODO doivent être résolu
    
    Ceci dit, la fonction à dejà fait ses preuve malgré, bien qu'irrégulièrement
"""
def rechercheHauteur2(y0, cibleHauteur):
    """Renvoie la hauteur initiale nécessaire d'une balle dont les autres composantes sont fixées pour qu'elle atteigne une hauteur cible donnée au niveau de la ligne de fond, après rebond.
#TODO documentation à finir
    Paramètres
    ----------
    y0 : vecteur des composantes # preciser les y0[0 à 8]
    op : callable
        This should be either numpy.amin or `numpy.amax` or `numpy.sum`.

    Returns
    -------
    result : float or ndarray
        If `x` is 2-D, the return values is a float.
        Otherwise, it is an array with ``x.ndim - 2`` dimensions.
        The return values are either the minimum or maximum or sum of the
        singular values of the matrices, depending on whether `op`
        is `numpy.amin` or `numpy.amax` or `numpy.sum`.
    """
    def pos_x0_a_cibleHauteur_apres_rebond(h_init):
        #TODO: Ne prends pas en compte la hauteur après rebond lorsque la balle retombe vu les complications au niveau de la continuité => eventuellement voir si c'est possible
        y0[2] = h_init
        h_cible_events = solve_ivp(oderhs,
                                   [0, t_max],
                                   y0,
                                   events = [pseudo_ev_rebond, ev_h_cible],
                                   max_step = 0.01
                                   ).y_events[1]
        if(h_cible_events.size > 0):
            return h_cible_events[0][0] - ligne_fond #TODO attention: là on vérifie pas si il y a un rebond
        else:
            return 42 #TODO: à améliorer mais en gros la balle n'atteint pas la hauteur cible demandée
    def pseudo_ev_rebond(t, y):
        if(y[2] <= 0 and y[5] < 0):
            y[5] *= -coef_restitution
            ev_h_cible.terminal = True
        return y[2]
    pseudo_ev_rebond.direction = -1
    pseudo_ev_rebond.terminal = False
    
    def ev_h_cible(t, y):
        if(y[5] < 0):#TODO: patch temporaire à retirer à terme
            return -1
        return y[2] - cibleHauteur
    ev_h_cible.direction = +1
    ev_h_cible.terminal = False
    
    return bissection(pos_x0_a_cibleHauteur_apres_rebond,
                  y0[2] / marge_ratio,
                  y0[2] * marge_ratio,
                  tol)