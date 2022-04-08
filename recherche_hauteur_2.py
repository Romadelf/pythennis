from scipy.integrate import solve_ivp

from RechercheRacine import bissection
from Trajectoire import oderhs

marge_ratio = 10
tol = 0.001 # 1 mm
ligne_fond = 11.89
coef_restitution = 0.7
t_max = 10  #TODO essayer sys.maxsize ? # Ce temps doit être suffisant pour garantir que la balle atteigne la ligne de fond ou rebondisse deux fois

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

    n_bounces_0_ = [0]
    
    def h_fond_apres_rebond(h_init):
        y0[2] = h_init
        n_bounces_0_[0] = 0
        y_fond_0_ = solve_ivp(oderhs,
                                   [0, t_max],
                                   y0,
                                   events = ev_ligne_fond,
                                   max_step = 0.01
                                   ).y_events[0]
        if(y_fond_0_.size > 0):
            if(n_bounces_0_[0] == 1):
                return y_fond_0_[0][2] - cibleHauteur
            else:
                return -cibleHauteur

    def ev_ligne_fond(t, y):

        # Rebond
        if(y[2] <= 0 and y[5] < 0):
            y[5] *= -coef_restitution
            n_bounces_0_[0] += 1

        return y[0] - ligne_fond
    ev_ligne_fond.terminal = True
    
    return bissection(h_fond_apres_rebond,
                  y0[2] / marge_ratio,
                  y0[2] * marge_ratio,
                  tol)