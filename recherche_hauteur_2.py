import copy

from scipy.integrate import solve_ivp

from RechercheRacine import bissection
from Trajectoire import oderhs
from Trajectoire import hauteur_filet

marge_ratio = 1.25
tol = 0.001 # 1 mm
ligne_fond = 11.89
coef_restitution = 0.7
t_max = 10  #TODO essayer sys.maxsize ? # Ce temps doit être suffisant pour garantir que la balle atteigne la ligne de fond ou rebondisse deux fois

"""
#FIXME La logique générale est vérifiée sauf la notion des bornes pour bissection.
Voir le FIXME au niveau du return bissection(...)
    
    Ceci dit, la fonction fait déja ses preuves dans certaines conditions initiales favorables
"""
def rechercheHauteur2(y0, cibleHauteur):
    """Renvoie la hauteur initiale nécessaire d'une balle dont les autres composantes sont fixées pour qu'elle atteigne une hauteur cible donnée au niveau de la ligne de fond, après rebond.
#TODO documentation à finir (la suite vient d'un copier-collé pour la structure)
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

    y0 = copy.deepcopy(y0) # Assure que le vecteur initial passé en argument ne sera pas muté
    n_bounces_0_ = [0]
    
    def h_fond_apres_rebond(h_init):
        y0[2] = h_init
        n_bounces_0_[0] = 0
        events = solve_ivp(oderhs,
                           [0, t_max],
                           y0,
                           events = [ev_filet,
                                     ev_rebonds,
                                     ev_ligne_fond],
                           max_step = 0.01
                           ).y_events

        if(len(events[2]) == 1 and n_bounces_0_[0] == 1):
            return events[2][0][2] - cibleHauteur
        else:
            return -cibleHauteur

    def ev_filet(t, y):
        if(y[0] <= 0):
            return y[2] - hauteur_filet(y[1])
        return 1
    ev_filet.terminal = True
    ev_filet.direction = -1
    
    def ev_rebonds(t, y):
        if(y[2] <= 0 and y[5] < 0):
            y[5] *= -coef_restitution
            n_bounces_0_[0] += 1
        return n_bounces_0_[0] - 2
    ev_rebonds.terminal = True
    
    def ev_ligne_fond(t, y):
        return y[0] - ligne_fond
    ev_ligne_fond.terminal = True
    
    return bissection(h_fond_apres_rebond,
                  y0[2] / marge_ratio, #FIXME erreur de logique : rien ne garantis que les bornes initiales donneront des images de signes opposées
                  y0[2] * marge_ratio, # idem
                  tol)