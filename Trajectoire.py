import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt

coef = 0.7 # coef de changement de vz
distance_maximal_terrain = 11.89

def oderhs(t, instant_ball_data):
    
    # definition des parametre du systeme
    d = 0.065
    m = 0.058
    p = 1.2
    C_d = 0.65
    g = 9.81
    
    # ici on initialise notre repere orthonormé
    vec_u_x1 = [1, 0, 0] # x
    vec_u_x2 = [0, 1, 0] # y
    vec_u_x3 = [0, 0, 1] # z
   
    # puis on calcule les normes
    vitesse_translation = [instant_ball_data[3],
                           instant_ball_data[4],
                           instant_ball_data[5]]
    vitesse_angulaire = [instant_ball_data[6],
                         instant_ball_data[7],
                         instant_ball_data[8]]
    norme_de_v = np.linalg.norm(vitesse_translation)
    norme_de_w = np.linalg.norm(vitesse_angulaire)
   
    # norme des forces
    C_m = 1 / (2 + 1.96 * norme_de_v / (norme_de_w * d))
    force_frottement = C_d * p * np.pi * d**2 / 8 * norme_de_v**2    
    force_magnus     = C_m * p * np.pi * d**2 / 8 * norme_de_v**2
    
    # ici on a un vecteur unitaire de la vitesse de translation
    direction_vitesse_translation = vitesse_translation / norme_de_v
    
    # mainenant on determine les composantes grace au produit scalaire des deux 
    # vecteur unitaire  pour obtenir nos trois composante
    # on fait (-) la direction car c'est l'inverse de la vitesse
    fdx1 = force_frottement * np.dot(-direction_vitesse_translation, vec_u_x1)
    fdx2 = force_frottement * np.dot(-direction_vitesse_translation, vec_u_x2)
    fdx3 = force_frottement * np.dot(-direction_vitesse_translation, vec_u_x3)
    
    # fmx1? fmx2? fmx3? en fonction de leur direction
    # ici on fait le produit vectoriel v * w puis on le norme
    produit_vectoriel = np.cross(vitesse_angulaire, vitesse_translation)
    norme_produit_vectoriel = np.linalg.norm(produit_vectoriel)
    direction_force_magnus = produit_vectoriel / norme_produit_vectoriel
    # Autre manière de le faire :
    # direction_force_magnus =  np.cross(
    #     vitesse_angulaire / np.linalg.norm(vitesse_angulaire),
    #     vitesse_translation / np.linalg.norm(vitesse_translation))
    
    fmx1 = force_magnus * np.dot(direction_force_magnus, vec_u_x1)
    fmx2 = force_magnus * np.dot(direction_force_magnus, vec_u_x2)
    fmx3 = force_magnus * np.dot(direction_force_magnus, vec_u_x3)
    
    # Poids
    fgx1 = 0
    fgx2 = 0
    fgx3 = -m * g

    # definition du systeme
    # on pose :
    # X = dx1_dt   ainsi on a  dX_dt = acceleration selon x1
    # Y = dx2_dt               dY_dt = acceleration selon x2
    # Z = dx3_dt               dZ_dt = acceleration selon x3
    
    # vecteur instant_ball_data (x1,x2,x3,dx1_dt,dx2_dt,dx3_dt,w1,w2,w3)
    # indice                     0  1  2  3      4       5     6  7  8
    
    dy = np.zeros(9)
    
    # on met les trois composantes de la vitesse, 
    dy[0] = instant_ball_data[3]
    dy[1] = instant_ball_data[4]
    dy[2] = instant_ball_data[5]
    
    # les trois composantes de l'acceleration
    dy[3] = (fgx1 + fdx1 + fmx1) / m
    dy[4] = (fgx2 + fdx2 + fmx2) / m
    dy[5] = (fgx3 + fdx3 + fmx3) / m
    
    # pas necessaire de remvoyer omega car tableau initialisé a zero
    return dy

def bouing(t, instant_ball_data): # event pour verifier quand x3 vaut zero
    return instant_ball_data[2]
bouing.direction = -1
bouing.terminal = True

def trajectoireFiletHorizontal(initial_ball_data, t_f):
    
    # définition des paramètres :
    
    frequence = 10**(2)
    instants_a_evaluer = np.linspace(0, t_f, int(frequence * t_f))
    
    # verification des donnés :
    # si vitesse initial de translation nulle, il n'y a pas de mouvements
    # comme cela pose des problèmes de division par 0 dans oderhs,
    # on exclut ce cas, lequel est de toute façon anecdotique
    # => à voir si on garde
    if(abs(initial_ball_data[3]) <= 0 and
       abs(initial_ball_data[4]) <= 0 and
       abs(initial_ball_data[5]) <= 0):
        return [0,0,0]
    
    pre_bounce_solve = solve_ivp(
        oderhs,
        [0, t_f],
        initial_ball_data,
        t_eval = instants_a_evaluer,
        events = bouing,
        rtol = 10**(-5),
        atol = 10**(-8))

    ball_data_timetable = pre_bounce_solve.y
    
    if pre_bounce_solve.status == 1 : # = si il rebondi
        # manière d'aller chercher les dernières variables dans ball_data_timetable
        # TODO_LOW possiblement plus clair d'utilisation via pre_bounce_solve.t.shape[0]
        after_bounce_ball_data = ball_data_timetable[0:ball_data_timetable.shape[0], ball_data_timetable.shape[1] - 1]
        after_bounce_ball_data[5] *= -coef
        
        # reinitialisation des données :
        # frequence * (delta_t) == on cree un tableau de frequence fixe allant du nouveau t_i au t_f donné 
        # ce qui assure une precision a chaque rebond quel que soit t_f donné
        
        t_i = pre_bounce_solve.t[-1] # avec -1, on lit le tableau à l'envers, donc on trouve le dernier élément
        delta_t = t_f - t_i
        instants_a_evaluer = np.linspace(t_i, t_f, int(frequence * delta_t))
        
        # on relance solve_ivp avec les nouvelles cond 
        
        post_bounce_solve = solve_ivp(
            oderhs,
            [t_i, t_f],
            after_bounce_ball_data,
            t_eval = instants_a_evaluer,
            events = bouing,
            rtol = 10**(-5),
            atol = 10**(-8))
           
        ball_data_timetable = np.concatenate((ball_data_timetable, post_bounce_solve.y), axis=1) # on regroupe les tableaux

    # debut :
    
    x=ball_data_timetable[0]  #longeur
    y=ball_data_timetable[1]  #largeur
    z=ball_data_timetable[2]  #hauteur
    
    plt.plot(x,z)
    plt.title('filet') #Ajout d'un titre
    plt.ylabel('Axe des Y') # Labélisation de l'axe des ordonnées
    plt.xlabel('Axe des X')
    
    # passe le filet ou pas ?
    # => on va parcourir les valeurs jusqu'au filet et sans dépasser l'indice maximal
    
    index_max = np.shape(x)[0] - 1 # car shape(x) = shape(y) = shape(z)
    i = 0
    while i <= index_max and x[i] <= 0: # tant que  on est avant le filet 
        if  z[i] <= hauteur_filet(y[i]) :  # et si la hauteur de la balle est inferieur a celle du filet 
                return [0,0,0] # hors-jeu 
        i += 1
        
    return [
        x[index_max],
        y[index_max],
        z[index_max]
    ] # on aurait aussi pu utiliser -1 (e.g.: x[-1])

def hauteur_filet(dist_centre):
    bord = 8.23 / 2 + 0.914 # = 5.029
    h_nominale = 1.07
    delta_h = h_nominale - 0.914
    
    if dist_centre >= bord or dist_centre <= -bord :
        return h_nominale
    else :
        return h_nominale - np.cos(dist_centre * (np.pi / 2) / bord) * delta_h