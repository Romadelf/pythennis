import numpy as np
from scipy.integrate import solve_ivp

tol = 10**(-5) #tolerance 
coef = 0.7 # coef de changement de vz
distance_maximal_terrain = 11.89

def oderhs(t, instant_ball_data):
    
    # definition des parametre du systeme
    d = 0.065
    m = 0.058
    rho = 1.2
    Cd = 0.65
    g = 9.81
    
    v = np.array([instant_ball_data[3],
                  instant_ball_data[4],
                  instant_ball_data[5]])
    r = np.array([instant_ball_data[6],
                  instant_ball_data[7],
                  instant_ball_data[8]])
    norm_r = np.linalg.norm(r)
    norm_v = np.linalg.norm(v)
    pi = np.pi
    Cm = 1 / (2 + 1.96 * norm_v / (norm_r * d))
    
    Fg = np.array([0, 0, -m * g])
    Fd = Cd * rho * pi * d**2 * -v**2 / 8
    ev = np.cross(r, v)
    n = ev / np.linalg.norm(ev)
    Fm = Cm * rho * pi * d**2 * np.linalg.norm(v)**2 / 8 * n 
    
    a = (Fm + Fd + Fg)/m
    dy = np.array([instant_ball_data[3],
                   instant_ball_data[4],
                   instant_ball_data[5],
                   a[0],
                   a[1],
                   a[2],
                   0,
                   0,
                   0])
    
    return dy

def bouing(t, instant_ball_data): # event pour verifier quand x3 vaut zero
    return instant_ball_data[2]
bouing.direction = -1
bouing.terminal = True

def trajectoireFiletHorizontal(initial_ball_data, t_f):
    
    # définition des paramètres :
        
    frequence = 10**(3)
    instants_a_evaluer = np.linspace(0, t_f, int(frequence * t_f))
    
    # verification des donnés :
    # si vitesse initial de translation nulle, il n'y a pas de mouvements
    # comme cela pose des problèmes de division par 0 dans oderhs,
    # on exclut ce cas, lequel est de toute façon anecdotique
    # => à voir si on garde
    if(abs(initial_ball_data[3]) <= tol and
       abs(initial_ball_data[4]) <= tol and
       abs(initial_ball_data[5]) <= tol):
        return [0,0,0]
    
    pre_bounce_solve = solve_ivp(
        oderhs,
        [0, t_f],
        initial_ball_data,
        t_eval = instants_a_evaluer,
        events = bouing)

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
            events = bouing)
           
        ball_data_timetable = np.concatenate((ball_data_timetable, post_bounce_solve.y), axis=1) # on regroupe les tableaux

    # debut :
    
    x=ball_data_timetable[0]  #longeur
    y=ball_data_timetable[1]  #largeur
    z=ball_data_timetable[2]  #hauteur
    
    # passe le filet ou pas ?
    # => on va parcourir les valeurs jusqu'au filet et sans dépasser l'indice maximal
    
    index_max = np.shape(x)[0] - 1 # car shape(x) = shape(y) = shape(z)
    i = 0
    while i <= index_max and x[i] <= tol: # tant que  on est avant le filet 
        if  z[i] <= hauteur_filet(y[i]) :  # et si la hauteur de la balle est inferieur a celle du filet 
                return [0,0,0] # hors-jeu 
        i += 1
        
    return [
        x[index_max],
        y[index_max],
        z[index_max]
    ] # on aurait aussi pu utiliser -1 (e.g.: x[-1])

def hauteur_filet(dist_centre):
    return 1 # temporaire, pour le milestone 2
#    bord = 8.23 / 2 + 0.914 # = 5.029
#    h_nominale = 1.07
#    delta_h = h_nominale - 0.914
#    
#    if dist_centre >= bord or dist_centre <= -bord :
#        return h_nominale
#    else :
#        return h_nominale - np.cos(dist_centre * (np.pi / 2) / bord) * delta_h
