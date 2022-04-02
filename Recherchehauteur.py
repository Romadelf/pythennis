 
import Trajectoire as Tj
import RechercheRacine as R
import cible as ci    

 
def rechercheHauteur(y0,cibleRebond): 
    #argument de fonction : 
    def f_hauteur(h):
        #Definition d'event
      

    
        y0[2] =  h 
        solve = Tj.solve_ivp(
            Tj.oderhs,
            [ci.t_i, ci.t_f],
            y0,
            t_eval = ci.instants_a_evaluer,
            events =ci.event_hauteur ,
            rtol = 10**(-5),
            atol = 10**(-8))
    
        timetable = solve.y
        x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
        return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
    return R.bissection(f_hauteur, ci.haut_min, ci.haut_max , ci.tol) #bissection car se comporte mieux 


def rechercheHauteur2(initial_ball_data,cibleHauteur):

      
    def hauteur2(h):
        # définition des paramètres :
        initial_ball_data[2] = h
        frequence = 10**(5)
        instants_a_evaluer = Tj.np.linspace(0, ci.t_f, int(ci.frequence * ci.t_f))
        

        
        pre_bounce_solve = Tj.solve_ivp(
            Tj.oderhs,
            [0, ci.t_f],
            initial_ball_data,
            t_eval = instants_a_evaluer,
            events = Tj.bouing,
            rtol = 10**(-5),
            atol = 10**(-8))

        ball_data_timetable = pre_bounce_solve.y
        
        x_rebond = ball_data_timetable[0,-1]
        #print ('la balle rebondit a ', x_rebond)
        
        if pre_bounce_solve.status == 1 : # = si il rebondi
            # manière d'aller chercher les dernières variables dans ball_data_timetable
            # TODO_LOW possiblement plus clair d'utilisation via pre_bounce_solve.t.shape[0]
            after_bounce_ball_data = ball_data_timetable[0:ball_data_timetable.shape[0], ball_data_timetable.shape[1] - 1]
            after_bounce_ball_data[5] *= Tj.coef
            
            # reinitialisation des données :
            # frequence * (delta_t) == on cree un tableau de frequence fixe allant du nouveau t_i au t_f donné 
            # ce qui assure une precision a chaque rebond quel que soit t_f donné
            
            t_i = pre_bounce_solve.t[-1] # avec -1, on lit le tableau à l'envers, donc on trouve le dernier élément
            delta_t = ci.t_f - t_i
            instants_a_evaluer = Tj.np.linspace(t_i, ci.t_f, int(frequence * delta_t))
            
            # on relance solve_ivp avec les nouvelles cond 
            
            post_bounce_solve =Tj.solve_ivp(
                Tj.oderhs,
                [t_i, ci.t_f],
                after_bounce_ball_data,
                t_eval = instants_a_evaluer,
                events = [ci.stop],
                rtol = 10**(-5),
                atol = 10**(-8))
               
            ball_data_timetable = Tj.np.concatenate((ball_data_timetable, post_bounce_solve.y), axis=1) # on regroupe les tableaux


        z=ball_data_timetable[2]  #hauteur
            


        return  z[-1]-cibleHauteur
    return R.bissection( hauteur2,ci.haut_min,ci.haut_max, ci.tol)

