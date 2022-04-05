#E = (1/2)*m*v**2  + (1/8)*m*(d**2)*(w**2)

 
import Trajectoire as Tj
import RechercheRacine as R
import RechercheVitesse as RV
import RechercheOmega as RW
import cible as ci    
 

def methode (E , vec):
    
    #condition initial 
    x = [-11.89 ,0 ,  2]
    frequence = 10**(5)
    instants_a_evaluer = Tj.np.linspace(0, ci.t_f, int(ci.frequence * ci.t_f))
    
   # simpl = 2*E/(Tj.m) #version simplifieer
    #ainsi on a  : simpl = v**2 + 0.25*d*(w**2)
    
    
        
    #parametre 
        
    H_max , H_inter = 0 ,0 
    norm_v, norm_w = 0 , 0 
    coef= 0.1 #m
    
    #debut :
    #idee on raisonne en terme de pourcentage
    #pour chaque pourcentage on verifie la hauteur atteint apres un rebond
    #qu'on stoke dans une variable Hmax par l'intermediaire de Hinter 
    
    while coef<1 :
       norm_v = Tj.np.sqrt(2*coef*E/Tj.m)                 #norme
       norm_w = Tj.np.sqrt((1-coef)*8*E/(Tj.m*(Tj.d**2)))   #norme
       v = norm_v*vec #vecteur
       w = norm_w*vec #vecteur
    
       y = [x , v , w]
    
      
       pre_bounce_solve = Tj.solve_ivp(
            Tj.oderhs,
            [0, ci.t_f],
            y,
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
    
       H_inter = z[-1]


       if H_inter > H_max :
           H_max = H_inter
           bon_coef = coef 
    
   #si il rentre dans la condition if alors cela veut dire que le coef 
   #present a une hauteur plus elevé 
   #a la fin de la boucle soit 9iteration boncoef aura stoque le dernier 
   #coef pour lequel Hinter est egal Hmax
       coef+= 0.1
    
    
    
    return bon_coef 
    