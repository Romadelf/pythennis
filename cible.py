#angle, statut = rechercheAngle(y0,cibleRebond)
#angle, statut = rechercheAngle2(y0,cibleHauteur)
#normeVitesse, statut = rechercheVitesse(y0,cibleRebond)
#normeVitesse, statut = rechercheVitesse2(y0,cibleHauteur)
#hauteur, statut = rechercheHauteur(y0,cibleRebond)
#hauteur, statut = rechercheHauteur2(y0,cibleHauteur)
#omega, statut = rechercheOmega(y0,cibleRebond)
#omega, statut = rechercheOmega2(y0,cibleHauteur)

import Trajectoire as Tj
import RechercheRacine as R

#Condition initial :     d   

y0 = [-11.89,0,2,50,1,0,30,15,0]    
cond =  [-11.89,0,2,50,1,0,30,15,0]    
#parametre : 
    
t_i , t_f = 0, 4  #t_span 
n_point = 200  #nombre de point 
n_point_apres = 150  #apres rebond
instants_a_evaluer = Tj.np.linspace(t_i, t_f,n_point ) 
  
statut  , error = 0 , 0
tol =  1/1000 #tolerance au millimetre 


#intervalle 

haut_min,haut_max =  1 , 10 # m
theta_min, theta_max = 0 , (Tj.np.pi)/2
vitesse_min, vitesse_max = 1 , 500
omega_min, omega_max = 1 , 500

cibleHauteur = 1 #arbitraire



#utilise dans toutes les  fonctions cibleRebond donc general
def event_hauteur (t,y0):
      return y0[2]
event_hauteur.direction = -1
event_hauteur.terminal = True



#utilise dans toutes les  fonctions cibleHauteur donc general

       



#convertie coordonne cartesienne en cylindrique:
  
def convert_Car_to_Cy(x , y , z):
    
    if x == 0:
       return [0,0,0] 
    r = Tj.np.sqrt((x**2)+(y**2))
    theta = Tj.np.arctan(y/x)
    z = z 
    return [r , theta , z]

#convertie coordonne  cylindrique  en cartesienne :
  
def convert_Cy_to_Car(r , theta , z ):
    
    x = r * Tj.np.cos(theta)
    y = r * Tj.np.sin(theta)
    z = z
       
    return  [x , y , z ]

# normalise l'angle pour qu'il reste entre [0 , pi/2]
def normalise(angle):
    
  #  if angle>=0 and angle<= Tj.np.pi:
  #      angle_normalise = angle 
        
    if angle < 0:
        
        while(angle < 0):
           angle +=Tj.np.pi 
           
    elif  angle > Tj.np.pi  :
        
        while(angle > Tj.np.pi):
           angle +=Tj.np.pi 
    angle_normalise =angle 
    
    return angle_normalise 





      
def rechercheHauteur(y0,cibleRebond): 
    #argument de fonction : 
    def f_hauteur(h):
        #Definition d'event
      

    
        y0[2] =  h 
        solve = Tj.solve_ivp(
            Tj.oderhs,
            [t_i, t_f],
            y0,
            t_eval = instants_a_evaluer,
            events =event_hauteur ,
            rtol = 10**(-5),
            atol = 10**(-8))
    
        timetable = solve.y
        x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
        return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
    return R.bissection(f_hauteur, haut_min, haut_max , tol) #bissection car se comporte mieux 



        
def rechercheAngle(y0,cibleRebond) :
   
#on convertie d'abord pour changer theta sans change la norme  de la VITESSE 
   pos_cylin = convert_Car_to_Cy(y0[3],y0[4],y0[5])
  
  
    
   def f_angle(theta) :
       
      theta = normalise(theta ) # on le normalise pour etre dans le bon intervalle 
      
      pos_cylin [1] = theta  #ici je change theta a chaque appel de fonction
      
      #ensuite on reconverti en cartesienne 
      
      pos_carte = convert_Cy_to_Car(pos_cylin[0],pos_cylin[1], pos_cylin [2])#
      
      # on reinitialise la position dans y0
      y0[3] = pos_carte[0]
      y0[4] = pos_carte[1]
      y0[5] = pos_carte[2]
      
      #on refait la meme chose que dans recherche hauteur 
      solve = Tj.solve_ivp(
          Tj.oderhs,
          [t_i, t_f],
          y0,
          t_eval = instants_a_evaluer,
          events =event_hauteur ,
          rtol = 10**(-5),
          atol = 10**(-8))
  
      timetable = solve.y
      x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
      return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
      
      
   return R.bissection(f_angle, theta_min, theta_max , tol) #bissection car se comporte mieux 




def rechercheVitesse(y0,cibleRebond) : 
    
    
    #on convertie d'abord pour changer la norme  sans change la DIRECTION de la VITESSE (theta)
    
    
    pos_cylin = convert_Car_to_Cy(y0[3],y0[4],y0[5])
    
    r =    pos_cylin[0]
    theta= pos_cylin[1]  # fix donc pas reinitialise dans la fonction 
    z =    pos_cylin[2]
    
    def f_vitesse(norm_v) : 
       #erreur 
        #TODO probleme r
       # r =  ((norm_v**2)-z)/(Tj.np.cos(theta)+Tj.np.sin(theta)) #retrouve r a partir de la norme
        r = Tj.np.sqrt((norm_v**2)-(z**2))
        
        pos_carte = convert_Cy_to_Car(r,theta ,z)#on reconvertie en cartesienne
        
          
        # on reinitialise la position dans y0
        y0[3] = pos_carte[0]
        y0[4] = pos_carte[1]
        y0[5] = pos_carte[2]
        
        
        #on refait la meme chose que dans recherche hauteur 
        solve = Tj.solve_ivp(
            Tj.oderhs,
            [t_i, t_f],
            y0,
            t_eval = instants_a_evaluer,
            events =event_hauteur ,
            rtol = 10**(-5),
            atol = 10**(-8))
    
        timetable = solve.y
        x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
        return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
        
        
    return R.bissection(f_vitesse, vitesse_min, vitesse_max , tol) #bissection car se comporte mieux 


def rechercheOmega(y0,cibleRebond) :
        
    #on convertie d'abord pour changer la norme  sans change la DIRECTION de la VITESSE (theta)
    
    
    pos_cylin = convert_Car_to_Cy(y0[6],y0[7],y0[8])
    
    r =    pos_cylin[0]
    theta= pos_cylin[1]  # fix donc pas reinitialise dans la fonction 
    z =    pos_cylin[2]
    
    def f_omega(norm_w) : 
       
       # r =  ((norm_w**2)-z)/(Tj.np.cos(theta)+Tj.np.sin(theta)) #retrouve r a partir de la norme
        r = Tj.np.sqrt((norm_w**2)-(z**2))
        
        pos_carte = convert_Cy_to_Car(r,theta ,z)#on reconvertie en cartesienne
        
          
        # on reinitialise la position dans y0
        y0[6] = pos_carte[0]
        y0[7] = pos_carte[1]
        y0[8] = pos_carte[2]
        
        
        #on refait la meme chose que dans recherche hauteur 
        solve = Tj.solve_ivp(
            Tj.oderhs,
            [t_i, t_f],
            y0,
            t_eval = instants_a_evaluer,
            events =event_hauteur ,
            rtol = 10**(-5),
            atol = 10**(-8))
    
        timetable = solve.y
        x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
        return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
        
        
    return R.bissection(f_omega, omega_min, omega_max , tol) #bissection car se comporte mieux 

        
        
      
    
#deuxieme partie recherche hauteur apres un Rebond  

#Definition d'event
def stop (t , y) :
    if y[0]>Tj.distance_maximal_terrain:
        return 0
    else:
        return int(y[0]-Tj.distance_maximal_terrain)
stop.terminal = True
#stop.direction = 1
    
def rechercheHauteur2(initial_ball_data,cibleHauteur):

      
    def hauteur2(h):
        # définition des paramètres :
        initial_ball_data[2] = h
        frequence = 10**(5)
        instants_a_evaluer = Tj.np.linspace(0, t_f, int(frequence * t_f))
        

        
        pre_bounce_solve = Tj.solve_ivp(
            Tj.oderhs,
            [0, t_f],
            initial_ball_data,
            t_eval = instants_a_evaluer,
            events = Tj.bouing,
            rtol = 10**(-5),
            atol = 10**(-8))

        ball_data_timetable = pre_bounce_solve.y
        
        x_rebond = ball_data_timetable[0,-1]
        print ('la balle rebondit a ', x_rebond)
        
        if pre_bounce_solve.status == 1 : # = si il rebondi
            # manière d'aller chercher les dernières variables dans ball_data_timetable
            # TODO_LOW possiblement plus clair d'utilisation via pre_bounce_solve.t.shape[0]
            after_bounce_ball_data = ball_data_timetable[0:ball_data_timetable.shape[0], ball_data_timetable.shape[1] - 1]
            after_bounce_ball_data[5] *= Tj.coef
            
            # reinitialisation des données :
            # frequence * (delta_t) == on cree un tableau de frequence fixe allant du nouveau t_i au t_f donné 
            # ce qui assure une precision a chaque rebond quel que soit t_f donné
            
            t_i = pre_bounce_solve.t[-1] # avec -1, on lit le tableau à l'envers, donc on trouve le dernier élément
            delta_t = t_f - t_i
            instants_a_evaluer = Tj.np.linspace(t_i, t_f, int(frequence * delta_t))
            
            # on relance solve_ivp avec les nouvelles cond 
            
            post_bounce_solve =Tj.solve_ivp(
                Tj.oderhs,
                [t_i, t_f],
                after_bounce_ball_data,
                t_eval = instants_a_evaluer,
                events = [stop],
                rtol = 10**(-5),
                atol = 10**(-8))
               
            ball_data_timetable = Tj.np.concatenate((ball_data_timetable, post_bounce_solve.y), axis=1) # on regroupe les tableaux

        # debut :
        
       # x=ball_data_timetable[0]  #longeur
       # y=ball_data_timetable[1]  #largeur
        z=ball_data_timetable[2]  #hauteur
            

       # index_max =  Tj.np.shape(x)[0] - 1 # car shape(x) = shape(y) = shape(z)
  
        return  z[-1]-cibleHauteur
    return R.bissection( hauteur2,haut_min,haut_max, tol)
 
      
    

                               

