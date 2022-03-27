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
  
statut  , error = 0 , 42
tol =  1/1000 #tolerance au millimetre 


#intervalle 

haut_min,haut_max =  0 , 80 # m
theta_min, theta_max = 0 , (Tj.np.pi)/2
vitesse_min, vitesse_max = 1 , 500
omega_min, omega_max = 1 , 1000

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
    
    
    
def rechercheHauteur2(y0,cibleHauteur):
    
    #Definition d'event
    def event_apres_rebond(t,m):
       
              return m[2]-cibleHauteur
             
    event_apres_rebond.direction = 1
    event_apres_rebond.terminal = True
   
    m = Tj.np.zeros(9)
    
    #argument de fonction : 
    def f_hauteur(h):
       
          
        y0[2] =  h 
        solve = Tj.solve_ivp(
            Tj.oderhs,
            [t_i, t_f],
            y0,
            t_eval = instants_a_evaluer,
            events =event_hauteur,
            rtol = 10**(-5),
            atol = 10**(-8))
    
        timetable = solve.y

       
        if solve.status == 1:
            
           
            m[0] = timetable [0,-1]
            m[1] = timetable [1,-1]
            m[2] = timetable [2,-1]
            m[3] = timetable [3,-1]
            m[4] = timetable [4,-1]
            m[5] = Tj.coef*timetable [5,-1] #le, fameux coefficient egal a -0.7
            m[6] = timetable [6,-1]
            m[7] = timetable [7,-1]
            m[8] = timetable [8,-1]


            t_in = solve.t[-1] #on continue a partir de l'instant precedent 
            
            instants_a_evaluer_2 = Tj.np.linspace(t_in, t_f,n_point_apres)#on recree un tableau
            solve_apres_rebond = Tj.solve_ivp(
                Tj.oderhs,
                [t_in, t_f],
                m,
                t_eval = instants_a_evaluer_2,
                events =event_apres_rebond,
                rtol = 10**(-5),
                atol = 10**(-8))
            
            
            if solve_apres_rebond.status == 1:
               
                Nouveau_timetable  = solve_apres_rebond.y
                timetable = Nouveau_timetable
           
      
        x3_hauteur = timetable [2,-1]  #x3 est la hauteur 
        
        return x3_hauteur - cibleHauteur
        
           
        
      
    return R.bissection(f_hauteur, haut_min, haut_max , tol) #bissection car se comporte mieux 


    

   
      
              
    

        
      
    

                               

