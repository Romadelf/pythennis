 
import Trajectoire as Tj
import RechercheRacine as R
import cible as ci    

def rechercheAngle(y0,cibleRebond) :
   
#on convertie d'abord pour changer theta sans change la norme  de la VITESSE 
   pos_cylin = ci.convert_Car_to_Cy(y0[3],y0[4],y0[5])
  
  
    
   def f_angle(theta) :
      
      
      pos_cylin [1] = theta  #ici je change theta a chaque appel de fonction
      
      #ensuite on reconverti en cartesienne 
      
      pos_carte = ci.convert_Cy_to_Car(pos_cylin[0],pos_cylin[1], pos_cylin [2])#
      
      # on reinitialise la position dans y0
      y0[3] = pos_carte[0]
      y0[4] = pos_carte[1]
      y0[5] = pos_carte[2]
      
      #on refait la meme chose que dans recherche hauteur 
      solve = Tj.solve_ivp(
          Tj.oderhs,
          [ci.t_i,ci.t_f],
          y0,
          t_eval = ci.instants_a_evaluer,
          events =ci.event_hauteur ,
          rtol = 10**(-5),
          atol = 10**(-8))
  
      timetable = solve.y
      x1_rebond = timetable [0,-1]  #x1 est l'abcisse 
      return x1_rebond - cibleRebond #on retourn la difference pour savoir quand c'est zero 
      
      
   
   return R.bissection(f_angle, ci.theta_min, ci.theta_max , ci.tol) #bissection car se comporte mieux 

def rechercheAngle2(y0,cibleHauteur) :
   
#on convertie d'abord pour changer theta sans change la norme  de la VITESSE 
   pos_cylin = ci.convert_Car_to_Cy(y0[3],y0[4],y0[5])
   
  
    
   def f_angle2(theta) :
      
      
      pos_cylin [1] = theta  #ici je change theta a chaque appel de fonction
      
        # définition des paramètres :
      frequence = 10**(5) 
      instants_a_evaluer = Tj.np.linspace(0, ci.t_f, int(frequence * ci.t_f))
      #theta = normalise(theta ) # on le normalise pour etre dans le bon intervalle 
     
        
    
      #ensuite on reconverti en cartesienne 
      
      pos_carte = ci.convert_Cy_to_Car(pos_cylin[0],pos_cylin[1], pos_cylin [2])#
      
      # on reinitialise la position dans y0
      y0[3] = pos_carte[0]
      y0[4] = pos_carte[1]
      y0[5] = pos_carte[2]
      
        
      pre_bounce_solve = Tj.solve_ivp(
            Tj.oderhs,
            [0, ci.t_f],
            y0,
            t_eval = instants_a_evaluer,
            events = Tj.bouing,
            rtol = 10**(-5),
            atol = 10**(-8))

      ball_data_timetable = pre_bounce_solve.y
        
      x_rebond = ball_data_timetable[0,-1]
     # print ('la balle rebondit a ', x_rebond)
        
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

        #  debut :
        
       # x=ball_data_timetable[0]  #longeur
       # y=ball_data_timetable[1]  #largeur
      z=ball_data_timetable[2]  #hauteur
            

       # index_max =  Tj.np.shape(x)[0] - 1 # car shape(x) = shape(y) = shape(z)
  
      
      return  z[-1]-cibleHauteur
   t = R.secante( f_angle2,ci.theta_min, ci.theta_max , ci.tol)
   t[0] = ci.normalise(t[0])
   return t


