#E = (1/2)*m*v**2  + (1/8)*m*(d**2)*(w**2)

 
import Trajectoire as Tj
import cible as ci    

#☺essai
vec3 = [1 , 2 ,0]

#donné static
x = [-11.89 ,0 ,  2] 

   # vec = [4,(Tj.np.pi)/4, 6]  #cylindrique
vec1 = [4,(Tj.np.pi)/2, 6]
    # vecteur de norme 4          #choix arbitraire
    # angle de pi/2 avec le sol
    #hauteur de 6
   # vec = ci.convert_Cy_to_Car(vec[0], vec[1], vec[2])   #cartesienne
vec1 = ci.convert_Cy_to_Car(vec1[0], vec1[1], vec1[2])
    

def methode (E , vec):
    
    #condition initial 
   
    frequence = 10**(5)
    instants_a_evaluer = Tj.np.linspace(0, ci.t_f, int(frequence * ci.t_f))
    
   # simpl = 2*E/(Tj.m) #version simplifieer
    #ainsi on a  : simpl = v**2 + 0.25*d*(w**2)
    #vecteur omega et vitesse ne doivent pas etre parrallele sinon erreur
    


    #verification :
        
    Prod_vec = Tj.np.cross(vec,vec1)
    norm_Prod_vec = Tj.np.linalg.norm(Prod_vec)
    if norm_Prod_vec == 0 :
        print('la direction choisi est parallele a la direction de omega ')
        print("veuillez modifié votre vecteur")
        return [42 , -1]
    #parametre 
    pat = 0.05    
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
       
       
       vx1 = Tj.np.dot(Tj.vec_u_x1,vec) #vecteur v
       vx2 = Tj.np.dot(Tj.vec_u_x2,vec)
       vx3 = Tj.np.dot(Tj.vec_u_x3,vec)
       v = [vx1*norm_v , vx2*norm_v , vx3*norm_v]
       
       
       wx1 = Tj.np.dot(Tj.vec_u_x1,vec1) #vecteur w
       wx2 = Tj.np.dot(Tj.vec_u_x2,vec1)
       wx3 = Tj.np.dot(Tj.vec_u_x3,vec1)
       w = [wx1*norm_w , wx2*norm_w , wx3*norm_w]
     
       y = Tj.np.zeros(9)
       
       y[0] = x[0]
       y[1] = x[1]
       y[2] = x[2]
       
       y[3] = v[0]
       y[4] = v[1]
       y[5] = v[2]
       
       y[6] = w[0]
       y[7] = w[1]
       y[8] = w[2]
       
    
      
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
       coef+= pat
    
    
    
    return bon_coef 

def verification (coef , E , vec , vec1):
    
    
    norm_v = Tj.np.sqrt(2*coef*E/Tj.m)                 #norme
    norm_w = Tj.np.sqrt((1-coef)*8*E/(Tj.m*(Tj.d**2)))   #norme
    
    
    vx1 = Tj.np.dot(Tj.vec_u_x1,vec) #vecteur v
    vx2 = Tj.np.dot(Tj.vec_u_x2,vec)
    vx3 = Tj.np.dot(Tj.vec_u_x3,vec)
    v = [vx1*norm_v , vx2*norm_v , vx3*norm_v]
    
    
    wx1 = Tj.np.dot(Tj.vec_u_x1,vec1) #vecteur w
    wx2 = Tj.np.dot(Tj.vec_u_x2,vec1)
    wx3 = Tj.np.dot(Tj.vec_u_x3,vec1)
    w = [wx1*norm_w , wx2*norm_w , wx3*norm_w]
  
    y = Tj.np.zeros(9)
    
    y[0] = x[0]
    y[1] = x[1]
    y[2] = x[2]
    
    y[3] = v[0]
    y[4] = v[1]
    y[5] = v[2]
    
    y[6] = w[0]
    y[7] = w[1]
    y[8] = w[2]
    
    Tj.trajectoireFiletHorizontal(y, ci.t_f)
 
    #faut que le coef soit tres important parce que sinon ca mene a des absurdite 
    #en effet si w represente ne serais-ce que 1 pour 100 alors la vitesse de rotation
    #sera tres importante 
    return 
    