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
import numpy as np
#Condition initial : 

y0 = [-11.89,0,2,50,1,0,30,15,0]    

#parametre : 
    
t_i , t_f = 0, 4  #t_span 
n_point = 200  #nombre de point 
instants_a_evaluer = np.linspace(t_i, t_f,n_point ) 
  
statut  , error = 0 , 42
tol =  1/1000 #tolerance au millimetre 
pas = 1/100 #1/100 DE METRE donc du centimetre 
haut_min =2 # m
haut_max =3 # m

z0 , z1 = 0 , 10



      
def rechercheHauteur(y0,cibleRebond): 
    #argument de fonction : 
    def f_trajectoire(h):
        #Definition d'event
        def event_hauteur (t,y0):
            return y0[2]
        event_hauteur.direction = -1
        event_hauteur.terminal = True

    
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
    return R.bissection(f_trajectoire, z0, z1, tol) #bissection car se comporte mieux 
        
        

    
  
        
      
    

                               

