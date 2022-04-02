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

haut_min,haut_max =  1 , 100 # m
theta_min, theta_max = 0, (Tj.np.pi)/2
theta_min2, theta_max2 = -(Tj.np.pi)/2 , (Tj.np.pi)/2
vitesse_min, vitesse_max = 1 , 2000
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
           angle -=Tj.np.pi 
    angle_normalise =angle 
    
    return angle_normalise 



#deuxieme partie recherche hauteur apres un Rebond  

#Definition d'event
def stop (t , y) :
    if y[0]>Tj.distance_maximal_terrain:
        return 0
    else:
        return int(y[0]-Tj.distance_maximal_terrain)
stop.terminal = True
#stop.direction = 1
    



 

        

      
    

                               

