import numpy as np
from scipy.integrate import solve_ivp

#definition des conditions initiales
x0_x1 =-11.89
x0_x2 = 0
x0_x3 =2

dx01_dt= 50
dx02_dt= 1
dx03_dt= 0

w0_x1= 30 
w0_x2= 15
w0_x3= 0

Cond= [x0_x1,x0_x2,x0_x3,dx01_dt,dx02_dt,dx03_dt,w0_x1,w0_x2,w0_x3]

t0,T=0,1
tol = 1/10000  #tolerance 
coef = 0.7   # coef de changement de vz
compteur_max = 50
nombre_iteration = np.linspace(t0,T,compteur_max)


def oderhs(t, y):
    
    #definition des parametre du systeme
    d=0.0065
    m=0.058
    p=1.2
    C_d=0.065
    g=9.81
    
    #ici on initialise notre repere orthonormé
    Vec_u_x1=[1,0,0] #x
    Vec_u_x2=[0,1,0] #y
    Vec_u_x3=[0,0,1] #z
    
    #x1= y[0]
    #x2= y[1]      # ici on intialise la position car y est un vecteur   d'indice 0,1,2
    x3= y[2] 
    
    
    #les forces  fdx et fm  a savoir les projection :
    
    #ici on prend les trois composante de la vitesse (angulaire et de translation)
    v_x1 = y[3]
    v_x2 = y[4]
    v_x3 = y[5]
    
#   if y[5] < 0 and x3 <= tol:                  # ici on créé les rebonds
#                                               # si la hauteur est inferieur a 0 
#      v_x3 = (-coef)  * y[5]     # alors on inverse a vitesse selon z
#   else :
#     v_x3 = y[5]
        
    
    w_x1 = y[6]
    w_x2 = y[7]
    w_x3 = y[8]
   
   
   #puis on calcule leurs normes
   
   
    vitesse_translation=[v_x1, v_x2 ,v_x3]
    vitesse_angulaire  =[w_x1, w_x2 ,w_x3]
    
    
    norme_de_v= np.linalg.norm(vitesse_translation)
    norme_de_w= np.linalg.norm(vitesse_angulaire)
   
    #NORME DES FORCES 
    
    Force_frottement= C_d*p*(np.pi)*((d*d)/8)*((norme_de_v)**2)
    
    C_m=1/(2+1.96*((norme_de_v)/(norme_de_w)*d))
    
    Force_magnus=C_m*p*(np.pi)*((d*d)/8)*((norme_de_v)**2)
    
    #ici on a un vecteur unitaire de la vitesse de translation
    
    direction_vitesse_translation= (vitesse_translation)/(norme_de_v)
    
    
    # mainenant on determine les composantes grace au produit scalaire des deux 
    #vecteur unitaire  pour obtenir nos trois composante
    #ON FAIT (-)  la direction car c'est l'inverse de la vitesse
    
    aux= np.dot(- direction_vitesse_translation,Vec_u_x1)
    fdx1= (Force_frottement)*(aux)
    aux= np.dot(- direction_vitesse_translation,Vec_u_x2)
    fdx2= (Force_frottement)*(aux)
    aux= np.dot(- direction_vitesse_translation,Vec_u_x3)
    fdx3= (Force_frottement)*(aux)
    
    
    
    #fmx1? fmx2? fmx3? en fonction de leur direction
    
    #ici on fait le produit vectorielle w vectoriel v
    #puis on le norme 
    
    Produit_vectoriel = np.cross(vitesse_translation,vitesse_angulaire)
    
    Norme_Produit_vectoriel=np.linalg.norm(Produit_vectoriel)
    
    direction_force_magnus=(Produit_vectoriel)/(Norme_Produit_vectoriel)                 
    
    aux1=np.dot(direction_force_magnus,Vec_u_x1)
    fmx1= (Force_magnus)*aux1
    aux1=np.dot(direction_force_magnus,Vec_u_x2)
    fmx2= (Force_magnus)*aux1
    aux1=np.dot(direction_force_magnus,Vec_u_x3)
    fmx3= (Force_magnus)*aux1
    
    #Poid 
    
    fgx1 = 0
    fgx2 = 0
    fgx3 =-m*g
    

    #definition du systeme
    
    # on pose  X= dx1_dt              dX_dt = acceleration selon x1
             # Y= dx2_dt   ainsi on a dY_dt = acceleration selon x2
             # Z= dx3_dt              dZ_dt = acceleration selon x3
   # vecteur y (x1,x2,x3,dx1_dt,dx2_dt,dx3_dt,w1,w2,w3)
   # indice     0  1  2  3      4       5     6  7  8

   
   
  
    dx1_dt = v_x1  # ici on  reprend les variables deja pris au debut pour les normes
    dx2_dt = v_x2  # on prend les dx_dt pour juste respecté les conventions
    dx3_dt = v_x3
    
    
    
    dy = np.zeros(9)
    
    
    dy[0]= dx1_dt
    dy[1]= dx2_dt      # on met les trois composantes de la vitesse, 
    dy[2]= dx3_dt
    
    dy[3]=(fgx1+fdx1+fmx1)*1/m
    dy[4]=(fgx2+fdx2+fmx2)*1/m      # les trois composante de l'acceleration 
    dy[5]=(fgx3+fdx3+fmx3)*1/m
    
    #pas necessaire de remvoyer omega car tableau initialisé a zero
    return dy

def bouing(t, y):
    return y[2]
bouing.direction = -1
bouing.terminal = True


def main():

    
    variables_0 = solve_ivp(oderhs,[t0,T],Cond ,t_eval = nombre_iteration, events = bouing)
    arr_0 = variables_0.y
    
    if variables_0.status == 1 : # = si il rebondi
    
        new_cond = arr_0[0:arr_0.shape[0], arr_0.shape[1] - 1] # manière d'aller chercher les dernières variables dans arr
        new_cond[5] = - coef * new_cond[5]
        
        variables_1 = solve_ivp(oderhs,[0,1], new_cond , t_eval = np.linspace(0,1,50))
        arr_1 = variables_1.y
        
        return np.concatenate((arr_0, arr_1), axis=1)
        
        

        
def trajectoireFiletHorizontal (yInit , T ):
    
    #initialisation des condition initial :
        
    x0_x1  = yInit[0]
    x0_x2  = yInit[1]
    x0_x3  = yInit[2]
    dx01_dt= yInit[3]
    dx02_dt= yInit[4]
    dx03_dt= yInit[5]
    w0_x1  = yInit[6]
    w0_x2  = yInit[7]
    w0_x3  = yInit[8]
    
    
    #definition de parametre :
        
    t0=0
    compteur_max = 200
    i = 0
    haut_filet = 1 # cas ou filet droit 
    tol = 0.001
    D = 11.89 # en m
    
    Cond= [x0_x1,x0_x2,x0_x3,dx01_dt,dx02_dt,dx03_dt,w0_x1,w0_x2,w0_x3]

    #verification des donnés :
        
        #autrement dit si ya pas de vitesse(translation et rotation) initial  ya pas de mouvements
        #on renvoit la position initial
          
    if (abs(dx01_dt)<=tol) and (abs(dx02_dt)<=tol) and (abs(dx03_dt)<=tol):
        
        return  [x0_x1,x0_x2,x0_x3] #position initial  au cas ou la vitesse initial est nul
    
    if (abs(w0_x1)<=tol) and (abs(w0_x2)<=tol) and (abs(w0_x3)<=tol): #cas ou on w0 = 0 est evite les problemes
        
        return  [0,0,0] #valeur  d'erreur au cas ou la vitesse initial est nul
  
    

    
    # debut :
        
        
    nombre_iteration = np.linspace(t0,T,compteur_max) #on cree un tableau de temp 

    variable = solve_ivp(oderhs,[t0,T],Cond) #pas besoin de t_eval car tout fait dans TSPAN

    #temp = variable.t pas necessaire 
    
    position = variable.y
  
    x1=position[0]  #longeur
    x2=position[1]  #largeur
    x3=position[2]  #hauteur
    
    
    
    #passe le filet   ou pas 
    
    
    #voir event de solve ivp
    
    for i in range(compteur_max):
        
        if x1[i]<=tol:  #si on est derriere le filet 
            
            if x3[i]<=(haut_filet + tol): #et que la hauteur de la balle est inferieur a cekke du filet 
                                         #tol pour eviter les erreur de calcul
                return [0,0,0] #erreur 
        else :  #ici ca veut dire qu'on est a droite du filet  x1[i]>tol
          
            if x1[i]> D+tol: #on est hors du terrain
                
                return[0,0,0]
            else :
                
                
                Dernier_indice = compteur_max-1  #parce que indice va de 0 a n-1
                
                return [x1[Dernier_indice],x2[Dernier_indice],x3[Dernier_indice]]
          
           
  







