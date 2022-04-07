import Trajectoire as Tj
import matplotlib.pyplot as plt

def trajectoireFiletHorizontal2(y0, T):
    
    #parametre :
    h = 1/100     #définition du pas
    t0  = 0                #temps initial
    tf = T                #temps final
    i = 0 
    #nombre d'itérations max
    
    #round arrondit a l'entier le plus proche
    N = round((tf - t0 ) / h)
    
    #initialisation des tableaux de positions à la premiere case avec les valeurs initiales
    x1, x2, x3 = [y0[0]], [y0[1]], [y0[2]]
    
    
    #initialisation d'un tableau de temps
    t = [t0]
    #calcul de trajectoire grâce à Euler
    for k in range(N):
        t0 += h
        #formule d'Euler
        y0 += h * Tj.oderhs(t0, y0)

        #ajout des nouvelles valeurs de yInit au tableau Lx1, Lx2, Lx3 créé précédement
        x1.append(y0[0])
        x2.append(y0[1])
        x3.append(y0[2])
        #ajout du temps ou tableau Lt
        t.append(t0)

        #vérification la balle est dans le filet

        if y0[2] <= 0 and y0[5] < 0:
          #  print("Rebond")
            #changement de la vitesse initial selon z par un coefficient de rebond
            y0[5] = Tj.coef*y0[5]
            #ca fait tout les rebond essaye avec coef egal -0.8  
        index_max = Tj.np.shape(x1)[0] - 1 # car shape(x) = shape(y) = shape(z)    
        while i <= index_max and x1[i] <= 0: # tant que  on est avant le filet 
              if  x3[i] <= Tj.hauteur_filet(x2[i]) :  # et si la hauteur de la balle est inferieur a celle du filet 
                 #  print(' la balle touche le filet ')   
                   return [0,0,0] # hors-jeu 
              i += 1
          
    plt.plot(x1,x3)
    plt.title('trajectoire') #Ajout d'un titre
    plt.ylabel('Axe des Y') # Labélisation de l'axe des ordonnées
    plt.xlabel('Axe des X')
   
    return [y0[0], y0[1], y0[2]] 
