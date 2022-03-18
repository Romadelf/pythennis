#angle, statut = rechercheAngle(y0,cibleRebond)
#angle, statut = rechercheAngle2(y0,cibleHauteur)
#normeVitesse, statut = rechercheVitesse(y0,cibleRebond)
#normeVitesse, statut = rechercheVitesse2(y0,cibleHauteur)
#hauteur, statut = rechercheHauteur(y0,cibleRebond)
#hauteur, statut = rechercheHauteur2(y0,cibleHauteur)
#omega, statut = rechercheOmega(y0,cibleRebond)
#omega, statut = rechercheOmega2(y0,cibleHauteur)
import Trajectoire as Tj
import RechercheRacine as sol
import numpy as np

def rechercheHauteur(y0,cibleRebond): 
    t_max = 10 
    haut_init = y0[2]
    position =  Tj.trajectoireFiletHorizontal(y0,t_max )
    pas = 0.1
    haut_max =5
    tol = 1/100
    i =1
    n_max=100
  
    
    while  abs(position[0]-cibleRebond)>tol :
    
        if position[0]<=cibleRebond:
            
            y0[2] += pas
            position =  Tj.trajectoireFiletHorizontal(y0,t_max )
        else :
            y0[2] -= pas
            position =  Tj.trajectoireFiletHorizontal(y0,t_max ) 
            
    hauteur =position[2]
            
    
    return hauteur
    
#event = lambda t, y: y[2]