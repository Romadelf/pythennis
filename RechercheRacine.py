import numpy as np
def secante(f, x0, x1, tol):
   
    nmax = 1000
    niter=0
    statut=0
    
    # Est-ce que x0 et x1 sont reel?
    if np.isreal(x0)==False or np.isreal(x1)==False  :
        statut=1
        return [42, statut]
    
    fx0= f(x0)  #fait en sorte dappele le moin possible la fonction
    fx1= f(x1)
    #verification sur l'existance de x0 et x1
    if np.isreal(fx0)==False or np.isreal(fx1)==False:

        if(np.isreal(fx0) and abs(fx0)<=tol):  #ici on verifie si l'intervalle de  
            x=x0                               # depart est la solution
            liste= [x,statut] 
            return liste
        if(np.isreal(fx1) and abs(fx1)<=tol):
            x=x1
            liste= [x,statut] 
            return liste
        
        statut=1
        return [42, statut]  
    
    
    if x1 == x0 and abs(fx1) <= tol: # cas particulier ou l'intervalle est un point et que ce point et la solution
        return [x1,statut] 
    
    
    
    while (abs(fx1-fx0) > tol) and (niter < nmax):
        
        if abs(fx0-fx1) <= tol:  #G pour eviter davoir fx1-fx0 = 0 lors du calcul de z
               statut=1
               return [42, statut]
        
        z = x1-(fx1*(x1-x0)/(fx1-fx0))
        x0 = x1
        x1= z
        if (not(np.isreal(x1))): # parce que le nouveau x0 c'est l'ancien x1
            statut=1
            return [42, statut] 

        fx0= fx1
        fx1=f(x1)
        
        if (not(np.isreal(fx1))):
            statut=1
            return [42, statut] 
        
        niter += 1

    if niter==nmax:
      x= 42 # valeur derreur
      statut=-1
      return [x,statut]
    else : 
      x=x0
      return [x,statut] 
                  
#cas ou le programme a correctement fonctionné statut=0
#si x0 et x1 ne satisfont pas les hypothese statut=1
   #dans ce cas x contient n'importe quelle valeur 
#si la fonction ne possede pas de racine alors statut=-1 
   #dans ce cas x contient n'importe quelle valeur 
def bissection(f,x0, x1 , tol):
    
    statut=0
    
    if np.isreal(x0)==False or np.isreal(x1)==False :
        statut=1
        return [42, statut] 
    

    fx0= f(x0)  #fait en sorte dappele le moin possible la fonction
    fx1= f(x1)
    
    
    if np.isreal(fx0)==False or np.isreal(fx1)==False:

        if(np.isreal(fx0) and abs(fx0)<=tol):   
            x=x0
            liste= [x,statut] 
            return liste
        if(np.isreal(fx1) and abs(fx1)<=tol):
            x=x1
            liste= [x,statut] 
            return liste
        
        statut=1
        return [42, statut]
    

    if x1 == x0 and abs(fx1) <= tol: # cas particulier ou l'intervalle est un point et que ce point et la solution
        return [x1,statut] 
   
    
    
    if(fx1*fx0>0): # verification des signes oppose
        statut=1
        return [42, statut]

    
    
    while (abs(x0-x1)>=tol):# continuer tant que x0 et x1 ne sont pas suffisamment proche 
    
        c=(x0+x1)/2
        if f(x0)*f(c)< 0:
            x1=c
        else :
            x0=c
       
    if(abs(x0-x1)>=tol):
        statut=1
        return [42, statut]
    else :
        x=(x0+x1)/2  # pour un maximum de precision on renvoi la moyenne
        return [x,statut] 
        