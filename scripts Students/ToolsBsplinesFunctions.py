#ToolsBsplinesFunctions

"""========================================================================
Bibliothèque
    de fonctions pour la reconstruction géométrique
    des données issues de la kinect :

--> fonctions des coefficients:
    def s00(t0,t1,t) _ s10/s11(t0,t1,t2,t) et s20/21/22(t0,t1,t2,t3,t)
    et si0/i1/i2/i3(a,b,c,d,e,t)

--> fonctions des différentes splines selon les noeuds:
    def Bspline_one(t,a,b,cond) et Bspline_two(t,a,b,c,cond) et
    Bspline_three(t,a,b,c,d,cond) et Bspline_ieme(t,a,b,c,d,e) 
     
--> Evaluation d'une Spline :
    def Bspline(tk,i,t)

========================================================================"""

import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------
#défini les coefficients
def s00(t0,t1,t):
    y = ((t1-t)/(t1-t0))**3
    return y

def s10(t0,t1,t2,t):
    y = ((t-t0)*((t1-t)**2))/((t1-t0)**3) \
        + ((t2-t)*(t1-t)*(t-t0))/((t2-t0)*((t1-t0)**2)) \
        + (((t2-t)**2)*(t-t0))/(((t2-t0)**2)*(t1-t0))
    return y

def s11(t0,t1,t2,t):
    y = ((t2-t)**3)/(((t2-t0)**2) * (t2-t1))
    return y

def s20(t0,t1,t2,t3,t):
    y = (((t-t0)**2)*(t1-t))/((t2-t0)*((t1-t0)**2)) \
        + (((t-t0)**2)*(t2-t))/(((t2-t0)**2)*(t1-t0)) \
        + ((t3-t)*((t-t0)**2))/((t3-t0)*(t1-t0)*(t2-t0))
    return y

def s21(t0,t1,t2,t3,t):
    y = ((t-t0)*((t2-t)**2))/(((t2-t0)**2)*(t2-t1)) \
        + ((t3-t)*(t-t0)*(t2-t))/((t3-t0)*(t2-t0)*(t2-t1)) \
        + (((t3-t)**2)*(t-t1))/((t3-t0)*(t3-t1)*(t2-t1))
    return y

def s22(t0,t1,t2,t3,t):
    y = ((t3-t)**3)/((t3-t0)*(t3-t1)*(t3-t2))
    return y

def si0(a,b,c,d,e,t):
    y = ((t-a)**3)/((d-a)*(c-a)*(b-a))
    return y

def si1(a,b,c,d,e,t):
    y = (((t-a)**2)*(c-t))/((d-a)*(c-a)*(c-b)) \
        + ((t-a)*(t-b)*(d-t))/((d-a)*(d-b)*(c-b)) \
        + (((t-b)**2)*(e-t))/((e-b)*(d-b)*(c-b))
    return y

def si2(a,b,c,d,e,t):
    y = ((t-a)*((d-t)**2))/((d-a)*(d-b)*(d-c)) \
        + ((t-b)*(e-t)*(d-t))/((e-b)*(d-b)*(d-c)) \
        + ((t-c)*((e-t)**2))/((e-b)*(e-c)*(d-c))
    return y

def si3(a,b,c,d,e,t):
    y = ((e-t)**3)/((e-b)*(e-c)*(e-d))
    return y

#-------------------------------------------------------------------
#défini les fonctions pour calculer les courbes
def Bspline_one(t,a,b,cond):
    """
    défini la première et dernière spline
    Input: t (valeur réelle comprise entre tk[0] et tk[n-1]), a et b (des noeuds)
        et cond (0 ou 1 pour savoir si c'est une spline du début ou de la fin)
    Output: une valeur réelle
    """
    
    if (a <= t <=b):###########
        if (cond == 0):
            #première courbe
            return s00(a,b,t)
        if (cond==1):
            #dernière courbe
            return s00(b,a,t)
    else:
        #entre tous les autres noeuds la courbe est nulle
        return 0

def Bspline_two(t,a,b,c,cond):
    """
    défini la seconde et l'avant dernière spline
    Input: t (valeur réelle comprise entre tk[0] et tk[n-1]), a et b (des noeuds)
        et cond (0 ou 1 pour savoir si c'est une spline du début ou de la fin)
    Output: une valeur réelle
    """
    if (a <= t < b):
        #première partie de la courbe
        if (cond ==0):
            #seconde
            return s10(a,b,c,t)
        if (cond ==1):
            #avant dernière
            return s11(c,b,a,t)
    elif (b <= t < c):
        #deuxième partie de la courbe
        if (cond == 0):
            #seconde
            return s11(a,b,c,t)
        if (cond == 1):
            #avant dernière
            return s10(c,b,a,t)
    else:
        #entre tous les autres noeuds la courbe est nulle
        return 0

def Bspline_three(t,a,b,c,d,cond):
    """
    défini la troisième et l'avant-avant dernière spline
    Input: t (valeur réelle comprise entre tk[0] et tk[n-1]), a et b (des noeuds)
        et cond (0 ou 1 pour savoir si c'est une spline du début ou de la fin)
    Output: une valeur réelle
    """
    if (a <= t < b):
        #première partie de la courbe
        if (cond == 0):
            #troisième
            return s20(a,b,c,d,t)
        elif (cond == 1):
            #avant-avant dernière
            return s22(d,c,b,a,t)
    elif (b <= t < c):
        #deuxième partie de la courbe
        if (cond == 0):
            #troisième
            return s21(a,b,c,d,t)
        elif (cond == 1):
            #avant-avant dernière
            return s21(d,c,b,a,t)
    elif (c <= t < d):
        #troisième partie de la courbe
        if (cond == 0):
            #troisième
            return s22(a,b,c,d,t)
        elif (cond == 1):
            #avant-avant dernière
            return s20(d,c,b,a,t)
    else:
         #entre tous les autres noeuds la courbe est nulle
        return 0

def Bspline_ieme(t,a,b,c,d,e):
    """
    défini la i-ème spline
    Input: t (valeur réelle comprise entre tk[0] et tk[n-1]), a et b (des noeuds)
        et cond (0 ou 1 pour savoir si c'est une spline du début ou de la fin)
    Output: une valeur réelle
    """
    if a <= t < b:
        #première partie
        return si0(a,b,c,d,e,t)
    elif b <= t < c:
        #deuxième partie
        return si1(a,b,c,d,e,t)
    elif c <= t < d:
        #troisième paartie
        return si2(a,b,c,d,e,t)
    elif d <= t < e:
        #quatrième partie
        return si3(a,b,c,d,e,t)
    else:
        #entre tous les autres noeuds la courbe est nulle
        return 0
    

#-------------------------------------------------------------------
#défini la Bspline associé
def Bspline(tk,i,t):
    """
        Calcul de la fonction B-spline N^3_i(t) pour un scalaire t
        Input : 
            tk = vecteur des noeuds = np.array([]) de taille n = len(tk)
            RMQ : il faut n >= 6
            t valeur réelle comprise entre tk[0] et tk[n-1]
            i = entier = numéro de la fonction B-spline concernée 
            avec 0 <= i <= n+1
        Output :
            valeur en t de la fonction B spline cubique numéro i 
            associé à la suite de noeuds tk
    """
    n = len(tk)
    #pour définir la première et dernière courbe
    if ((i == 0) or (i == n+1)):
        if (i == 0) :
            t0 = tk[0]
            t1 = tk[1]
            return Bspline_one(t,t0,t1,0)
        if (i == n + 1) :
            t0 = tk[i-3]
            t1 = tk[i-2]
            return Bspline_one(t,t0,t1,1)

    # pour définir la deuxième et avant dernière courbe
    if ((i == 1) or(i == n)):
        if (i == 1):
            t0 = tk[0]
            t1 = tk[1]
            t2 = tk[2]
            return Bspline_two(t,t0,t1,t2,0)
        if (i == n):
            t0 = tk[i-3]
            t1 = tk[i-2]
            t2 = tk[i-1]
            return Bspline_two(t,t0,t1,t2,1)

    # pour définir la troisième et avant-avant dernière courbe
    if ((i == 2) or (i == n-1)):
        if (i == 2):
            #pour la troisième
            t0 = tk[0]
            t1 = tk[1]
            t2 = tk[2]
            t3 = tk[3]
            return Bspline_three(t,t0,t1,t2,t3,0)
        if (i == n-1):
            #pour l'avant-avant dernière
            t0 = tk[i-3]
            t1 = tk[i-2]
            t2 = tk[i-1]
            t3 = tk[i]
            return Bspline_three(t,t0,t1,t2,t3,1)
        
    #pour les autres courbes
    if ((i >= 3) and (i < n-1)):
        t0 = tk[i-3]
        t1 = tk[i-2]
        t2 = tk[i-1]
        t3 = tk[i]
        t4 = tk[i+1]
        return Bspline_ieme(t,t0,t1,t2,t3,t4)


