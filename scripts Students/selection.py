#selection.py

"""========================================================================
Bibliothèque 
    decoupage et récupération des données/points/angles

--> donne les points de la zone choisie et l'angle associé
    ainsi que les points complémentaires et leur angle associé
    def anglepoints(x2,y2,z2,thetaMin,thetaMax,ym,yM):

========================================================================"""

import numpy as np


def anglepoints(x2,y2,z2,thetaMin,thetaMax,ym,yM):
    """
    permet de récuperer les points de la zone choisie, leur indice 
    ainsi que leur angle associé
    de plus permet de recuperer les points complémentaires, leur indice,
    et leur angle associé
    Input : x2,y2,z2 (les coordonnées de tous les points)
            thetaMin,thetaMax (angles minimum et maximum)
            ym,yM (les hauteurs min et max)
    Output : x3,y3,z3 (les coordonnées de la zone choisie)
             t3 (les indices de ces points pour le retour dans Meshlab)
             th3 (les angles theta associés)
             x3c,y3c,z3c (les points complémentaires)
             t3c (les indices de ces points pour le retour dans Meshlab)
    """
    # angles des points (x2,y2,z2):
    # on détermine l'angle theta associé à chaque point (x2,y2,z2)
    # danger divZ par0!!!! si la figure a été mal choisie
    # selection des points selon la coordonnées y et l'angle 
    #
    x3 = []
    y3 = []
    z3 = []
    t3 = []
    th3 = []
    x3c = []
    y3c = []
    z3c = []
    t3c = []

    for i in range(len(x2)):
        theta = np.arctan(x2[i]/z2[i])
        if thetaMin < theta < thetaMax and ym < y2[i] < yM and z2[i] > 0:
            x3.append(x2[i])
            y3.append(y2[i])
            z3.append(z2[i])
            t3.append(i)
            th3.append(theta)
        else:
            x3c.append(x2[i])
            y3c.append(y2[i])
            z3c.append(z2[i])
            t3c.append(i)
    return (x3,y3,z3,x3c,y3c,z3c,t3,th3,t3c)