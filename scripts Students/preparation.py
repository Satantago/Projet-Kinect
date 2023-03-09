# preparation.py

"""========================================================================
Bibliothèque
    de fonctions pour l'acquisition et la mise en place des données

--> Récupération et centralisation des données du fichier Mesh
    def acquisition(nom_fichierPly):

--> barycentre à l'origine pour centrer le visage sur nos axes
    def barycentre (x,y,z,NbPts):
        
--> Selection des points de devant (de la figure)
    def pointsdevant(X0,Y0,Z0):

--> rotation d'un angle autour de l'axe x sur les données de devant
    def rotationX(xu,yu,xA,yA,zA,X0,Y0,Z0):
        
--> rotation d'un angle autour de l'axe y sur les données en avant
    def rotationY(xu,yu,xA,yA,zA,X1,Y1,Z1):    
========================================================================"""

from ToolsMesh import *


"""-----------------------------------------------------
MISE EN PLACE
-----------------------------------------------------"""

def acquisition(nom_fichierPly):
    """ 
    On récupère la tête entière, avec probablement des trous à l'arrière
    la sélection de la zone à approximer est faite dans Python
    Input : Nom d'un fichier
    Output : X0,Y0,Z0,Nbtps les points du fichier recentrés
    """

    FilePly = str(nom_fichierPly) # attention au chemin
    x,y,z = CoordonneesFromPly(FilePly)
    NbPts = len(x)
    print("Nombre de points : ",NbPts)
    return (barycentre(x,y,z,NbPts))


def barycentre (x,y,z,NbPts):
    """
    barycentre à l'origine pour centrer le visage sur nos axes
    Input : x,y,z,Nbtps : des points et le nombre de points associé
    Output : X0,Y0,Z0,NbPts : les points recentrés et le nombre de points
    """
    
    # Barycentre (xb,yb,zb)
    xb = sum(x)/NbPts
    yb = sum(y)/NbPts
    zb = sum(z)/NbPts
    # translation pour ramener le barycentre à l'origine
    X0 = [xi - xb for xi in x]
    Y0 = [yi - yb for yi in y]
    Z0 = [zi - zb for zi in z]
    return (X0,Y0,Z0,NbPts)
    

"""-----------------------------------------------------
SELECTION DES POINTS DE DEVANT
-----------------------------------------------------"""

def pointsdevant(X0,Y0,Z0):
    """
    Pour mieux visualiser et sélectionner les rotations à effectuer 
    on sélectionne les points de "devant", cad tels que zA > 0.
    Input :  X0,Y0,Z0 = les points 
    Output : xA,yA,zA = les points de "devant"
    """
    NbPts = len(X0)
    
    # sélection des points de "devant" (pour l'affichage 2D)
    xA = []
    yA = []
    zA = []
    for i in range(NbPts):
        if Z0[i] > 0 :
            xA.append(X0[i])
            yA.append(Y0[i])
            zA.append(Z0[i])
    return (xA,yA,zA)


"""-----------------------------------------------------
FONCTIONS DE ROTATIONS
-----------------------------------------------------"""

def rotationX(xu,yu,xA,yA,zA,X0,Y0,Z0):
    """
    rotation autour de l'axe x
    Input : xu,yu (vecteur qui défini l'angle de rotation : 
                   on veut que ce vecteur soit horizontal)
            xA,yA,zA  (les points de devant pour l'affichage)
            X0,Y0,Z0  (l'entièreté des points à faire bouger)
    Output : Xdisp,Ydisp,Zdisp (les points d'affichage)
             X1,Y1,Z1 (l'entiereté des points après rotation)
    """
    
    Vx = xu[1]-xu[0]
    Vy = yu[1]-yu[0]
    # Vx supposé non nul, sinon pas besoin de rotation
    if Vy == 0 :
        return (xA,yA,zA,X0,Y0,Z0)
    else :
        costeta = Vx / (Vx*Vx + Vy*Vy)**.5
        sinteta = -Vy / (Vx*Vx + Vy*Vy)**.5
    # on fait d'abord  tourner uniquement la sélection (xA,yA,zA) pour controle
    Xdisp = xA
    Ydisp = [yA[i]*costeta - zA[i]*sinteta for i in range(len(yA))]
    Zdisp = [yA[i]*sinteta + zA[i]*costeta for i in range(len(yA))]
    #on affichera ces données
    
    #On fait la rotation maintenant sur l'ensemble des données
    X1 = X0
    Y1 = [Y0[i]*costeta - Z0[i]*sinteta for i in range(len(Y0))]
    Z1 = [Y0[i]*sinteta + Z0[i]*costeta for i in range(len(Y0))]
    
    return(Xdisp,Ydisp,Zdisp,X1,Y1,Z1)
    

def rotationY(xu,yu,xA,yA,zA,X1,Y1,Z1):
    """
    rotation autour de l'axe y
    Input : xu,yu (vecteur qui défini l'angle de rotation : 
                   on veut que ce vecteur soit vertical)
            xA,yA,zA  (les points de devant pour l'affichage)
            X1,Y1,Z1  (l'entièreté des points à faire bouger)
    Output : Xdisp,Ydisp,Zdisp (les points d'affichage)
             X2,Y2,Z2 (l'entiereté des points après rotation)
    """
    
    Vx = xu[1]-xu[0]
    Vy = yu[1]-yu[0]
    # Vx supposé non nul, sinon pas besoin de rotation
    if Vx == 0 :
        return (xA,yA,zA,X1,Y1,Z1)
    else :
        costeta = Vy / (Vx*Vx + Vy*Vy)**.5
        sinteta = Vx / (Vx*Vx + Vy*Vy)**.5

    # on fait d'abord  tourner uniquement la sélection (xA,yA,zA) pour controle
    Ydisp = yA
    Xdisp = [xA[i]*costeta - zA[i]*sinteta for i in range(len(yA))]
    Zdisp = [xA[i]*sinteta + zA[i]*costeta for i in range(len(yA))]
    #On fait maintenant la rotation sur l'ensemble des elements
    Y2 = Y1
    X2 = [X1[i]*costeta - Z1[i]*sinteta for i in range(len(Y1))]
    Z2 = [X1[i]*sinteta + Z1[i]*costeta for i in range(len(Y1))]
    return(Xdisp,Ydisp,Zdisp,X2,Y2,Z2)

