# reconstruction.py
import random

"""========================================================================
Bibliothèque 
    de fonctions pour la reconstruction géométrique
    des données issues de la kinect :

--> Moindres carrées sans contraintes
    def moindres_carres (x3,y3,z3,nb_noeudsT,nb_noeudsY,th3,VnoeudsT,VnoeudsY):

--> Moindres carrées avec contraintes
    def moindres_carres_contrainte(ATA,ATx,ATy,ATz,Bx,By,Bz,H,nb_noeudsT,nb_noeudsY):
        
--> Réseau des points de contrôle
    def point_de_controle(pcxx,pcyy,pczz,nbPC,nb_noeudsY,nb_noeudsT):

========================================================================"""

from ToolsBsplines3D import *


def moindres_carres (x3,y3,z3,nb_noeudsT,nb_noeudsY,th3,VnoeudsT,VnoeudsY):
    """
    approximation moindres carrée sans contraintes:
    Construction du systeme principal : "normal equations"
    Input : x3,y3,z3 (les points sur lesquels on fait les calculs)
            nb_noeudsT,nb_noeudsY (nombre de noeuds)
            th3 (angles associés aux points)
            VnoeudsT,VnoeudsY (vecteur des noeuds)
    Output : ATA,ATx,ATy,ATz (matrices des equations normales)
             nbC (nombre de points de contrôle)
    """
    # matrice moindre carré
    nbData3 = len(y3)
    nbPC    = (nb_noeudsT+2)*(nb_noeudsY+2)
    A = np.zeros((nbData3,nbPC))
    # construction de la matrice A pour les moindres carrés
    
    for i in range(nbData3):
        for j in range(nbPC):
            A[i, j] = Bspline(VnoeudsT, j%(nb_noeudsT+2), th3[i])*Bspline(VnoeudsY, j//(nb_noeudsT+2), y3[i])
    
    AT = A.transpose()
    
    
    
    x3 = np.array(x3)
    y3 = np.array(y3)
    z3 = np.array(z3)
    
    # Equations normales :
    ATA = np.dot(AT, A)
    ATx = np.dot(AT, x3)
    ATy = np.dot(AT, y3)
    ATz = np.dot(AT, z3)
    
    #Résolution directe des équations normales ( avec linalg.solve)
    pcx = np.linalg.solve(ATA, ATx)
    pcy = np.linalg.solve(ATA, ATy)
    pcz = np.linalg.solve(ATA, ATz)
    return(ATA,ATx,ATy,ATz,nbPC)





def moindres_carres_contrainte(ATA,ATx,ATy,ATz,Bx,By,Bz,H,nb_noeudsT,nb_noeudsY):
    """
    approximation moindres carrés avec contraintes:
    Construction de la matrice des moindres carrés avec
    des contraintes linéaires
    Input : ATA,ATx,ATy,ATz (equations normales principales)
            Bx,By,Bz,H (matrices des contraintes)
            nb_noeudsT,nb_noeudsY (nombre de noeuds)
    Output : pcxx,pcyy,pczz (les points de controle, mais non correctement structurés...)
    """
    # matrice principale AHH
    AHH = np.zeros((len(ATA) + len(H), len(ATA[0]) + len(H)))
    for i in range(len(ATA)):
        for j in range(len(ATA[0])):
            AHH[i, j] = ATA[i, j]
    for i in range(len(ATA), len(ATA) + len(H)):
        for j in range(len(ATA[0])):
            AHH[i, j] = H[i-len(ATA), j]
    for i in range(len(ATA)):
        for j in range(len(ATA[0]), len(ATA[0]) + len(H)):
            AHH[i, j] = H[j-len(ATA[0]), i]
    
    
    # deuxieme membre
    ATBx = np.array([x for x in ATx] + [x for x in Bx])
    ATBy = np.array([y for y in ATy] + [y for y in By])
    ATBz = np.array([z for z in ATz] + [z for z in Bz])
    
    # Resolution directe (with linalg.solve)
    pcxx = np.linalg.solve(AHH, ATBx)
    pcyy = np.linalg.solve(AHH, ATBy)
    pczz = np.linalg.solve(AHH, ATBz)
    return(pcxx,pcyy,pczz)


def point_de_controle(pcxx,pcyy,pczz,nbPC,nb_noeudsY,nb_noeudsT):
    """
    Fabrication du Réseau des points de contrôle
    Input : pcxx,pcyy,pczz (les points)
            nbPC (nombre points de controle)
            nb_noeudsY,nb_noeudsT (nombres de noeuds)
    Output : PCX,PCY,PCZ (les points de controle structurés)
    """
    # on récupère des points de contrôles
    pcx = pcxx[:nbPC]
    pcy = pcyy[:nbPC]
    pcz = pczz[:nbPC]
    
    PCX = pcx.reshape(nb_noeudsY+2,nb_noeudsT+2)
    PCY = pcy.reshape(nb_noeudsY+2,nb_noeudsT+2)
    PCZ = pcz.reshape(nb_noeudsY+2,nb_noeudsT+2)
    PCX = PCX.T
    PCY = PCY.T
    PCZ = PCZ.T
    return(PCX,PCY,PCZ)
