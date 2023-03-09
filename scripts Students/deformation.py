# deformation.py

"""========================================================================
Bibliothèque 
    deformation active du visage

--> mise à jour des variables globales n0 et m0
    def shape(PCX):

--> cherche un point à proximité du clic
    def LookForPointToMove(CPu,CPv,xm,ym,m0,epsilon):
        
--> gère l'entièreté des fonctions de la déformation interactive
    def RealTime(num, but ,ii ,jj, n0, m0, epsilon, ax1, ax2, ax3, PCX,PCY,PCZ,VnoeudsT,VnoeudsY,nbTEval,nbYEval):

        --> Cherche un point à proximité du clic détécté
            def onPress(event):
                
        --> Déplace le point choisi et recalcule les points associés
            def onMove(event):
                
        --> Met à jour les variables de la souris quand les bouttons sont relachés
            def onReleased(event)::
========================================================================"""

from ToolsBsplines3D import *

"""-------------------------------------------------
 Initialisation variables
-------------------------------------------------"""

epsilon = 0.001   # distance au point choisi
num = -1        # numéro ou index du point en mouvement
but = 0         # bouton souris (1 : gauche, 3 : droit)
ii = -1         # (ii,jj) = indices du point choisi en mouvement
jj = -1
(n0,m0) = (-1,-1) # initialisation des variables globales n0 et m0 par des valeurs absurdes

def shape(PCX):
    """
    met à jour les variables globales n0 et m0
    """
    global n0, m0
    (n0,m0) = np.shape(PCX)
    return


def LookForPointToMove(CPu,CPv,xm,ym,m0,epsilon):
    """
    On cherche si le clic de la souris effectué est à proximité d'un point de
    controle
    """
    number = -1
    k = 0
    while (k < m0 and number < 0):
        # On cherche si un point est à distance "epsilon" de la souris :
        if (abs(CPu[k] - xm) < epsilon) and (abs(CPv[k] - ym) < epsilon):
            number = k
        else:
            k = k+1
    return number


def RealTime(num, but ,ii ,jj, n0, m0, epsilon, ax1, ax2, ax3, PCX,PCY,PCZ,VnoeudsT,VnoeudsY,nbTEval,nbYEval):
    """
    gère l'entièreté des fonctions de la déformation interactive
    """
    def onPress(event):
        """
        cherche si à proximité du clic gauche de la souris effectué il y a un point à déplacer
        """
        global num, but, ii, jj,n0, m0, epsilon
        
        if event.inaxes:
            ax = event.inaxes   
            liste = ax.get_children()
            # on recherche un point(avec pour nombre 'num') à bouger:
            for i in range(n0):
                CPu = liste[i].get_xdata()
                CPv = liste[i].get_ydata()
                num = LookForPointToMove(CPu, CPv, event.xdata, event.ydata,m0,epsilon)
                # si on en trouve un :
                if num >= 0:
                    but = 1
                    ii = i
                    jj = num 
                    break   # pour sortir de la boucle for    
        if event.button == 3:
            print('Stop Moving points')
            plt.disconnect(cidPress)
            plt.disconnect(cidMove)
            plt.disconnect(cidReleased)
        
    
    def onMove(event):
        """
        définie les actions de déplacemnts d'un point et recalcule les nouvelles splines associées'
        """
        global num, but, ii, jj, n0, m0
        # la souris bouge avec "clique gauche" 
        if ((event.inaxes == ax1) or (event.inaxes == ax3)) and but == 1:
            if event.inaxes == ax1 and but == 1:
                axW = ax1   # fenetre de travail
                axS = ax3   # fenetre subsidiaire 
            if event.inaxes == ax3 and but == 1:
                axW = ax3   # fenetre de travail
                axS = ax1   # fenetre subsidiaire 
            # Stockage des informations des differentes fenetres graphiques 
            # dans des listes : travail / subsidiaire
            listeW = axW.get_children()
            listeS = axS.get_children()
            
            #--------------------------
            # POINTS DE CONTROLE
            # on récupère les points de contrôle actuel des axes:
            pcx = listeW[ii].get_xdata()
            pcy = listeW[ii].get_ydata()
            # on met à jour les coordonnées du point en mouvement:
            pcx[jj] = event.xdata
            pcy[jj] = event.ydata
            # on met à jour les points de contrôle dans les axes principaux:
            listeW[ii].set_xdata(pcx)
            listeW[ii].set_ydata(pcy) 
            plt.draw()
            # on met à jour dans les axes secondaires:
            listeS[ii].set_xdata(pcx)
            plt.draw()
            
            #--------------------------
            # POLYGONES DE CONTROLE
            # On met à jour le polygône de controle dans la figure:
            # lignes horizontales dans axW
            listeW[ii+n0].set_xdata(pcx)
            listeW[ii+n0].set_ydata(pcy)
            # lignes horizontales dans axS
            listeS[ii+n0].set_xdata(pcx)
            # lignes verticales dans axW
            vlx = listeW[n0+n0+jj].get_xdata()
            vly = listeW[n0+n0+jj].get_ydata()
            vlx[ii] = event.xdata
            vly[ii] = event.ydata
            listeW[n0+n0+jj].set_xdata(vlx)
            listeW[n0+n0+jj].set_ydata(vly)
            # lignes verticales dans axS           
            listeS[n0+n0+jj].set_xdata(vlx)  
            plt.draw()
            
            #--------------------------
            # LA SURFACE
            # Update des matrices des points de controle en x,y,z 
            Sx, Sy, Sz = BsplineSurfaceParamEvaluation(PCX,PCY,PCZ,
                                                       VnoeudsT,VnoeudsY,
                                                       nbTEval,nbYEval)
            # on met à jour la surface dans l'axe 3D ax2
            liste2 = ax2.get_children()
            Nb_u = nbTEval
            Nb_v = nbYEval
            Line3D = []
            for i in range(Nb_u):
                L3D = []
                for j in range(Nb_v):
                    L3D.append((Sx[i,j],Sy[i,j],Sz[i,j]))
                Line3D.append(L3D)
    
            for j in range(Nb_v):
                L3D = []
                for i in range(Nb_u):
                    L3D.append((Sx[i,j],Sy[i,j],Sz[i,j]))
                Line3D.append(L3D)
    
            liste2[0].set_segments(Line3D) # on injecte les nouvelles coordonnées
            plt.draw()
                
    def onReleased(event):
        """
        Met à jour les variables de la souris quand les bouttons sont relachés
        """
        global num, but
        num = -1
        but = 0
        
    cidMove= plt.connect('motion_notify_event', onMove)
    cidPress  = plt.connect('button_press_event', onPress)
    cidReleased = plt.connect('button_release_event', onReleased)
    return (PCX,PCY,PCZ)
