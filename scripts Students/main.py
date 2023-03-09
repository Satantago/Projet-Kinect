
from ToolsGeometry import *
from ToolsBorderCurves import *

from preparation import *
from selection import *
from reconstruction import *
from deformation import *
from affichage import *


"""======================================================================
A) MISE EN PLACE
========================================================================="""

"""----------------------------
A1) Chargement des données
-------------------------------"""
# création d'une figure 3D vide
ax,fig1 = initfigure3D("x axis","y axis","z axis","Initial data")
ax.get_figure()

# chargement et placement du buste 
# attention à marquer le bon chemin en plus du nom du fichier sans ".ply"
FilePly = "../Final_Face"
(X0,Y0,Z0,NbPts) = acquisition(FilePly)
# affichage du buste avec barycentre à l'origine
ax.scatter(X0,Y0,Z0, c='c', marker='o', s=0.1)   

#---------------------------------
# Fin Chargement des données
#---------------------------------

"""-------------------------------
A2) ROTATIONS
-------------------------------"""
# Sélection des points de devant pour choisir les rotations
(xA,yA,zA) = pointsdevant(X0, Y0, Z0)

# Fenêtres graphiques 2D pour sélectionner les rotations
(ax1,ax2),fig2 = two_viewsfigure('y-axis','z-axis',"affichage (y,z)",'x-axis','z-axis',"affichage (x,z)") 
ax1.plot(yA,zA,'co',markersize=0.2)
ax2.plot(xA,zA,'co',markersize=0.2)

"""
ROTATION X
"""
# choix du vecteur "d'horizontalité souhaitée" (de gauche à droite) 
# à la souris dans la figure 2D (y,z)
xu,yu = ChoixVecteur(ax1)
(Xdisp,Ydisp,Zdisp,X1,Y1,Z1) = rotationX(xu,yu,xA,yA,zA,X0,Y0,Z0)
ax1.plot(Ydisp,Zdisp,'bo',markersize=0.2)

"""
ROTATION Y
"""
# choix du vecteur de "verticalité souhaitée" à la souris 
# dans la figure 2D (x,z) : 
# origine = point (0,0)
# extrémité = nez
xu,yu = ChoixVecteur(ax2)
(Xdisp,Ydisp,Zdisp,X2,Y2,Z2) = rotationY(xu,yu,xA,yA,zA,X1,Y1,Z1)
ax2.plot(Xdisp,Zdisp,'bo',markersize=0.2)

#---------------------------------
# Fin des ROTATIONS
#---------------------------------

x2 = np.array(X2)
y2 = np.array(Y2)
z2 = np.array(Z2)
plt.close(fig1) # fermeture buste initial 3D

#--------------------------------------------------------------
# FIN MISE EN PLACE
#--------------------------------------------------------------


"""======================================================================
B) SELECTION des points pour la reconstruction : coordonnées cylindriques
   --> entre ym et yM : hauteurs min et max de la partie cylindrique
   --> et entre thetaMin et thetaMax :
    angle autour de l'axe y
    l'angle nul est dirigé vers les z croissants
    donc theta < 0 correspond à x < 0
========================================================================="""

"""-------------------------------------
B1) Sélection des hauteurs ymin et ymax
----------------------------------------"""
# sélection des points de devant (pour le choix de la partie à reconstruire)
(xx2,yy2,zz2) = pointsdevant(x2, y2, z2)
aa = min(xx2)
bb = max(xx2)
cc = min(yy2)
dd = max(yy2)

# Affichage 2D des points de "devant" (triangulés)
ax5,fig5 = initfigure2D(3,aa,bb,cc,dd,"x-axis","y-axis","choix ymin et ymax pour la reconstruction")
DisplayTriangulation2D(xx2,yy2,ax5)

# Sélection à la souris des hauteurs min et max (commencer par le bas)
ym, yM = ChoixHauteurs(ax5,aa,bb)
#-----------------------------------------------------
# Fin selection des hauteurs ymin et ymax
#-----------------------------------------------------

"""------------------------------------------
B2) Sélection des angles thetaMin et thetaMax
---------------------------------------------"""
# affichage des points de devant : vue de dessous => axe x vers la droite
(x2angle,y2angle,z2angle) = pointsdevant(x2, y2, z2)

ax5bis,fig5bis = initfigure2D(4,min(x2angle),max(x2angle),min(z2angle),max(z2angle),"x-axis","z-axis","choix theta")
ax5bis.plot(x2angle,z2angle,'or',markersize=0.4)
plt.axis("equal")

# choix des angles thetaMin et thetaMax
# commencer par un point à gauche (correspond à des x < 0)
thetaMin,thetaMax = ChoixAngles(ax5bis)
#-----------------------------------------------------
# Fin selection des angles thetaMin et thetaMax
#-----------------------------------------------------

plt.close(fig2) # fermeture fenetre des rotations


"""======================================================================
C) DETERMINATION et SEPARATION des points pour la reconstruction
   entre (ym et yM)  et entre (thetaMin et thetaMax)
    x3,y3,z3 : points conservés pour la reconstruction
    x3c,y3c,z3c : points complémentaires pour affichage (facultatif)
                  mais surtout pour le retour dans MeshLab
========================================================================="""

(x3,y3,z3,x3c,y3c,z3c,t3,th3,t3c) = anglepoints(x2,y2,z2,thetaMin,thetaMax,ym,yM)
# graphique avec 2 axes pour afficher la zone sélectionnée et son complémentaire
(ax6,ax7) = two_viewsfigure3D('X','Y','Z',"selected face",'X','Y','Z',"complementary part")

# affichage des points selectionnés selon 2 couleurs
ax6.scatter(x3, y3, z3, color='r', marker='o', s=0.3)
ax7.scatter(x3c, y3c, z3c, color='g', marker='o', s=0.1)
#affichage de la zone selectionné sur la tête complète
ax7.scatter(x3, y3, z3, color='r', marker='o', s=0.1)

#-----------------------------------------------------
#FIN DETERMINATION des points pour la reconstruction
#-----------------------------------------------------
plt.close(fig5)     # fermeture sélection hauteur
plt.close(fig5bis)  # fermeture sélection angles



"""======================================================================
D) DEBUT DE L'APPROXIMATION CYLINDRIQUE par B-spline paramétrique
========================================================================="""

"""--------------------------------------------
D1) CHOICE of B-spline knots (Noeuds uniformes)
-----------------------------------------------"""
nb_noeudsT = 12   # nombre de noeuds angulaire (T for Theta)
nb_noeudsY = 12   # nombre de noeuds selon y (verticaux)

# vecteurs des noeuds
VnoeudsT = np.linspace(thetaMin,thetaMax,nb_noeudsT)
VnoeudsY = np.linspace(ym,yM,nb_noeudsY)
#-----------------------------------------------------
# Fin choix noeuds B-spline
#-----------------------------------------------------

"""---------------------------------------------------
D2) DETECTION ET APPROXIMATION DES 4 COURBES FRONTIERES
------------------------------------------------------"""
err1 = 0.05         # erreur angulaire
err2 = 0.005        # erreur verticale
disp = True         # variable Booléenne pour l'affichage
nbLR = nb_noeudsY   # nombre de noeuds pour les courbes droite et gauches
nbTB = nb_noeudsT   # nombre de noeuds pour les courbes haute et basse
LC,RC,TC,BC,ax0 = BorderCurvesM4(x3,y3,z3,thetaMin,thetaMax,
                                 ym,yM,err1,err2,nbLR,nbTB,disp)

#-----------------------------------------------------
# Fin Detection et Approximation des 4 courbes frontières
#-----------------------------------------------------

"""-----------------------------------------------------
D3) approximation par moindres carrés sans contraintes:
   construction du systeme principal : equations normales
--------------------------------------------------------"""
(ATA,ATx,ATy,ATz,nbPC) = moindres_carres (x3,y3,z3,nb_noeudsT,nb_noeudsY, th3,VnoeudsT,VnoeudsY)

#"""------------------------------------------------------------
#FIN approximation par moindre carré sans contraintes
#---------------------------------------------------------------"""

"""--------------------------------------------------------
D4) COUTURE : approximation moindres carres avec contraintes 
   Construction du systeme MC avec des contraintes linéaires
-----------------------------------------------------------"""
# Construction du système avec contraintes linéaire
H, Bx, By, Bz = ConstraintsSystem(nb_noeudsT,nb_noeudsY,LC,RC,TC,BC)

(pcxx,pcyy,pczz) = moindres_carres_contrainte(ATA, ATx, ATy, ATz, Bx, By, Bz, H, nb_noeudsT, nb_noeudsY)


#calcul des points de controle
(PCX,PCY,PCZ) = point_de_controle(pcxx, pcyy, pczz,nbPC, nb_noeudsY, nb_noeudsT)

#---------------------------------------------------------------
# FIN DE L'APPROXIMATION CYLINDRIQUE par B-spline paramétrique
#---------------------------------------------------------------


"""======================================================================
E) AFFICHAGES
========================================================================="""

# affichage de l'approximation paramétrique
ax,fig5 = initfigure3D("angle axis","y axis","radius axis","parametric approximation")

nbTEval = 50
nbYEval = 40
# approximation cylindrique : r = r(theta,y)
Sx, Sy, Sz = BsplineSurfaceParamEvaluation(PCX, PCY, PCZ, VnoeudsT, VnoeudsY, nbTEval, nbYEval)

ax.plot_wireframe(Sx, Sy, Sz, color='b', linewidth=0.5)


# AFFICHAGE de la reconstruction B-spline AVEC LES COURBES DE BORDURE :
ax0.plot_wireframe(Sx, Sy, Sz, color='b', linewidth=0.5)

plt.show()

# affichage des points de contrôles:
#PolyhedronPlotV2(ax,PCX,PCY,PCZ, 1, 25)

# les points (xi,yi,zi) approximés par la Bspline cylindrique :
ax.scatter(x3, y3, z3, color='r', marker='o', s=0.2)
# et la partie complémentaire : (Facultatif : on voit plus rien)
# ax.scatter(x3c, y3c, z3c, color='g', marker='o', s=0.1)

#"""------------------------------------------------------------
#FIN AFFICHAGE
#------------------------------------------------------------"""


"""======================================================================
F) DEFORMATION
========================================================================="""

"""---------------------------------
F1) AFFICHAGES POUR DEFORMATIONS
---------------------------------"""
# bounding box de données
mx = min(PCX.flatten()) ; Mx = max(PCX.flatten()) ; hx = Mx - mx
my = min(PCY.flatten()) ; My = max(PCY.flatten()) ; hy = My - my
mz = min(PCZ.flatten()) ; Mz = max(PCZ.flatten()) ; hz = Mz - mz

# Un Graphique avec 3 axes
(ax1,ax2,ax3) = three_viewsfigure(mx-hx/10,Mx+hx/10, my-hy/10,  My+hy/10, mz-hz/10, Mz+1.5*hz)

# affichage dans ax2 : vue 3D
# polyhedron de controle
PolyhedronPlot3D(ax2,PCX,PCY,PCZ, 2, 50)

# Construction de la surface Bspline
Sx, Sy, Sz = BsplineSurfaceParamEvaluation(PCX,PCY,PCZ,VnoeudsT,VnoeudsY,nbTEval,nbYEval)
# affichage :
ax2.plot_wireframe(Sx, Sy, Sz, color='b', linewidth=0.6)

# dans ax1 : X,Y
PolyhedronPlot2D(ax1,PCX,PCY, 1, 50)

# dans ax3 : X,Z
PolyhedronPlot2D(ax3,PCX,PCZ, 1, 50)  

#"""------------------------------------------------------------
#FIN AFFICHAGE
#------------------------------------------------------------"""

"""---------------------------------
F2) INITIALISATION VARIABLES
---------------------------------"""
#variables globale pour la déformation définies dans le fichier deformation.py
shape(PCX)

#"""-------------------------------------------------
# FIN  Initialisation variables
#-------------------------------------------------"""

"""---------------------------------
F3) DEFORMATION INTERACTIVE
---------------------------------"""
# Nous allons bouger les points de controles, soit dans l'axe ax1 soit dans
# ax3 à l'aide des fonctions d'événements à la souris
# en se servant de l'axe 2 pour observer les changements en temps réels
 
# Regroupement de l'ensemble des fonctions d'interaction :
(PCX,PCY,PCZ) = RealTime(num, but ,ii ,jj, n0,m0, epsilon,ax1,ax2,ax3,PCX,PCY,PCZ,VnoeudsT,VnoeudsY,nbTEval,nbYEval)

#"""------------------------------------------------------------
#FIN DEFORMATION intéractive
#---------------------------------------------------------------"""

"""======================================================================
G) RETOUR dans MeshLab
========================================================================="""
#- pour la partie complémentaire (points verts)
#- pour la partie reconstruite (à partir des points rouges
#  modifié en bleu par la déformation intéractive) => B-spline 

# Optionel : pour une meilleure résolution, on ré-évalue la B-spline
nbTEval = 70
nbYEval = 60

# Reconstruction complete du Mesh (format PLY) après l'approximation 3D :
    
# 1) on enleve les sommets avec l'index t3 dans les données initiales
#    on met à jour le visage
# 2) on crée le Mesh ply associé avec la surface B-spline
#    mais d'abord on met à jour les 3 matrices Sx,Sy,Sz
# RMQ : le fichier final Mesh sera appélé: FilePly + '4.ply'
Sx, Sy, Sz = BsplineSurfaceParamEvaluation(PCX,PCY,PCZ,VnoeudsT,VnoeudsY,nbTEval,nbYEval)

WholeMeshToPly(FilePly,x2,y2,z2,t3,Sx, Sy, Sz)

#"""------------------------------------------------------------
# Fin RETOUR dans MeshLab
#------------------------------------------------------------"""






