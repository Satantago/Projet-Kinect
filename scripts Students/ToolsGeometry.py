# ToolsGeometry.py

"""========================================================================
Bibliothèque (1) 
    de fonctions pour la reconstruction géométrique
    des données issues de la kinect :

--> Tracé d'un polyèdre de contrôle (Bézier ou spline)
    PERMET de passer en paramètre les tailles d'affichage
    def PolyhedronPlotV2(ax,CPx,CPy,CPz, tailleLW, tailleS):        
        
--> Tracé d'un polyèdre de contrôle (Bézier ou spline) en 3D
    no scatter 
    def PolyhedronPlot3D(ax,CPx,CPy,CPz, tailleLW, tailleS):
    
--> Tracé d'un polyèdre de contrôle (Bézier ou spline) en 2D
    no scatter
    def PolyhedronPlot2D(ax,CPu,CPv, tailleLW, tailleS):

--> Choix et tracé à la souris des hauteurs ymin et ymax 
    dans une fenêtre graphique 2D
    def ChoixHauteurs(ax,xmin,xmax):
        
--> Choix et tracé à la souris d'un vecteur dans une fenêtre graphique 2D
    def ChoixVecteur(ax):
        
--> Choix des angles thetaMin et thetaMax par 2 points définissant 
    les extrémités de 2 vecteurs partant de l'origine
    def ChoixAngles(ax):

--> Determine the smallest gap between 2 consecutive points
    in matrix P, either horizontally, or vertically
    def smallestGap(P):

--> Determine the smallest gap between 2 consecutive points
    in 3 matrices, either horizontally, or vertically
    def eps(PCX,PCY,PCZ):

========================================================================"""
import numpy as np
import matplotlib.pyplot as plt


# Tracé d'un polyèdre de contrôle (Bézier ou spline)
# PERMET de passer en paramètre les tailles d'affichage
def PolyhedronPlotV2(ax,CPx,CPy,CPz, tailleLW, tailleS):
    """ Plotting of the control polyhedron of the surface
        defined by control points Cpx,Cpy,Cpz
        in the window with subplot "ax"
        Input:
        CPx, CPy, CPz : 2D arrays of the coordinates of the control points
        tailleLW et tailleS sont des paramètres d'affichage
    """
    (n0,m0) = np.shape(CPx)
    # n0 = udegree + 1
    # m0 = vdegree + 1
    # Plot of the Bézier control points :
    ax.scatter(CPx,CPy,CPz, c='c', marker='o', s=tailleS)
    # Plot of the Bézier control polyhedron :
    for i in range(n0):
            sx = CPx[i,:]
            sy = CPy[i,:]
            sz = CPz[i,:]
            ax.plot(sx,sy,sz,'-c',lw=tailleLW)
    for j in range(m0):
            sx = CPx[:,j]
            sy = CPy[:,j]
            sz = CPz[:,j]
            ax.plot(sx,sy,sz,'-c',lw=tailleLW)


# Tracé d'un polyèdre de contrôle (Bézier ou spline) en 3D
# no scatter 
def PolyhedronPlot3D(ax,CPx,CPy,CPz, tailleLW, tailleS):
    """ Plotting of the control polyhedron of the surface
        defined by control points Cpx,Cpy,Cpz
        in the window with subplot "ax"
        Input:
        CPx, CPy, CPz : 2D arrays of the coordinates of the control points
        tailleLW et tailleS sont des paramètres d'affichage
    """
    (n0,m0) = np.shape(CPx)
    # Plot of the Bézier control points :
    # ax.scatter(CPx,CPy,CPz, c='c', marker='o', s=tailleS)
    # Plot of the Bézier control polyhedron :
    for i in range(n0):
            sx = CPx[i,:]
            sy = CPy[i,:]
            sz = CPz[i,:]
            ax.plot(sx,sy,sz,'c-',lw=tailleLW)
            ax.plot(sx,sy,sz,'oc')
    for j in range(m0):
            sx = CPx[:,j]
            sy = CPy[:,j]
            sz = CPz[:,j]
            ax.plot(sx,sy,sz,'-c',lw=tailleLW)


# Tracé d'un polyèdre de contrôle (Bézier ou spline) en 2D
# no scatter
def PolyhedronPlot2D(ax,CPu,CPv, tailleLW, tailleS):
    """ Plotting the control polyhedron of the surface in 2D
        defined by control points CPu, CPv
        in the window with subplot "ax"
        Input:
        CPu, CPv : 2D arrays of the coordinates of the control points
        tailleLW et tailleS sont des paramètres d'affichage
    """
    (n0,m0) = np.shape(CPu)
    # Plot of the Bézier control polyhedron :
    for i in range(n0):     # Points
        su = CPu[i,:]
        sv = CPv[i,:]
        ax.plot(su,sv,'oc')
            
    for i in range(n0):     # Segments
        su = CPu[i,:]
        sv = CPv[i,:]
        ax.plot(su,sv,'-c',lw=tailleLW)
            
    for j in range(m0):     # Segments
        su = CPu[:,j]
        sv = CPv[:,j]
        ax.plot(su,sv,'-c',lw=tailleLW)
            

# Choix et tracé à la souris des hauteurs ymin et ymax 
# dans une fenêtre graphique 2D
def ChoixHauteurs(ax,xmin,xmax):
    """
        Choix et tracé des hauteurs ymin et ymax dans une fenêtre graphique 2D
        (commencer par le bas)
        Input : ax = pointeur vers la fenêtre graphique
                xmin,xmax = valeurs min et max pour le tracé des hauteurs
        Output : ymin,ymax = hauteurs min et max
    """
    yu = []
    coord = 0
    while len(yu) < 2 :
        coord = plt.ginput(1, mouse_add=1, mouse_stop=3, mouse_pop=2)
        # coord is a list of tuples : coord = [(xx,yy)]
        # if right click ==> coord = []
        if coord != []:
            xx = coord[0][0]
            yy = coord[0][1]
            ax.plot(xx,yy,'ob',markersize=10)
            yu.append(yy)
            ax.plot([xmin,xmax],[yy,yy],'b--',lw=2)
            plt.draw()
    ymin = np.array(yu[0]) 
    ymax = np.array(yu[1])
    return ymin,ymax


# Choix et tracé à la souris d'un vecteur dans une fenêtre graphique 2D
def ChoixVecteur(ax):
    """
        Choix et tracé à la souris d'un vecteur dans une fenêtre graphique 2D
        afin d'effectuer ensuite une rotation
        (pour un visage, commencer par le menton)
        Input : ax = pointeur vers la fenêtre graphique
        Output : xu,yu = array représentant les coordonnées 
                 des extrémités de ce vecteur
    """
    xu = []  
    yu = []
    coord = 0
    while len(xu) < 2 :
        coord = plt.ginput(1, mouse_add=1, mouse_stop=3, mouse_pop=2)
        # coord is a list of tuples : coord = [(xx,yy)]
        # if right click ==> coord = []
        if coord != []:
            xx = coord[0][0]
            yy = coord[0][1]
            ax.plot(xx,yy,'or',markersize=6)
            xu.append(xx)
            yu.append(yy)
            plt.draw()
            if len(xu) == 2 :
                ax.plot([xu[-2],xu[-1]],[yu[-2],yu[-1]],'r--',lw=2)
    
    xu = np.array(xu) 
    yu = np.array(yu)
    return xu,yu


# Choix des angles thetaMin et thetaMax par 2 points définissant 
# les extrémités de 2 vecteurs partant de l'origine
def ChoixAngles(ax):
    """
        Choix et tracé à la souris de 2 points représentant les extrémités 
        de 2 vecteurs partant de l'origine dans une fenêtre graphique 2D
        afin de déterminer les angles Min et Max d'acquisition
        Input : ax = pointeur vers la fenêtre graphique
        Output : thetaMin,thetaMax = angles min et max 
    """
    xu = []  
    yu = []
    coord = 0
    while len(xu) < 2 :
        coord = plt.ginput(1, mouse_add=1, mouse_stop=3, mouse_pop=2)
        # coord is a list of tuples : coord = [(xx,yy)]
        # if right click ==> coord = []
        if coord != []:
            xx = coord[0][0]
            yy = coord[0][1]
            ax.plot(xx,yy,'or',markersize=6)
            xu.append(xx)
            yu.append(yy)
            ax.plot([0.,xx],[0.,yy],'r--',lw=2)
            plt.draw()
    xu = np.array(xu) 
    yu = np.array(yu)    
    thetaMin = min(np.arctan(xu[0]/yu[0]), np.arctan(xu[1]/yu[1]))
    thetaMax = max(np.arctan(xu[0]/yu[0]), np.arctan(xu[1]/yu[1]))
    print("thetaMin = ",thetaMin*180/np.pi," thetaMax = ",
          thetaMax*180/np.pi," en degrés")
    return thetaMin,thetaMax


def smallestGap(P):
    """
    Determine the smallest gap between 2 consecutive points
    in matrix P, either horizontally, or vertically
    Input : P une matrice 
    Output : une valeur  qui est le minimum de 2 points consécutifs
    """
    eH = P[:,1:] - P[:,:-1]
    minH = min(abs(eH.flatten()))
    eV = P[1:,:] - P[:-1,:]
    minV = min(abs(eV.flatten()))
    return min(minH,minV)



def eps(PCX,PCY,PCZ):
    """
    Determine the smallest gap between 2 consecutive points
    in 3 matrices, either horizontally, or vertically
    Input : 3 matrices
    Output : valeur min (de 2 points consécutifs) entre 3 matrices 
    """
    eX = smallestGap(PCX)
    eY = smallestGap(PCY)
    eZ = smallestGap(PCZ)
    return min(eX,eY,eZ)

