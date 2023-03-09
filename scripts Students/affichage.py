# affichage.py
"""========================================================================
Bibliothèque 
    de fonctions pour la création de fenêtre et de gaphique 
    à l'aide de matplotlib

--> figure 3D
    def initfigure3D(namex, namey, namez, title):

--> 2 figures 2D côte à côte
    def two_viewsfigure(fig1x, fig1y, titre1, fig2x, fig2y, titre2): 
        
--> figure 2D avec les limites données
    def initfigure2D(num, aa, bb, cc, dd, namex, namey, title):

--> 2 figures 3D côte à côte
    def two_viewsfigure3D(fig1x,fig1y,fig1z,titre1,fig2x,fig2y,fig2z,titre2):
        
--> 3 figures 2 2D, 1 3D côte à côte
    def three_viewsfigure(a,b,c,d,e,f):        

--> Affichage 2D (dans le plan (X,Y)) des points de devant triangulés 
    pour choisir les hauteurs ymin et ymax
    def DisplayTriangulation2D(xx2,yy2,ax):
========================================================================"""


import matplotlib.pyplot as plt

"""------------------------------------------------------
GESTION DE L'AFFICHAGE
------------------------------------------------------"""

def initfigure3D(namex, namey, namez, title):
    """
    initialise et affiche un espace 3d avec les noms donnés 
    axe X   puis   axe Y    puis axe Z  et enfin le titre du graphique
    Input : namex, namey, namez, title (noms des axes et le titre de la figure)
    Output : ax (l'axe), fig (la figure)
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel(namex)
    ax.set_ylabel(namey)
    ax.set_zlabel(namez)
    ax.set_title(title)
    return ax,fig


def two_viewsfigure(fig1x,fig1y,titre1,fig2x,fig2y,titre2):  
    """
    initialise et affiche 2 espaces 2d avec les noms donnes
    d'abord pour la premiere figure axe 1 , axe 2 puis titre
    et ensuite pour la deuxieme figure axe1 axe 2 puis titre
    Input : fig1x,fig1y,titre1,fig2x,fig2y,titre2 (noms des axes figure 1 +titre 
            puis figure 2 +titre )
    Output : (ax1,ax2),fig2 = les 2 axes et le pointeur vers la figure
    """
    fig2 = plt.figure(2,(12,5))
    #figure 1
    ax1 = fig2.add_subplot(121) 
    ax1.set_xlabel(fig1x)
    ax1.set_ylabel(fig1y)
    ax1.set_title(titre1) 
    #figure 2
    ax2 = fig2.add_subplot(122) 
    ax2.set_xlabel(fig2x)
    ax2.set_ylabel(fig2y)
    ax2.set_title(titre2) 
    return (ax1,ax2),fig2


def initfigure2D(num,aa,bb,cc,dd,namex,namey,title):
    """
    initialise et affiche un espace 2d avec les noms donnés 
    axe X   puis   axe Y   et  le titre du graphique
    et en ayant pour limite en x (aa bb)  et en y (cc dd)
    Input : num = numéro de la figure
            aa,bb,cc,dd,namex,namey,title (les limites des axes puis les noms 
            et le titre de la figure)
    Output : ax (l'axe), fig (la figure)
    """
    fig = plt.figure(num,(8,8))
    ax = fig.add_subplot()
    ax.set_title(title)
    ax.set_xlim((aa,bb))
    ax.set_ylim((cc,dd))
    ax.set_xlabel(namex)
    ax.set_ylabel(namey)
    return ax, fig


def two_viewsfigure3D(fig1x,fig1y,fig1z,titre1,fig2x,fig2y,fig2z,titre2):
    """
    initialise et affiche 2 espaces 3d avec les noms donnes
    d'abord pour la premiere figure axe 1 , axe 2, axe 3 puis titre
    et ensuite pour la deuxieme figure axe1 axe 2, axe 3 puis titre
    Input : fig1x,fig1y,fig1z,titre1,fig2x,fig2y,fig2z,titre2 (noms des axes figure 1 +titre 
            puis figure 2 +titre)
    Output : ax6,ax7 (les 2 figures)
    """
    fig = plt.figure(1,(12,6))
    #figure 1
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.set_title(titre1)
    ax1.set_xlabel(fig1x)
    ax1.set_ylabel(fig1y)
    ax1.set_zlabel(fig1z)
    #figure 2
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.set_title(titre2)
    ax2.set_xlabel(fig2x)
    ax2.set_ylabel(fig2y)
    ax2.set_zlabel(fig2z)
    return(ax1,ax2)


def three_viewsfigure(a,b,c,d,e,f):
    """
    initialise et affiche 2 espaces 2d et un espace 3d avec les limites données
    Input : a,b,c,d,e,f (les limites des différents axes)
    Output : ax1,ax2,ax3 (les 3 figures)
    """ 
    
    fig11 = plt.figure(11,figsize=(18,6))
    
    #première figure : affichage dans X,Y
    ax1 = fig11.add_subplot(131)
    ax1.set_xlabel('x axis')
    ax1.set_ylabel('y axis')
    ax1.set_title('X,Y axis')
    # optionel (deformations future) :
    ax1.axis([a, b, c, d])
    
    #deuxième figure : vue 3D
    ax2 = fig11.add_subplot(132,projection='3d')
    ax2.set_xlabel('x axis')
    ax2.set_ylabel('y axis')
    ax2.set_zlabel('z axis')
    ax2.set_title('B-spline surface')
     # optionel (deformations future) :
    ax2.set_xlim3d(a, b)
    ax2.set_ylim3d(c, d)
    ax2.set_zlim3d(e, f)
    #on aggrandit les axes de ax2
    pos1 = ax2.get_position()
    pos2 = [pos1.x0 - 0.045, pos1.y0 - 0.15,
            pos1.width * 1.35,
            pos1.height * 1.35]
    ax2.set_position(pos2)
    
    #troisième figure : affichage dans X,Z
    ax3 = fig11.add_subplot(133)
    ax3.set_xlabel('x axis')
    ax3.set_ylabel('z axis')
    ax3.set_title('X,Z axis')
    # optionel (deformations future) :
    ax3.axis([a, b, e, f])

    return(ax1,ax2,ax3)


# NEW
from scipy.spatial import Delaunay
import numpy as np

def DisplayTriangulation2D(xx2,yy2,ax):
    """
        Affichage 2D (dans le plan (X,Y)) des points de devant triangulés 
        pour choisir les hauteurs ymin et ymax
        input : xx2, yy2 = points 2D 
                ax = pointeur axe
    """
    # triangulation :
    points = np.column_stack((xx2,yy2))
    tri = Delaunay(points)
    # affichage des (X,Y) projetés pour choisir les hauteurs
    ax.triplot(xx2, yy2, tri.simplices.copy(),linewidth=0.2)
    plt.axis('equal')
    ax.plot(xx2,yy2,'or',markersize=0.4)
    plt.draw()
    # no return


