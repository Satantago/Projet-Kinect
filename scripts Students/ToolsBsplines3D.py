# ToolsBsplines3D.py

"""========================================================================
Bibliothèque
    de fonctions pour la reconstruction géométrique
    des données issues de la kinect :

--> Evaluation of a cubic Bspline 3D curve
    def Bspline3DCurveEvaluation(cpx,cpy,cpz,tk):

    Evaluation d'une surface B-spline : 
--> CAS PARAMETRIQUE
    def BsplineSurfaceParamEvaluation(CPx,CPy,CPz,knotsU,knotsV,Nbu,Nbv):

========================================================================"""
from ToolsBsplinesFunctions import *


# Evaluation of a cubic Bspline 3D curve
def Bspline3DCurveEvaluation(pcx,pcy,pcz,tk):
    """
        evaluation of a cubic Bspline 3D curve
        cpx,cpy,cpz : control points
        tk : sequence of knots
        Output : xs ,ys,zs (coordonnées des points de la cubic splines 3D)
    """
    nb_noeuds = len(tk)
    a = tk[0]
    b = tk[-1]
    # Bspline evaluation :
    nbval = 500
    t = np.linspace(a,b,nbval)
    xs = np.zeros(nbval)
    ys = np.zeros(nbval)
    zs = np.zeros(nbval)
    for k in range(nbval):
        tc = t[k]
        for i in range(nb_noeuds+2):
            xs[k] += pcx[i]*Bspline(tk,i,tc)
            ys[k] += pcy[i]*Bspline(tk,i,tc)
            zs[k] += pcz[i]*Bspline(tk,i,tc)
    return xs, ys, zs


# CAS PARAMETRIQUE
def BsplineSurfaceParamEvaluation(CPx,CPy,CPz,knotsU,knotsV,Nbu,Nbv):
    """ Evaluation of a tensor B-spline surface from its control net.
        Parametric case
        Input:
        CPx, CPy, CPy : 2D arrays of the coordinates of the control points
        knotsU, knotsV : vectors of u-knots and v-knots
        Nbu, Nbv : number of points for the u and v discretization 
        return Sx,Sy,Sz = 2D arrays of surface points
    """
    NbKnotsU = len(knotsU)
    NbKnotsV = len(knotsV)
    a = knotsU[0]
    b = knotsU[-1]
    c = knotsV[0]
    d = knotsV[-1]
    # discrétisation pour l'évaluation de la surface 
    ui = np.linspace(a,b,Nbu)
    vj = np.linspace(c,d,Nbv)
    
    Sx = np.zeros((Nbu,Nbv))
    Sy = np.zeros((Nbu,Nbv))
    Sz = np.zeros((Nbu,Nbv))
    
    Sku = np.zeros(Nbu)
    Slv = np.zeros(Nbv)
    for k in range(NbKnotsU + 2):
        for i in range (Nbu):
            Sku[i] = Bspline(knotsU,k,ui[i])    
        for l in range(NbKnotsV + 2):
            for j in range (Nbv):
                Slv[j] = Bspline(knotsV,l,vj[j])
            Suv = np.outer(Sku,Slv)
            Sx = Sx + CPx[k,l] * Suv
            Sy = Sy + CPy[k,l] * Suv
            Sz = Sz + CPz[k,l] * Suv
    return Sx,Sy,Sz 
