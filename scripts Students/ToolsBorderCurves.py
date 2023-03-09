# ToolsBorderCurves.py

"""========================================================================
Bibliothèque
    Detection and B-spline approximation of the 4 border curves 
    of the selected area

--> B-spline approximation of 2D-data
    def curve_approximation(xx,yy,zz,tt,nb_noeuds):

--> Detection and B-spline approximation of the border curves (method M4)
    BE CAREFUL : corner points are not necessary identical
    def BorderCurvesM4(x2,y2,z2,thetaMin,thetaMax,ym,yM,
                       err1,err2,nbLR,nbTB,disp):

--> Definition of the linear system of constraints 
    associated with the border curves of the B-spline surfaces
    def ConstraintsSystem(nb_noeudsT,nb_noeudsY,LC,RC,TC,BC):
========================================================================"""
from ToolsBsplines3D import *

#-------------------------------------------------------------------
# B-spline approximation of 2D-data (in 3D space)
def curve_approximation(xx,yy,zz,tt,a,b,nb_noeuds):
    """
    Approximation of data (xx,yy,zz) according to parameters tt
    on the interval [a,b]
    by a uniform cubic Bspline curve with "nb_noeuds" knots 
    Output : pcx, pcy, pcz : Bspline control points coordinates
             xs, ys, zs : sampling of the Bspline solution
    """
    N = len(xx)
    # sequence of uniform Bspline knots :
    tk = np.linspace(a,b,nb_noeuds)

    # least squares approximation matrix :
    A = np.zeros((N,nb_noeuds+2))
    for i in range(N):
        for j in range(nb_noeuds+2):
            A[i, j] = Bspline(tk, j, tt[i])
    AT = A.transpose()
    
    xi = np.array(xx)
    yi = np.array(yy)
    zi = np.array(zz)
    
    # Normal equations :
    ATA = np.dot(AT, A)
    ATX = np.dot(AT, xi)
    ATY = np.dot(AT, yi)
    ATZ = np.dot(AT, zi)

    # resolution with numpy :
    pcx =  np.linalg.solve(ATA, ATX)    # Bspline x-control points
    pcy =  np.linalg.solve(ATA, ATY)    # Bspline y-control points
    pcz =  np.linalg.solve(ATA, ATZ)    # Bspline z-control points
    
    xs, ys, zs = Bspline3DCurveEvaluation(pcx,pcy,pcz,tk)
    return pcx, pcy, pcz, xs, ys, zs


#-------------------------------------------------------------------
# Detection and B-spline approximation of the border curves (method M4)
# BE CAREFUL : corner points are not necessary identical
def BorderCurvesM4(x2,y2,z2,thetaMin,thetaMax,ym,yM,err1,err2,nbLR,nbTB,disp):
    """
        Detection and B-spline approximation of the border curves 
        of the selected area
        --> for method M4
        (x2,y2,z2) : coordinates of all points
        thetaMin,thetaMax,ym,yM : parameters of the selected area
        err1 : angular error
        err2 : vertical error (according y-coordinates)
        nbLR : number of knots for Left & Right curves
        nbTB : number of knots for Top & Bottom curves
        disp = boolean : for display (or not) of the reconstructed curves
        Return control points of each border curve
            and the figure axes 'ax0'
        --> BE CAREFUL : corner points are not necessary identical !!!
    """
    #----------------------------------------------
    #graphic windows for display
    #----------------------------------------------
    if disp :
        fig = plt.figure()
        ax0 = plt.axes(projection='3d')
        ax0.set_xlabel('x axis')
        ax0.set_ylabel('y axis')
        ax0.set_zlabel('z axis')
        ax0.set_title('border curves')
    
    #------------------------------------
    # Cylindrical coordinates
    # on détermine la distance r et l'angle theta 
    # associé à chaque point (x2,y2,z2)
    #------------------------------------

    theta, r = [None for _ in range(len(x2))], [None for _ in range(len(x2))]
    for i in range(len(x2)) :
        theta[i] = np.arctan(x2[i]/z2[i])
        r[i] = np.sqrt(x2[i]*x2[i] + z2[i]*z2[i])

    """------------------------------------------------------------
    DETECTION of points close to the borders for approximation
    ------------------------------------------------------------"""  
    error1 = err1 # angular error
    #--------------------------------------------------------------
    # LEFT CURVE (from observer)
    #--------------------------------------------------------------
    # kept data for approximation
    xLC = []
    yLC = []
    zLC = []
    for i in range(len(x2)) :
        if abs(theta[i] - thetaMin) < error1 :
            xLC.append(x2[i])
            yLC.append(y2[i])
            zLC.append(z2[i])
    
    if disp :
        ax0.scatter(xLC, yLC, zLC, color='r', marker='o', s=0.5) 
    
    #--------------------------------------------------------------
    # RIGHT CURVE (from observer)
    #--------------------------------------------------------------
    # kept data for approximation
    xRC = []
    yRC = []
    zRC = []
    for i in range(len(x2)) :
        if abs(theta[i] - thetaMax) < error1 :
            xRC.append(x2[i])
            yRC.append(y2[i])
            zRC.append(z2[i])
    
    if disp :
        ax0.scatter(xRC, yRC, zRC, color='b', marker='o', s=0.5) 
        
    error2 = err2 # vertical error

    #--------------------------------------------------------------
    # TOP CURVE
    #--------------------------------------------------------------
    # kept data for approximation
    xTC = []
    yTC = []
    zTC = []
    thTC = [] # we keep angle theta associated with each point for reconstruction
    for i in range(len(x2)) :
        if abs(y2[i] - yM) < error2 :
            xTC.append(x2[i])
            yTC.append(y2[i])
            zTC.append(z2[i])
            thTC.append(theta[i])
    
    if disp :
        ax0.scatter(xTC, yTC, zTC, color='g', marker='o', s=0.5) 
        
    #--------------------------------------------------------------
    # BOTTOM CURVE
    #--------------------------------------------------------------
    # kept data for approximation
    xBC = []
    yBC = []
    zBC = []
    thBC = [] # we keep angle theta associated with each point for reconstruction
    for i in range(len(x2)) :
        if abs(y2[i] - ym) < error2 :
            xBC.append(x2[i])
            yBC.append(y2[i])
            zBC.append(z2[i])
            thBC.append(theta[i])
    
    if disp :
        ax0.scatter(xBC, yBC, zBC, color='c', marker='o', s=0.5) 
    
    """------------------------------------------------------------
    APPROXIMATION of Border curves
    ------------------------------------------------------------"""    
    #----------------------
    # Left curve
    #----------------------
    nb_noeudsG = nbLR
    a = ym
    b = yM
    pcxL, pcyL, pczL, xsL, ysL, zsL = curve_approximation(
        xLC, yLC, zLC, yLC, a, b, nb_noeudsG)
    LC = (pcxL, pcyL, pczL)
    if disp :
        # display of control polygon
        ax0.plot3D(pcxL,pcyL,pczL, 'co--', lw=1, ms=5) 
        # B-spline
        ax0.plot3D(xsL,ysL,zsL,'r',lw=1)
    
    #----------------------
    # Right curve
    #----------------------
    nb_noeudsD = nbLR
    a = ym
    b = yM
    pcxR, pcyR, pczR, xsR, ysR, zsR = curve_approximation(
        xRC, yRC, zRC, yRC, a, b, nb_noeudsD)
    RC = (pcxR, pcyR, pczR)
    if disp :
        # display of control polygon
        ax0.plot3D(pcxR,pcyR,pczR, 'co--', lw=1, ms=5) 
        # B-spline
        ax0.plot3D(xsR,ysR,zsR,'r',lw=1)
    
    #----------------------
    # Top curve
    #----------------------
    nb_noeudsH = nbTB
    a = thetaMin
    b = thetaMax
    pcxT, pcyT, pczT, xsT, ysT, zsT = curve_approximation(
        xTC, yTC, zTC, thTC, a, b, nb_noeudsH)
    TC = (pcxT, pcyT, pczT)
    if disp :
        # display of control polygon
        ax0.plot3D(pcxT,pcyT,pczT, 'co--', lw=1, ms=5) 
        # B-spline
        ax0.plot3D(xsT,ysT,zsT,'r',lw=1)
    
    #----------------------
    # Bottom curve
    #----------------------
    nb_noeudsB = nbTB
    a = thetaMin
    b = thetaMax
    pcxB, pcyB, pczB, xsB, ysB, zsB = curve_approximation(
        xBC, yBC, zBC, thBC, a, b, nb_noeudsB)
    BC = (pcxB, pcyB, pczB)
    if disp :
        # display of control polygon
        ax0.plot3D(pcxB,pcyB,pczB, 'co--', lw=1, ms=5) 
        # B-spline
        ax0.plot3D(xsB,ysB,zsB,'r',lw=1) 

    return LC, RC, TC, BC, ax0


#-------------------------------------------------------------------
# Definition of the linear system of constraints : 'couture'
# associated with the border curves of the B-spline surfaces
def ConstraintsSystem(nb_noeudsT,nb_noeudsY,LC,RC,TC,BC):
    """
        Definition of the linear system of constraints 
        associated with the border curves of the B-spline surfaces
        Input :
        nb_noeudsT,nb_noeudsY = number of knots (T for Theta)
        LC = control points of Left curve
        RC = control points of Right curve
        TC = control points of Top curve
        BC = control points of Bottom curve
        Output :
        the matrix H and second members according to x,y,z
    """
    #----------------------------------
    # number of control points
    ni = nb_noeudsT + 2
    nj = nb_noeudsY + 2
    # total number of constraints :
    # (constraints at corners are taken just once)
    nbConstraints = 2*(ni + nj) - 4
    # number of variables (i.e., of control points)
    nbPC = ni*nj
    
    #----------------------------------
    # matrix H of constraints :
    #----------------------------------
    H = np.zeros((nbConstraints,nbPC)) 
    
    # control points according to Bottom, Left, Right and Top curves
    # (prise de tête...)
    """
    i = 0
    for j in range(nj - 1):
        H[i, j] = 1
        i += 1
    for j in range(1, ni):
        H[i, j*nj] = 1
        i += 1
    for j in range(1, ni):
        H[i, j*nj-1] = 1
        i += 1
    for j in range((ni-1)*nj + 1, ni*nj):
        H[i, j] = 1
        i += 1
    """

    for i in range(ni-1):
        H[i, i] = 1
    for i in range(ni-1, nj+ni-2):
        H[i, ni*(i-ni+2)] = 1
    for i in range(ni+nj-2, 2*nj+ni-3):
        H[i, ni*(i-ni-nj+3) - 1] = 1
    for i in range(2*nj+ni-3, 2*(nj+ni)-4):
        H[i, ni*nj + i - 2*(nj+ni) + 4] = 1

    #-------------------------------------------
    # second member :
    #----------------------------------
    # Corner points must be identical
    # so we consider the average value for each corner point
    (pcxL, pcyL, pczL) = LC # Left curve
    (pcxR, pcyR, pczR) = RC # Right curve
    (pcxT, pcyT, pczT) = TC # Top curve
    (pcxB, pcyB, pczB) = BC # Bottom curve
    
    # average at Bottom Left
    xBL = (pcxL[0] + pcxB[0])/2
    yBL = (pcyL[0] + pcyB[0])/2
    zBL = (pczL[0] + pczB[0])/2
    
    pcxL[0] = pcxB[0] = xBL
    pcyL[0] = pcyB[0] = yBL
    pczL[0] = pczB[0] = zBL

    # average at Bottom Right
    xBR = (pcxR[0] + pcxB[-1])/2
    yBR = (pcyR[0] + pcyB[-1])/2
    zBR = (pczR[0] + pczB[-1])/2

    pcxR[0] = pcxB[-1] = xBR
    pcyR[0] = pcyB[-1] = yBR
    pczR[0] = pczB[-1] = zBR
    
    # average at Top Left
    xTL = (pcxL[-1] + pcxT[0])/2
    yTL = (pcyL[-1] + pcyT[0])/2
    zTL = (pczL[-1] + pczT[0])/2

    pcxL[-1] = pcxT[0] = xTL
    pcyL[-1] = pcyT[0] = yTL
    pczL[-1] = pczT[0] = zTL
    
    # average at Top Right
    xTR = (pcxR[-1] + pcxT[-1])/2
    yTR = (pcyR[-1] + pcyT[-1])/2
    zTR = (pczR[-1] + pczT[-1])/2
    
    pcxR[-1] = pcxT[-1] = xTR
    pcyR[-1] = pcyT[-1] = yTR
    pczR[-1] = pczT[-1] = zTR


    """
    for i in range(3) :
        LC[i][-1] = (xBL, yBL, zBL)[i]
        BC[i][0] = (xBL, yBL, zBL)[i]
        RC[i][-1] = (xBR, yBR, zBR)[i]
        BC[i][-1] = (xBR, yBR, zBR)[i]
        LC[i][0] = (xTL, yTL, zTL)[i]
        TC[i][0] = (xTL, yTL, zTL)[i]
        RC[i][0] = (xTR, yTR, zTR)[i]
        TC[i][-1] = (xTL, yTL, zTL)[i]
    """

    # second member
    Bx = np.array([xBL] + [pcxB[i] for i in range(1, len(pcxT) - 1)] + [pcxL[i] for i in range(1, len(pcxL)-1)] + [xTL] + [xBR] + [pcxR[i] for i in range(1, len(pcxR)-1)] + [pcxT[i] for i in range(1, len(pcxB) - 1)] + [xTR])
    By = np.array([yBL] + [pcyB[i] for i in range(1, len(pcyT) - 1)] + [pcyL[i] for i in range(1, len(pcyL)-1)] + [yTL] + [yBR] + [pcyR[i] for i in range(1, len(pcyR)-1)] + [pcyT[i] for i in range(1, len(pcyB) - 1)] + [yTR])
    Bz = np.array([zBL] + [pczB[i] for i in range(1, len(pczT) - 1)] + [pczL[i] for i in range(1, len(pczL)-1)] + [zTL] + [zBR] + [pczR[i] for i in range(1, len(pczR)-1)] + [pczT[i] for i in range(1, len(pczB) - 1)] + [zTR])


    return H, Bx, By, Bz




