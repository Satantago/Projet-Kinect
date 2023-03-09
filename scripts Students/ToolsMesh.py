# ToolMesh.py

"""========================================================================
Bibliothèque "mesh""
de fonctions pour la reconstruction géométrique (projet Kinect)

--> Récupération des sommets d'un mesh à partir d'un fichier texte .ply
    def CoordonneesFromPly(FilePly):

--> Delete duplicate integers in a sorted list
    def DeleteDuplicates(r1):        

--> Creation of a PLY file associated with the Bspline surface
    def BsplineToPly(FilePly, Sx, Sy, Sz):

--> Creation of the whole mesh .ply after Bspline reconstruction :
    The mesh is created in the NEW FRAME (after translation and rotations)
    def WholeMeshToPly(FilePly,x2,y2,z2,t3V,Sx, Sy, Sz):

========================================================================"""
import numpy as np


#--------------------------------------------------------------------
# Récupération des sommets d'un mesh à partir d'un fichier texte .ply
def CoordonneesFromPly(FilePly):
    """ Récupération des coordonnées (x,y,z) des sommets 
        d'un mesh triangulaire à partir d'un fichier texte .ply
        Input : FilePly = string = fichier .ply (sans le suffixe .ply)
        Output : x,y,z = 3 listes de float
    """
    FilePly1 = FilePly + '.ply'
    file = open(FilePly1, "r")
    
    # lecture des 10 lignes d'entête :
    # (car fichier non texturé, initialement fichier STL)
    entete = []
    for k in range(10):
        ligne = file.readline()
        entete.append(ligne)
    
    # On récupère le nombre de vertex ligne 4 (3 en python)
    ligne3 = entete[3].split()
    NbVertex = int(ligne3[2])
    
    # On récupère les coordonnées x,y,z des sommets du maillage
    x = []
    y = []
    z = []
    for k in range(10,10+NbVertex):
        ligne = file.readline()
        coords = ligne.split()
        x.append(float(coords[0]))
        y.append(float(coords[1]))
        z.append(float(coords[2]))

    file.close()
    return x,y,z


#-------------------------------------------------------------
# Delete duplicate integers in a sorted list
def DeleteDuplicates(r1):
    """
        Delete duplicate integers in a sorted list
        Input :  
            r1 (= list) is a increasing sorted list of integers
        Output : 
            r2 (= list) is the list of strictly increasing values of r1
    """
    N = len(r1)
    r2 = [r1[0]]
    for i in range(1,N):
        if r1[i] != r1[i-1] :
            r2.append(r1[i])    # we keep this element
    return r2



#-------------------------------------------------------------
# Creation of a PLY file associated with the Bspline surface
# The mesh is created in the NEW FRAME (after translation and rotations)
def BsplineToPly(FilePly, Sx, Sy, Sz):
    """ 
        Creation of a PLY file associated with the Bspline surface
        Input : FilePly (string) = initial file.ply
                Sx, Sy, Sz : matrices of points coordinates (Bspline surface)
        Output : FilePly3 (string) = file.ply
    """
    nbTEval, nbYEval = np.shape(Sx)
    
    Vertices3 = []
    for i in range(nbTEval):
        for j in range(nbYEval):
            line = str(Sx[i,j]) + ' ' + str(Sy[i,j]) + ' ' + str(Sz[i,j]) + '\n'
            Vertices3.append(line)
    
    Faces3 = []
    for i in range(1,nbTEval):
        for j in range(nbYEval-1):
            v1 = i*nbYEval + j
            v2 = (i-1)*nbYEval + j + 1
            v3 = (i-1)*nbYEval + j
            line = '3 ' + str(v1) + ' ' + str(v2) + ' ' + str(v3) + '\n'
            Faces3.append(line)
            v2 = i*nbYEval + j + 1
            v3 = (i-1)*nbYEval + j + 1
            line = '3 ' + str(v1) + ' ' + str(v2) + ' ' + str(v3) + '\n'
            Faces3.append(line)

    nbVertices3 = nbTEval * nbYEval
    nbFaces3 = 2 * (nbTEval -1) * (nbYEval - 1)
    entete = [
    'ply\n'
    'format ascii 1.0\n'
    'comment VCGLIB generated\n'
    'element vertex ' + str(nbVertices3) + '\n'
    'property float x\n'
    'property float y\n'
    'property float z\n'
    'element face ' + str(nbFaces3) + '\n'
    'property list uchar int vertex_indices\n'
    'end_header\n'
    ]

    # Bspline mesh as a string :
    BsplineMesh = entete[0]
    # we concatenate string's elements of Vertices3 and Faces3
    for stg in Vertices3:
        BsplineMesh += stg
    for stg in Faces3:
        BsplineMesh += stg
    
    # we save the new mesh in a file
    FilePly3 = FilePly + '3.ply'
    newfile = open(FilePly3, "w")
    newfile.write(BsplineMesh)
    newfile.close()
    # nothing essential to return


##########################################################################
# Creation of the whole mesh .ply after Bspline reconstruction :
# The mesh is created in the NEW FRAME (after translation and rotations)
# 1) we remove vertices with index t3V in the initial file FilePly
#    and we update faces
# 2) we create the ply mesh associated with the B-spline surface
##########################################################################
def WholeMeshToPly(FilePly,x2,y2,z2,t3V,Sx, Sy, Sz):
    """ 
        Input : FilePly (string) = fichier.ply
                x2,y2,z2 : coordinates after translation & rotations
                t3V list of vertices to remove
                Sx, Sy, Sz : matrices of points coordinates (Bspline surf)
        Output : FilePly4 (string) = fichier.ply
    """
    FilePly1 = FilePly + '.ply'
    file = open(FilePly1, "r")
    
    """-------------------------------------------
    PART 1 : we remove vertices with index t3V in the initial file FilePly
             and we update faces
             ==> complementary part (associated with the green points)
    ----------------------------------------------"""
    # reading of the 10 header lines:
    # (non textured file, initially STL file)
    entete = []
    for k in range(10):
        line = file.readline()  # line is a string
        entete.append(line)     # entete is a list of strings
    
    # we read the number of vertex in line 4 (index 3 in python)
    line3 = entete[3].split()   # line3 is now a list of strings
    NbVertex = int(line3[2])
    # and the number of faces in line 8 (index 7 in python)
    line7 = entete[7].split()
    NbFace = int(line7[2])
    
    # We load all the following lines (Vertices + Faces) 
    # in the list meshVF (= list of strings)
    NbVF = NbVertex + NbFace
    meshVF = []
    for k in range(NbVF):
        line = file.readline()  # line is a string
        meshVF.append(line)     # meshVF is a list of strings
    file.close()

    # we replace coordinates of vertices by x2,y2,z2
    for k in range(NbVertex):
        meshVF[k] = str(x2[k]) + ' ' + str(y2[k]) + ' ' + str(z2[k]) + '\n'

    # Before removing vertices of t3V, we define new index of each vertex :
    # ind[k] is the new index of vertex k (for k not in t3V)
    t3V = sorted(t3V)
    ind = np.arange(NbVertex)
    for k in t3V:
        ind[k+1:] = ind[k+1:] - 1

    # we look for t3F = list of faces to remove
    # (i.e. faces which contain a vertex from t3V)
    t3F = [] 
    for i in range(NbVertex, NbVF):
        line = meshVF[i].split()
        v1 = int(line[1])
        v2 = int(line[2])
        v3 = int(line[3])
        for k in t3V:
            # if one vertex is in t3F we will delete the current face
            if v1 == k or v2 == k or v3 == k:
                t3F.append(i)
        # we update the vertex indices
        v1 = ind[v1]
        v2 = ind[v2]
        v3 = ind[v3]
        meshVF[i] = '3 ' + str(v1) + ' ' + str(v2) + ' ' + str(v3) + '\n'

    # we remove duplicates in t3F
    t3F = sorted(t3F)
    t3F = DeleteDuplicates(t3F)
    t3F = sorted(t3F, reverse=True)
    
    # there are no duplicates in t3V
    t3V = sorted(t3V, reverse=True)

    # we remove vertices (of t3V) and faces (of t3F)
    for k in t3F:
        del meshVF[k]
    for k in t3V:
        del meshVF[k]

    # number of vertices and faces remaining in this complementary part
    nbVertices3c = NbVertex - len(t3V)
    nbFaces3c = NbFace - len(t3F)
    

    """-------------------------------------------
    PART 2 : we create the ply mesh associated with the B-spline surface
    ----------------------------------------------"""
    nbTEval, nbYEval = np.shape(Sx)
    
    # vertices of the B-spline mesh :
    Vertices3 = []
    for i in range(nbTEval):
        for j in range(nbYEval):
            line = str(Sx[i,j]) + ' ' + str(Sy[i,j]) + ' ' + str(Sz[i,j]) + '\n'
            Vertices3.append(line)
    nbVertices3 = nbTEval * nbYEval

    # Faces of the B-spline mesh :
    # -- each rectangle is subdivided in 2 triangles
    # -- and indices of vertices are shifted from the value "nbVertices3c"
    #    which is the number of vertices of the complementary part
    Faces3 = []
    for i in range(1,nbTEval):
        for j in range(nbYEval-1):
            v1 = nbVertices3c + i*nbYEval + j
            v2 = nbVertices3c + (i-1)*nbYEval + j + 1
            v3 = nbVertices3c + (i-1)*nbYEval + j
            line = '3 ' + str(v1) + ' ' + str(v2) + ' ' + str(v3) + '\n'
            Faces3.append(line)
            v2 = nbVertices3c + i*nbYEval + j + 1
            v3 = nbVertices3c + (i-1)*nbYEval + j + 1
            line = '3 ' + str(v1) + ' ' + str(v2) + ' ' + str(v3) + '\n'
            Faces3.append(line)
    nbFaces3 = 2 * (nbTEval -1) * (nbYEval - 1)


    """-------------------------------------------
    PART 3 : we create the full mesh .ply
    ----------------------------------------------"""
    # total number of vertices and faces
    nbVertTotal = nbVertices3c + nbVertices3
    nbFaceTotal = nbFaces3c + nbFaces3

    # HEADER
    # we update the number of vertices and faces in the header
    line3 = entete[3].split()
    line3[2] = str(nbVertTotal)
    entete[3] = line3[0] + ' ' + line3[1] + ' ' + line3[2] + '\n'
    
    line7 = entete[7].split()
    line7[2] = str(nbFaceTotal)
    entete[7] = line7[0] + ' ' + line7[1] + ' ' + line7[2] + '\n'   

    # VERTICES
    # new mesh (is a list of strings):
    newMesh = entete + meshVF[0:nbVertices3c] + Vertices3
    
    # FACES
    newMesh += meshVF[nbVertices3c:] + Faces3

    # we concatenate string's elements of the list
    newMeshString = '' 
    for stg in newMesh:
        newMeshString += stg

    # we save the new mesh in a file
    FilePly4 = FilePly + '4.ply'
    newfile = open(FilePly4, "w")
    newfile.write(newMeshString)
    newfile.close()
    # nothing essential to return....






