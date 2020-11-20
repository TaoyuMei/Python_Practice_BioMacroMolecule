### apply RMSD to 2 point sets a and b ###
from numpy import *
from numpy.linalg import *

a = matrix([[18.92238689, 9.18841188, 8.70764463, 9.38130981, 8.53057997],
            [1.12391951, 0.8707568, 1.01214183, 0.59383894, 0.65155349],
            [0.46106398, 0.62858099, -0.02625641, 0.35264203, 0.53670857]], 'f')

b = matrix([[1.68739355, 1.38774297, 2.1959675, 1.51248281, 1.70793414],
            [8.99726755, 8.73213223, 8.86804272, 8.31722197, 8.9924607],
            [1.1668153, 1.1135669, 1.02279055, 1.06534992, 0.54881902]], 'f')

def centre_the_matrix(matr):
    """centre the 3 by n matrix matr to its centre of mass"""
    centre_of_mass = matr.sum(1) / matr.shape[1]
    centred_matr = matr - centre_of_mass
    return centred_matr


def RMSD(X, Y):
    """to find the optimal rotation matrix U and the value of RMSD
    using both the formula and the rotated coordinates"""
    X = centre_the_matrix(X)
    Y = centre_the_matrix(Y)
    R = Y * X.T
    V, S, W = svd(R)
    Z = diag([1, 1, -1])
    U = W.T * V.T
    if(det(U) == -1):
        U = W.T * Z * V.T
    E0 = 0
    for i in range(X.shape[1]):
        E0 = E0 + sum(array(X[:, i])**2) + sum(array(Y[:, i])**2)

    # use the formula to calculate the E(U) and the RMSD
    E_fo = E0 - 2 * trace(X.T * U * Y)
    RMSD_fo = sqrt(E_fo / X.shape[1])

    # use the coordinates to calculate the E(U) and the RMSD
    L_co = 0
    for i in range(X.shape[1]):
        L_co = L_co + X[:, i].T * U * Y[:, i]
    E_co = E0 - 2 * L_co
    RMSD_co = sqrt(E_co / X.shape[1])

    return U, RMSD_fo, RMSD_co[0, 0]


U, RMSD_fo, RMSD_co = RMSD(a, b)
print("apply RMSD to 2 point sets a and b", "\n")
print("the optimal rotation matrix U: \n", U, "\n")
print("RMSD calculated by the formula: ", RMSD_fo, "\n")
print("RMSD calculated by the coordinates: ", RMSD_co, "\n")


### apply RMSD to a real PDB file ###
from Bio.PDB import PDBParser

p = PDBParser()  # create a parser object
s = p.get_structure("1LCD", "1lcd.pdb")  # get the structure from a PDB file

def get_coor_matr(chain):
    """get the coordinates matrix (by column) from a peptide chain"""
    coor_matr = matrix([[], [], []])
    for residue in chain:
        if residue.has_id('CA'):
            coor = matrix(residue['CA'].get_coord())  # only collect the coordinate of CA atom
            coor = coor.reshape((3, 1))
            coor_matr = hstack((coor_matr, coor))  # add a new column of coordinate to the matrix
    return coor_matr


# get the coordinates of chain A in model 0 and model 1
m0_matr = get_coor_matr(s[0]['A'])
m1_matr = get_coor_matr(s[1]['A'])

U, RMSD_fo, RMSD_co = RMSD(m0_matr, m1_matr)
print("apply RMSD to a real PDB file", "\n")
print("the optimal rotation matrix U: \n", U, "\n")
print("RMSD calculated by the formula: ", RMSD_fo, "\n")
print("RMSD calculated by the coordinates: ", RMSD_co, "\n")