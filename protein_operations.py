from Bio.PDB import PDBParser
import Bio.PDB as PDB
import os
from Bio.PDB.Polypeptide import standard_aa_names
import random
import numpy as np
from numpy.linalg import *
import matplotlib.pyplot as plt

##  Part A: Define All functions needed
#   1. extract all 9-residue polypeptide segments (2 functions)

def extract_segments(k, dir_pdb):
    """extract all polypeptide segments of a given length from all proteins in a given directory"""
    segment_list = []  # to collect all segments (polypeptides with a length of 9)
    for maindir, subdir, file_name_list in os.walk(dir_pdb):
        for filename in file_name_list:
            try:
                apath = os.path.join(maindir, filename)
                p = PDBParser()
                s = p.get_structure(filename, apath)
                for chain in s[0]:
                    res_list = list(chain.get_residues())  # all residues in the chain
                    for i in range(len(res_list) - k + 1):
                        res = []  # collect 9 continuous residues to form a segment
                        for j in range(i, i + k):
                            res.append(res_list[j])
                            segment_list.append(res)  # collect these polypeptide chains in a list
            except:
                continue  # disregard any structures that cannot be parsed by the Bio.PDB parser
    return segment_list


def clean_segment(segment):
    """remove segments with non-standard amino acid or amino acid lacking any backbone atoms"""
    segments_clean = []
    for seg in segment:
        is_standard = True
        for res in seg:
            if(res.resname not in standard_aa_names) or not (res.has_id('N') and res.has_id('CA') and res.has_id('C')):
                is_standard = False
                break
        if(is_standard == True):
            segments_clean.append(seg)
    return(segments_clean)


#   2. select 500 random segment pairs (4 functions)

def filter_segment_with_cen(aa, segments):
    """collect all segments whose central amino acid is the given one"""
    aa = aa.upper()  # amino acid names in PDB files are upper case.
    seg_list = []  # to collect each segment whose centre is the provided amino acid
    for seg_ in segments:
        if (seg_[4].resname == aa):
            seg_list.append(seg_)
    return seg_list


def select_pairs(seg_list):
    """get 500 random segment pairs from the given segment list
    make sure that segments within each pair are from different proteins"""
    random.shuffle(seg_list)  # make sure it is random
    seg_pairs = []
    num_pair = 0
    while(num_pair < 500):
        for ii in range(1, len(seg_list)):  # get a pair of segments from different proteins
            if(seg_list[0][0].get_full_id()[0] == seg_list[ii][0].get_full_id()[0]):
                continue
            else:
                seg_pairs.append((seg_list[0], seg_list[ii]))  # each element of the list is a tuple
                del(seg_list[ii], seg_list[0]) #  delete used segments, the order matters
                break
        num_pair = num_pair + 1
    return(seg_pairs)


def get_coor_matr(seg):
    """get the coordinate matrix (by column) from a polypeptide segment
    the segment will be centred on the C-alpha atom of its central amino acid"""
    coor_matr = np.matrix([[], [], []])
    for residue in seg:
        # only collect the coordinates of backbone atoms
        coorN = np.matrix(residue['N'].get_coord()).T
        coorCA = np.matrix(residue['CA'].get_coord()).T
        coorC = np.matrix(residue['C'].get_coord()).T
        coor_matr = np.hstack((coor_matr, coorN, coorCA, coorC))  # add a new column of coordinate to the matrix
    CA_cent_coor = np.matrix(seg[4]['CA'].get_coord()).T
    coor_matr = coor_matr - CA_cent_coor  # centre the segment's coordinate matrix
    return coor_matr


def RMSD(X, Y):
    """to find the optimal rotation matrix for superposition and calculate the value of RMSD
    need 2 centred coordinate matrices"""
    R = Y * X.T
    V, S, W = svd(R)
    Z = np.diag([1, 1, -1])
    U = W.T * V.T
    if(det(U) == -1):
        U = W.T * Z * V.T
    E0 = 0
    for i in range(X.shape[1]):
        E0 = E0 + sum(np.array(X[:, i])**2) + sum(np.array(Y[:, i])**2)
    E = E0 - 2 * np.trace(X.T * U * Y)
    RMSD_value = np.sqrt(E / X.shape[1])
    return RMSD_value[0]  # extract the value from a list


#   3. draw a histogram of the RMSD value for each of the 7 amino acid (no function)

#   4. RMSD histogram for 500 randomly chosen superimposed 9-residue segment pairs (no function)

#   5. visualize the 6 segment pairs with lowest RMSD for Proline in PyMol (1 function)

def export_seg_pdb(segmt, k, dir_pdb, out_name):
    """export a given segment to a pdb file, the length of the segment is k,
    the original protein PDB file that the segment comes from is in the given directory"""
    full_id = []
    for a in range(k):  # collect the full ids of all residues in a list
        full_id.append(segmt[a].get_full_id())
    pro_name = segmt[0].get_full_id()[0]  # the protein that the segment comes from
    print("a segment from", pro_name)
    p = PDBParser()
    s = p.get_structure(pro_name, dir_pdb + "/" + pro_name)

    class ModelandResidueSelect(PDB.Select):
        """define a sub-class of the class Select
        discard alternative conformation, thus make sure each amino acid has only one CA atom when parsed by PyMol
        select residues with given id from the original protein that contains the segment"""
        def accept_model(self, model):
            if(model.id == 0):
                return True
            else:
                return False
        def accept_residue(self, residue):
            res_id = residue.get_full_id()
            if(res_id in full_id):
                return True
            else:
                return False

    io = PDB.PDBIO()  # create a PDBIO object
    io.set_structure(s)  # set the structure
    io.save(out_name, select=ModelandResidueSelect())


## Part B: Main Part

if __name__ == "__main__":

    #   1. extract all 9-residue polypeptide segments

    k = 9  # the length of each segment
    dir_pdb = "top100H"  # the directory of the 100 PDB files
    segment_list = extract_segments(k, dir_pdb)
    print("Number of segments extracted: ", len(segment_list))
    segments_clean = clean_segment(segment_list)
    print("Number of segments qualified: ", len(segments_clean))

    #   2. select 500 random segment pairs

    amino_acids = ['Ala', 'Gly', 'Pro', 'Phe', 'Asp', 'Arg', 'Leu']
    RMSD_list = []  # contain 7 elements corresponding to 7 amino acids, each element is a list of 500 RMSD values.
    for aa in amino_acids:
        print(aa, "begin")
        rmsd_list = []  # contain 500 RMSD values for the amino acid.
        sig_list = filter_segment_with_cen(aa, segments_clean)
        pair_list = select_pairs(sig_list)
        for pair in pair_list:
            x = get_coor_matr(pair[0])
            y = get_coor_matr(pair[1])
            rmsd_value = RMSD(x, y)
            rmsd_list.append(rmsd_value)
        print(aa, "finish")
        RMSD_list.append(rmsd_list)
    print(len(RMSD_list), "amino acids finish")

    #   3. draw a histogram of the RMSD value for each of the 7 amino acids

    colour = ['b', 'g', 'r', 'c', 'm', 'y', 'b']
    for i in range(len(RMSD_list)):
        plt.cla()  # remove existing plots to avoid plots overlapping
        plt.hist(np.array(RMSD_list[i]), bins=70, facecolor=colour[i], edgecolor='black', density=1)
        plt.title(amino_acids[i])
        plt.savefig('./' + amino_acids[i] + ".png", format='png')

    #   4. RMSD histogram for 500 randomly chosen superimposed 9-residue segment pairs

    ran_pairs_500 = select_pairs(segments_clean)
    rmsd_list_ran = []
    for pa in ran_pairs_500:
        x = get_coor_matr(pa[0])
        y = get_coor_matr(pa[1])
        rmsd_value_ran = RMSD(x, y)
        rmsd_list_ran.append(rmsd_value_ran)
    plt.cla()
    plt.hist(np.array(rmsd_list_ran), facecolor="blue", edgecolor="black", density=1, bins=70)
    plt.title("regardless of central aa type")
    plt.savefig('./regardless of central aa type.png', format='png')

    #   5. visualize the 6 segment pairs with lowest RMSD values for Proline in PyMol

    rmsd_pair = {}  # for this dictionary,a key is an RMSD value, a value is the corresponding segment pair
    sig_list = filter_segment_with_cen('Pro', segments_clean)
    pair_list = select_pairs(sig_list)
    for pair in pair_list:
        x = get_coor_matr(pair[0])
        y = get_coor_matr(pair[1])
        rmsd_value = RMSD(x, y)
        rmsd_pair[rmsd_value] = pair

    lowest_pairs = []  # collect 6 pairs of segments with lowest RMSD value
    for h in (range(6)):
        key_list_array = np.array(list(rmsd_pair.keys()))  # must convert to list first
        min_rmsd = np.min(key_list_array)
        print("RMSD value", min_rmsd, "selected")
        low_pair = rmsd_pair.pop(min_rmsd)  # record the current minimum RMSD value and remove it from the dictionary
        lowest_pairs.append(low_pair)  # collect the corresponding segment pair
        print("Collected segment pair: ", low_pair)

    for c in range(6):  # export the 6 pairs of (12) segments as PDB files.
        file_name_1 = str(c + 1) + "a" + ".pdb"
        file_name_2 = str(c + 1) + "b" + ".pdb"
        export_seg_pdb(lowest_pairs[c][0], k, dir_pdb, file_name_1)
        export_seg_pdb(lowest_pairs[c][1], k, dir_pdb, file_name_2)
        print("pair", c+1, "exported")


