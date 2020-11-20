import Bio.PDB as PDB
from Bio.PDB import PDBParser
import numpy as np

p = PDBParser()  # create a parser object
s = p.get_structure("2PTC", "2ptc.pdb")  # get the structure from a PDB file

def calc_cen_coord(chain):
    """calculate the coordinate of the enzyme's centre"""
    x = []
    y = []
    z = []
    for residue in chain:
        print(residue.get_id()[1], "begin")
        try:
            res_xyz = tuple(residue['CA'].get_coord())  # get_vector could be better
            x.append(res_xyz[0])  # collect the x-coordinate of all CA atom in a list
            y.append(res_xyz[1])
            z.append(res_xyz[2])
            print(residue.get_id()[1], "over")
        except:  # some "residues" do not have CA atom
            continue

    print("all residues are finished")
    cen_x = (min(x) + max(x)) / 2  # 老师说直接求平均数就行，不用求最大最小的平均数
    cen_y = (min(y) + max(y)) / 2
    cen_z = (min(z) + max(z)) / 2
    return cen_x, cen_y, cen_z


cen_x, cen_y, cen_z = calc_cen_coord(s[0]['E'])  # get the coordinate of trypsin's centre

class AtomSelect(PDB.Select):
    """define a sub-class of the class Select"""
    def accept_atom(self, atom):
        coord = atom.get_coord()
        dist = np.sqrt(np.square(coord[0] - cen_x) + np.square(coord[1] - cen_y) + np.square(coord[2] - cen_z))
        # calculate each atom's distance to the centre
        if (dist <= 4):
            return True
        else:
            return False


io = PDB.PDBIO()  # create a PDBIO object
io.set_structure(s)  # set the structure
io.save("atom_within_sphere_4Å.pdb", select=AtomSelect())
