import numpy as np
import random
import Bio.PDB as pdb
import BioHelpers_FABER.bio_mod as bm

import warnings

warnings.filterwarnings("ignore")


def getContactFromResFile(line: str) -> tuple[int]:
    parts = line.split(sep=" ")
    i = int(parts[1].split(sep="/")[1])
    j = int(parts[2].split(sep="/")[1])
    return (i, j)


def getRestraints(res_file: str) -> list[tuple[int]]:
    out = []
    with open(res_file, "r") as f:
        for line in f:
            out.append(getContactFromResFile(line))
            next(f)
            next(f)
    return out


def writeCylinder(restraint, atomnames, color):
    out = ""
    out += (
        'lassign [[atomselect top "resid '
        + str(restraint[0])
        + " and name "
        + atomnames[0]
        + '"] get {x y z}] coord1 \n'
    )
    out += (
        'lassign [[atomselect top "resid '
        + str(restraint[1])
        + " and name "
        + atomnames[1]
        + '"] get {x y z}] coord2 \n'
    )
    out += "draw color " + color + "\n"
    out += "draw cylinder $coord1 $coord2 radius 0.6 filled yes\n"
    return out


def getResInfo(restraint, filename):
    rna_structure = pdb.PDBParser().get_structure(" ", filename)
    rna_chain = next(rna_structure.get_chains())
    return (
        rna_chain[int(restraint[0])].get_resname(),
        rna_chain[int(restraint[1])].get_resname(),
    ), bm.calc_residue_dist(rna_chain[int(restraint[0])], rna_chain[int(restraint[1])])


def resnameToAtomname(resname):
    if resname in ["A", "G"]:
        return "N9"
    else:
        return "N1"


def disToColor(d):
    if d > 9.5:
        return "red2"
    else:
        return "green2"


def preamble(file: str = "preamble.tcl") -> str:
    with open(file, "r") as f:
        return "".join(f.readlines())


# def spins(filename: str, radius: float = 0.6) -> str:
#     out = []
#     rna_structure = pdb.PDBParser().get_structure(" ", filename)
#     rna_chain = next(rna_structure.get_chains())
#     residues = rna_chain.get_residues()
#     for i, residue in enumerate(residues):
#         out.append(f'set coord [atomselect top "resid {i+1}"]\n')
#         out.append(f"set coord_spin [geom_center $coord]\n")
#         out.append(f"set rot_mat [transabout {{1 1 0}} 30]\n")
#         out.append(f"draw color {resname_to_color(residue.get_resname())}\n")
#         out.append(
#             f"draw cylinder [vecadd $coord_spin [vectrans $rot_mat {{-4 0 0}}]] [vecadd $coord_spin [vectrans $rot_mat {{4 0 0}}]] radius {radius} filled yes\n"
#         )
#         out.append(
#             random.choice(
#                 [
#                     f"draw cone [vecadd $coord_spin [vectrans $rot_mat {{4 0 0}}]] [vecadd $coord_spin [vectrans $rot_mat {{6 0 0}}]] radius {radius*2}\n",
#                     f"draw cone [vecadd $coord_spin [vectrans $rot_mat {{-4 0 0}}]] [vecadd $coord_spin [vectrans $rot_mat {{-6 0 0}}]] radius {radius*2}\n",
#                 ]
#             )
#         )
#     return "".join(out)


def resname_to_color(resname: str) -> str:
    if resname == "A":
        return "yellow"
    elif resname == "C":
        return "green"
    elif resname == "G":
        return "blue"
    else:
        return "red"


def main():
    target = "4frg"
    res_file1 = "res_clusterDense/4frg.res"
    res_file2 = "res_gauss/4frg.res"
    solFileName = "../../RNA_Testset/PDB/" + target + ".pdb"
    noc = bm.numberOfResidues(solFileName)
    restraints1 = getRestraints(res_file1)
    restraints2 = getRestraints(res_file2)
    # with open("vmd_visualisation_" + target + ".tcl", "w") as f:
    with open("test.tcl", "w") as f:
        f.write(preamble())
        for r in restraints1:
            resNames, distance = getResInfo(r, solFileName)
            f.write(
                writeCylinder(
                    r,
                    (resnameToAtomname(resNames[0]), resnameToAtomname(resNames[1])),
                    "green",
                )
            )
        for r in restraints2:
            resNames, distance = getResInfo(r, solFileName)
            f.write(
                writeCylinder(
                    r,
                    (resnameToAtomname(resNames[0]), resnameToAtomname(resNames[1])),
                    "red",
                )
            )


if __name__ == "__main__":
    main()
