import numpy as np
import sys
import os
import math
import Bio.PDB as pdb
import BioHelpers.modules.bio_mod as bm
import BioHelpers.modules.gmap as gm
import BioHelpers.modules.contacts as con

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})


pathToPDB = "../../RNA_Testset/PDB"
all_pdbs = os.listdir(pathToPDB)
lam = 0.05 # Lambda Parameter

def write_res_file(filename: str, contacts: list):
    with open(filename, 'w') as f:
        for con in contacts:
            l = " A/" + str(con[0]+1)+"/N A/" + str(con[1]+1)+"/N "   
            f.write("WELL"+l+"3.5 9.5 1.0 \n")
            f.write("SLOPE"+l+"3.5 9.5 1.0 \n")
            f.write("SLOPE"+l+"0 25 -1.0 \n")



def main():
    # Create new folder
    if not os.path.exists("lambda_"+str(int(lam*100)).zfill(2)):
        os.mkdir("lambda_"+str(int(lam*100)).zfill(2))


    for p in all_pdbs:
        if ".pdb" in p:
            neigh = 2
            struct = pdb.PDBParser().get_structure("", pathToPDB+"/"+p)[0]
            struct = next(struct.get_chains())
            contactMap = bm.calc_contact_matrix(struct, struct, 9.5)
            length = len(bm.sequFromPDB(pathToPDB+"/"+p))
            contactMap = bm.triangularMatrix(contactMap)
            contactMap = bm.deleteNeighbours(contactMap, neigh)
            noc = math.floor(length/2)

            # All contacts
            # ---------------------------------------
            rna_contactMap = gm.arrToList(contactMap)
            # plot = np.array([[p[0],p[1]] for p in rna_contactMap if p[2]==1])
            # fig,ax = plt.subplots(figsize=(7.5,7.5))
            # ax.set_aspect(1)
            # ax.set_xlim([0,max(np.array(rna_contactMap)[:,0])])
            # ax.set_ylim([0,max(np.array(rna_contactMap)[:,1])])
            # ax.set_xlabel(r"$\mathrm{Residue}\ i\ \longrightarrow$")
            # ax.set_ylabel(r"$\mathrm{Residue}\ j\ \longrightarrow$")
            # ax.scatter(plot[:,0],plot[:,1], c='C0', marker='o'); 
            # ax.scatter(plot[:,1],plot[:,0],c='C0', marker='o'); 

            # Random Pick
            # ---------------------------------------
            #contacts = con.pickRandomFalse(contactMap, noc, lam, neigh)
            #write_res_file("res_cluster/"+p.replace(".pdb", ".res"), contacts)
            #ax.scatter([e[0] for e in contacts], [e[1] for e in contacts], c='C1', marker='x'); 


            # Gauss Contacts
            # ---------------------------------------
            contacts = con.pickBestGaussFalse(contactMap, noc, lam, neigh, 10000, 3)
            write_res_file("lambda_"+str(int(lam*100)).zfill(2)+"/"+p.replace(".pdb", ".res"), contacts)
            #ax.scatter([e[1] for e in contacts], [e[0] for e in contacts], color='C1', marker='x'); 

            #plt.show()
            #plt.savefig("Maps_png/"+p.replace(".pdb", ".png"), transparent=True)



if __name__ == "__main__":
    arc = len(sys.argv)
    if (arc > 2) and (sys.argv[1] == "-lam"):
        lam = float(sys.argv[2])
        main()
    else:
        print("No no no", arc, sys.argv)