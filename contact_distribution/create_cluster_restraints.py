import numpy as np
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

def write_res_file(filename: str, contacts: list):
    with open(filename, 'w') as f:
        for con in contacts:
            l = " A/" + str(con[0]+1)+"/N A/" + str(con[1]+1)+"/N "   
            f.write("WELL"+l+"3.5 9.5 1.0 \n")
            f.write("SLOPE"+l+"3.5 9.5 1.0 \n")
            f.write("SLOPE"+l+"0 25 -1.0 \n")




for p in all_pdbs:
    if ".pdb" in p:
        struct = pdb.PDBParser().get_structure("", pathToPDB+"/"+p)[0]
        struct = next(struct.get_chains())
        contactMap = bm.calc_contact_matrix(struct, struct, 9.5)
        length = len(bm.sequFromPDB(pathToPDB+"/"+p))
        contactMap = bm.triangularMatrix(contactMap)
        contactMap = bm.deleteNeighbours(contactMap, 2)

        # All contacts
        # ---------------------------------------
        rna_contactMap = gm.arrToList(contactMap)
        plot = np.array([[p[0],p[1]] for p in rna_contactMap if p[2]==1])
        fig,ax = plt.subplots(figsize=(7.5,7.5))
        ax.set_aspect(1)
        ax.set_xlim([0,max(np.array(rna_contactMap)[:,0])])
        ax.set_ylim([0,max(np.array(rna_contactMap)[:,1])])
        ax.set_xlabel(r"$\mathrm{Residue}\ i\ \longrightarrow$")
        ax.set_ylabel(r"$\mathrm{Residue}\ j\ \longrightarrow$")
        ax.scatter(plot[:,0],plot[:,1], c='C0', marker='o'); 
        ax.scatter(plot[:,1],plot[:,0],c='C0', marker='o'); 

        # Largest Cluster
        # ---------------------------------------
        noc = math.floor(length/2)
        contacts = con.pickInLargestCluster(contactMap, noc)
        #write_res_file("res_cluster/"+p.replace(".pdb", ".res"), contacts)

        #ax.scatter([e[0] for e in contacts], [e[1] for e in contacts], c='C1', marker='x'); 

        # Largest Cluster (with fillInLargestCluster algorithm)
        # ---------------------------------------
        #contacts = con.fillInLargestCluster(contactMap, noc)
        #write_res_file("res_clusterDense/"+p.replace(".pdb", ".res"), contacts)

        ax.scatter([e[0] for e in contacts], [e[1] for e in contacts], c='C1', marker='x'); 

        # Random
        # ---------------------------------------
        contacts = con.pickRandom(contactMap, noc)
        write_res_file("res_random/"+p.replace(".pdb", ".res"), contacts)

        # Gauss Contacts
        # ---------------------------------------
        #contacts = con.pickBestGauss(contactMap, noc, 10000, 3)
        #write_res_file("res_gauss/"+p.replace(".pdb", ".res"), contacts)
        ax.scatter([e[1] for e in contacts], [e[0] for e in contacts], color='C1', marker='x'); 

        #plt.show()
        #plt.savefig("Maps_png/"+p.replace(".pdb", ".png"), transparent=True)
