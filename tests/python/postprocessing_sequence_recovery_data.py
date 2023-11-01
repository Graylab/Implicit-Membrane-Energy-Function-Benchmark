from pyrosetta import *
from pyrosetta.teaching import *
init()


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os.path import exists
from os import path
import os

#Franklin2019
# scratch_dir ='/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2019/sequence-recovery/weights_from_test8_heavyatomcorrection_changedddG_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001'
#Franklin2023
scratch_dir = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/weights_from_test8_heavyatomcorrection_changedddG_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001'
# Read list of energy landscape test cases
list_of_targets = "/home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/"
# list_of_targets = list_of_targets + "targets/design/monomer_chains_design.list"

beta=False
if(beta==True):
    list_of_targets = list_of_targets + "targets/design/beta_monomer_chains.list"
else:
    list_of_targets = list_of_targets + "targets/design/monomer_chains_design.list"
    
with open( list_of_targets, 'rt' ) as f: 
    protein = f.readlines()
    protein = [ x.strip() for x in protein ]

#protein = ['4V1G','4X89','5IVA','4Y25']
#print(protein)
#sort command:
#sort -n -k2 example_score_file | awk '{print $2 "\t" $3}'
from itertools import islice
n_file = 30
count = 0
filename0 = []
filename1 = []
filename3 = []
filename4 = []
filename5 = []
filename6 = []
filename7 = []
#np.array()#np.zeros((n_file,),dtype=str)
print("pdb #residues")
for index in protein:
    print(index)
    os.chdir( scratch_dir +"/"+index )
    analysis_executable = "sort -n -k2 score.sc | head -n 100 | awk " + "'{print $25 " + '"\t"' + " $2}'" + " > score_rmsd"
    os.system(analysis_executable)
    # print("sort -n -k2 score.sc | awk '{print $26"+ " \"\t\" " +"$2}' > score_rmsd")
    # os.system( "sort -n -k2 score.sc | awk '{print $26"+ " " +"$2}' > score_rmsd" )

print("======================")
# for index in protein:
#     file = scratch_dir +"/"+index +"/score_rmsd"
#     # print(file)
#     if(not path.exists(file)):
#         print(index)
# sys.exit()

for index in protein:
    file = scratch_dir +"/"+index +"/score_rmsd"
    if(beta!=True):
        file_output = scratch_dir + "/" + "design.list"
    else:
        file_output = scratch_dir + "/" + "beta_design.list"
        
    if(not path.exists(file)):
        print(index)
        continue
    
    with open(file, 'r') as f:
        with open(file_output,'a') as o:
            lowest_score_structure = list(islice(f,5))
            if(lowest_score_structure==[]):
                continue
            #lowest_protien = lowest_score_structure.split("\t")
            # print(lowest_score_structure)
            if(not lowest_score_structure[1] == " "):
                lowest_protien = lowest_score_structure[0].split("\t")
            else:
                lowest_protien = lowest_score_structure[2].split("\t")
            
            pdb = lowest_protien[0].split("_")
            
            #print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            if(pdb[0]=='description' or pdb[0]==' '):
                continue
            pdb_file = scratch_dir + "/" + pdb[0] + "/" + lowest_protien[0] + ".pdb.gz" 
            print(pdb_file)
            if(not path.exists(pdb_file)):
                continue
            # pose_start = pose_from_file(pdb_file)
            
            #print(pdb[0] + " " + str(pose_start.total_residue()))
            # if(pose_start.total_residue()<100):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename0.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>100 and pose_start.total_residue()<=200):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename1.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>200 and pose_start.total_residue()<=300):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename3.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>300 and pose_start.total_residue()<=400):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename4.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>400 and pose_start.total_residue()<=500):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename5.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>500 and pose_start.total_residue()<=600):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            #     filename6.append(pdb[0]+"/"+lowest_protien[0]+".pdb")
            # elif(pose_start.total_residue()>600):
            #     print(pdb[0]+"/"+lowest_protien[0]+".pdb")
            
            
            filename7.append(pdb[0]+"/"+lowest_protien[0]+".pdb.gz")
            o.write(pdb[0]+"/"+lowest_protien[0]+".pdb.gz"+"\n")
    #count = count+1
            #o.write(lowest_protien[0]+"\n")
# print("proteins<100")
# print(filename0)
# print("100<proteins<200")
# print(filename1)
# print("200<proteins<300")
# print(filename3)
# print("300<proteins<400")
# print(filename4)
# print("400<proteins<500")
# print(filename5)
# print("500<proteins<600")
# print(filename6)
# print("600<proteins")
print(filename7)