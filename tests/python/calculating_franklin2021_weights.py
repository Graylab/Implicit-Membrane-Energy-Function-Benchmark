#!/usr/bin/env python
import sys
import os
import read_config
import glob
import math
import fileinput
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))


#note: for fa_w_b 0.966, before correction of nheavy atoms, the folder was fa_imm_modified_withgp
def main():
    config = read_config.read_config()
    #-0.114,1.941, 1.308
    #fa_water_to_bilayer_ompla = np.array([-0.114,-2.77,0.877,4.025,-2.22,4.727,-2.221,-1.357,2.872,2.25,3.058,1.548,0.802,-1.483,-1.001,-1.689])
    outdir = config.benchmark_path + "data/franklin2021/" + "ddG-of-mutation/afternheavyatomcorrection/fa_wb_1.0_felecbilayer_1.0_fimm_1.0/"
    ompla_file = outdir + "C1_OmpLA_canonical_ddGs/breakdown.sc" 
    df_ompla = pd.read_csv(ompla_file, delimiter=" ")
    df_ompla = df_ompla[(df_ompla['mutation']!='A181A') & (df_ompla['mutation']!='A181D') & (df_ompla['mutation']!='A181E') & (df_ompla['mutation']!='A181P')]
    df_ompla['total_except3'] = df_ompla.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    print(df_ompla)
    
    fa_water_to_bilayer_ompla_v2 = np.array(df_ompla['fa_water_to_bilayer:'])
    # print(fa_water_to_bilayer_ompla_v2)
    # fa_water_to_bilayer_ompla_v2 = np.array([-0.117,-2.783,0.874,4.013,-2.233,4.715,-2.233,-1.37,2.869,2.239,3.046,1.546,0.79,-1.496,-1.027,-1.702])
    # # print(fa_water_to_bilayer_ompla_v2)
    # sys.exit()
    
    # total_ompla_v2 = np.array([2.365,3.804,5.38,2.847,3.448,1.093,1.844,2.393,3.437,1.871,2.192,0.73,1.624,2.468,-9.3,3.73])
    total_ompla_v2 = np.array(df_ompla['total_except3'])
    
    #fa_imm_elec_ompla = np.array([-1.932,-3.766,-1.431,-5.186,-2.833,-1.296,-4.85,-2.146,-1.734,-1.822,-1.976,-1.193,-7.235,-2.274,-5.143,-4.214])
    fa_imm_elec_ompla_v2 = np.array(df_ompla['fa_imm_elec:'])
    # fa_imm_elec_ompla_v2 = np.array([-0.462,-0.387,-0.127,-0.813,-0.718,0.016,-0.669,-0.693,0.678,-0.703,0.046,-0.965,-0.996,-0.296,-0.027,-0.343])
    
    
    #fa_elec_ompla = np.array([-0.147,0.015,-0.34,-0.157,0.085,-0.161,-0.165,0.143,-0.033,0.125,-0.457,-0.648,-0.036,0.131,-0.053])
    
    
    #f_elec_lipidlayer_ompla = np.array([-0.644,-3.987,1.586,-2.017,-5.008,-2.013,-4.998,-3.28,-2.23,-5.621,-3.357,-1.636,-2.966,-3.271,-5.277,-2.996])
    f_elec_lipidlayer_ompla_v2 = np.array(df_ompla['f_elec_lipidlayer:'])
    # f_elec_lipidlayer_ompla_v2 = np.array([-0.607,-3.989,1.629,-2.108,-5.01,-2.015,-5,-3.282,-2.187,-5.702,-3.36,-1.593,-2.968,-3.273,-5.21,-2.998])
    
    experimental_ddG_ompla = np.array([0.49,-2.2,1.72,4.76,-1.56,5.39,-1.81,-0.76,3.47,3.01,3.71,1.83,1.78,-0.78,-0.38,-1.09])
            
    outdir = config.benchmark_path + "data/franklin2021/" + "ddG-of-mutation/afternheavyatomcorrection/fa_wb_1.0_felecbilayer_1.0_fimm_1.0/"
    pagp_file = outdir + "C2_PagP_canonical_ddGs/breakdown.sc" 
    df_pagp = pd.read_csv(pagp_file, delimiter=" ")
    df_pagp = df_pagp[(df_pagp['mutation']!='A104A') & (df_pagp['mutation']!='A104D') & (df_pagp['mutation']!='A104E') & (df_pagp['mutation']!='A104P')]
    df_pagp['total_except3'] = df_pagp.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    print(df_pagp)
    
    # total_pagp_v2 =np.array([0.555,-0.994,2.796,-0.763,-1.328,1.472,2.234,0.96,-0.179,0.171,-0.14,-0.776,-2.193,-2.053,-1.141,-1.884]) 
    total_pagp_v2 = np.array(df_pagp['total_except3'])
    
    f_elec_lipidlayer_pagp_v2 = np.array(df_pagp['f_elec_lipidlayer:'])
    # f_elec_lipidlayer_pagp_v2 = np.array([-0.346,-3.466,1.136,-1.944,-3.704,-1.7,-3.831,-2.337,-1.781,-4.599,-3.166,-0.993,-2.256,-2.426,-4.42,-3.092])
    print(f_elec_lipidlayer_pagp_v2)
    
    fa_water_to_bilayer_pagp_v2 = np.array(df_pagp['fa_water_to_bilayer:'])
    print(fa_water_to_bilayer_pagp_v2)
    # fa_water_to_bilayer_pagp_v2 = np.array([0.019,-2.025,0.116,2.716,-0.975,4.642,-1.431,-0.276,1.4,1.473,3.038,0.027,0.509,-0.297,0.566,-0.962])
    
    fa_imm_elec_pagp_v2 = np.array(df_pagp['fa_imm_elec:'])
    print(fa_imm_elec_pagp_v2)
    # fa_imm_elec_pagp_v2 = np.array([-0.195,0.33,0.095,0.012,-0.032,0.139,0.413,-0.319,0.023,-0.68,-0.014,-0.429,-0.38,0.115,-0.19,0.095])
    # sys.exit()
    #-0.72,2.49, 1.18
    experimental_ddG_pagp = np.array([-0.72,-2.44,1.64,3.32,-2.17,3.54,-2.01,-1.15,2.95,2.54,3.22,1.83,0.95,-1.75,-2.21,-1.02])
    #coefficients = np.array([0.0,0.0,0.0])
    for i in range(len(experimental_ddG_pagp)-2):
        if(experimental_ddG_pagp[i]==0):
            print(i)
            continue
        A = np.array([[fa_water_to_bilayer_pagp[i], fa_imm_elec_pagp[i], f_elec_lipidlayer_pagp[i]],
                      [fa_water_to_bilayer_pagp[i+1], fa_imm_elec_pagp[i+1], f_elec_lipidlayer_pagp[i+1]],
                      [fa_water_to_bilayer_pagp[i+2], fa_imm_elec_pagp[i+2], f_elec_lipidlayer_pagp[i+2]]])
       #               [fa_water_to_bilayer[i+3], fa_imm_elec[i+3], f_elec_lipidlayer[i+3]]])
        B = np.array([experimental_ddG_pagp[i]-total_pagp[i], experimental_ddG_pagp[i+1]-total_pagp[i+1], experimental_ddG_pagp[i+2]-total_pagp[i+2]])

        #print(A)
        #print(B)
        C = np.linalg.solve(A,B)
        if(i==0):
            coefficients = C 
        else:
            coefficients = np.vstack((coefficients, C))
        #print(C)
    coefficient_final = pd.DataFrame(coefficients, columns = ['fa_water_to_bilayer','fa_imm_elec','fa_elec_lipidlayer'])
    print(coefficient_final)
    coefficient_final.to_csv('Pagp_coefficient.csv',index=False)
    
    #hist = coefficient_final['fa_water_to_bilayer', 'fa_imm_elec', 'fa_elec_lipidlayer'].hist()
    #hist_faie = coefficient_final['fa_imm_elec'].hist()
    #hist_fael = coefficient_final['fa_elec_lipidlayer'].hist()
    #print(hist)
    #plt.show()
    #plt.savefig('hist_ompla_coefficients.png', bbox_inches = 'tight')

if __name__ == "__main__": main()
