#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import read_config
import math
# from sklearn.linear_model import LinearRegression
from random import seed
from random import randint
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import Bounds

import matplotlib.pyplot as plt
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))
import itertools as it

def skiprows_df(file_):
    with open(file_, 'r') as f:        
        lines = f.readlines()
        seqs = [t.rstrip() for t in lines if t.find('fa')==-1]
        seqsnp = np.array(seqs)
    return seqsnp

def removelines_todf(file_,startpattern):
    
    df = pd.read_csv(file_, delimiter=" ", header=None)
    index_list = df.index[df[0]==startpattern].tolist()
   
    df_pagp = pd.read_csv(file_, delimiter=" ", skiprows=index_list)
    df_pagp.columns=['fa_water_bilayer', 'fa_imm_elec','fa_elec_bilayer','abs_error_sum', 'corr_coeff', 'corr_ompla', 'corr_pagp']
    df_pagp_final = df_pagp.sort_values(by=['abs_error_sum'], ascending=True)
    df_pagp_final.to_csv(file_, sep=" ", index=False)
    

def main():
    # config = read_config.read_config()
    # wts=['0.938','0.32','0.185']
    # wts=['1.044','0.107','0.125']
    # wts=['0.91','0.121','0.046']
    # 1.044	0.107	0.125
    # ç	0.121	0.046
    
    # outdir = config.benchmark_path + "data/franklin2021/" + "ddG-of-mutation/afternheavyatomcorrection/fa_wb_"+wts[0]+"_felecbilayer_"+wts[2]+"_fimm_"+wts[1]+"/"
    # print(outdir)
    # sys.exit
    # outfile = outdir + 'coefficient_f2021_wts_basedonompla.dat'
    # removelines_todf(outfile,'fa_water_bilayer')
    
    # outfile = outdir + 'coefficient_f2021_wts_basedoncombined_ompla+heightcombined.dat'
    # removelines_todf(outfile,'fa_water_bilayer')
    
    # outfile = outdir + 'coefficient_f2021_wts_basedoncombined_pagpompla+heightcombined.dat'
    # removelines_todf(outfile,'fa_water_bilayer')
 
    # outfile = outdir + "C1_OmpLA_canonical_ddGs/ddG_franklin2021.dat"
    # df_ompla = pd.read_csv(outfile, delimiter=" ")
    # df_ompla = df_ompla[(df_ompla['Mut']!='P') & (df_ompla['Mut']!='A')]
    # print(df_ompla)
    # sys.exit()
    # df_ompla['diff'] = df_ompla['predicted_ddG']-df_ompla['experimental_ddG']
    # print('ompla sum of error: {}'.format(df_ompla['diff'].pow(2).sum()))
    # df_ompla.plot(kind='scatter',x='experimental_ddG',y='predicted_ddG', c='Black')
    # plt.tight_layout()
    # outpng = '{}/ompla_predicted_ddG.png'.format(outdir)
    # plt.savefig(outpng, transparent=True, dpi=600)
    # plt.close()
    
    # outfile = outdir + "C2_PagP_canonical_ddGs/ddG_franklin2021.dat"
    # df_pagp = pd.read_csv(outfile, delimiter=" ")
    # df_pagp = df_pagp[(df_pagp['Mut']!='P') & (df_pagp['Mut']!='A')]
    # df_pagp['diff'] = df_pagp['predicted_ddG']-df_pagp['experimental_ddG']
    # print('Pagp sum of error: {}'.format(df_pagp['diff'].pow(2).sum()))
    # df_pagp.plot(kind='scatter',x='experimental_ddG',y='predicted_ddG',c='Black')
    # plt.tight_layout()
    # outpng = '{}/pagp_predicted_ddG.png'.format(outdir)
    # plt.savefig(outpng, transparent=True, dpi=600)
    # plt.close()
    calculate_weights()
    # calculate_dg_rosetta_atom_type()
    
    
def calculate_weights():
    config = read_config.read_config()
    
    outdir = config.benchmark_path + "data/franklin2021/" + "ddG-of-mutation/afternheavyatomcorrection/fa_wb_1.0_felecbilayer_1.0_fimm_1.0/"
    ompla_file = outdir + "C1_OmpLA_canonical_ddGs/breakdown.sc" 
    df_ompla = pd.read_csv(ompla_file, delimiter=" ")
    df_ompla = df_ompla[(df_ompla['mutation']!='A181A') & (df_ompla['mutation']!='A181D') & (df_ompla['mutation']!='A181E') & (df_ompla['mutation']!='A181P')]
    df_ompla['total_except3'] = df_ompla.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    print(df_ompla)
    
    fa_water_to_bilayer_ompla_v2 = np.array(df_ompla['fa_water_to_bilayer:'])
    temp=np.array([fa_water_to_bilayer_ompla_v2.T])
    fa_water_to_bilayer_ompla_v2T = np.array(temp.T)
    # fa_water_to_bilayer_ompla_v2 = np.array([-0.114,-2.77,0.877,4.025,-2.22,4.727,-2.221,-1.357,2.872,2.251,3.058,1.548,0.802,-1.483,-1.001,-1.689])
    
    total_ompla_v2 = np.array([0.78,2.218,3.792,1.201,1.787,-0.573,0.236,0.732,1.856,0.21,0.532,-0.83,-0.043,0.804,2.402,2.139])
    total_ompla_v2 = np.array(df_ompla['total_except3'])
    temp=np.array([total_ompla_v2.T])
    total_ompla_v2T = np.array(temp.T)
    
    # fa_imm_elec_ompla_v2 = np.array([-0.447,-0.321,-0.082,-0.782,-0.678,0.042,-0.616,-0.658,0.717,-0.671,0.085,-0.938,-0.955,-0.259,-0.237,-0.294])
    fa_imm_elec_ompla_v2 = np.array(df_ompla['fa_imm_elec:'])
    temp=np.array([fa_imm_elec_ompla_v2.T])
    fa_imm_elec_ompla_v2T = np.array(temp.T)

    # f_elec_lipidlayer_ompla_v2 = np.array([-0.644,-3.987,1.585,-2.106,-5.008,-2.013,-4.998,-3.28,-2.23,-5.7,-3.358,-1.636,-2.966,-3.271,-5.277,-2.996])
    f_elec_lipidlayer_ompla_v2 = np.array(df_ompla['f_elec_lipidlayer:'])
    temp=np.array([f_elec_lipidlayer_ompla_v2.T])
    f_elec_lipidlayer_ompla_v2T = np.array(temp.T)
    
    experimental_ddG_ompla = np.array([0.49,-2.2,1.72,4.76,-1.56,5.39,-1.81,-0.76,3.47,3.01,3.71,1.83,1.78,-0.78,-0.38,-1.09])
    temp=np.array([experimental_ddG_ompla.T])
    experimental_ddG_omplaT = np.array(temp.T)

    #print(len(fa_water_to_bilayer_pagp))
    #print(len(total_pagp))
    #print(len(fa_imm_elec_pagp))
    #print(len(f_elec_lipidlayer_pagp))
    #print(len(experimental_ddG_pagp))
    temp_ompla = np.hstack(( fa_water_to_bilayer_ompla_v2T,fa_imm_elec_ompla_v2T ))       
    A_ompla = np.hstack((temp_ompla, f_elec_lipidlayer_ompla_v2T))
    

    pagp_file = outdir + "C2_PagP_canonical_ddGs/breakdown.sc" 
    df_pagp = pd.read_csv(pagp_file, delimiter=" ")
    df_pagp = df_pagp[(df_pagp['mutation']!='A104A') & (df_pagp['mutation']!='A104D') & (df_pagp['mutation']!='A104E') & (df_pagp['mutation']!='A104P')]
    df_pagp['total_except3'] = df_pagp.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    print(df_pagp)
    
    # total_pagp_v2 =np.array([0.585,-0.962,2.834,-0.726,-1.282,1.511,2.273,0.972,-0.146,0.212,-0.098,-0.736,-2.15,-2.009,-1.605,-1.767]) 
    total_pagp_v2 = np.array(df_pagp['total_except3'])
    temp=np.array([total_pagp_v2.T])
    total_pagp_v2T = np.array(temp.T)
    
    f_elec_lipidlayer_pagp_v2 = np.array(df_pagp['f_elec_lipidlayer:'])
    # f_elec_lipidlayer_pagp_v2 = np.array([-0.345,-3.465,1.137,-1.943,-3.703,-1.699,-3.83,-2.336,-1.78,-4.598,-3.165,-0.992,-2.256,-2.425,-4.42,-3.091])
    temp=np.array([f_elec_lipidlayer_pagp_v2.T])
    f_elec_lipidlayer_pagp_v2T = np.array(temp.T)
    
    # fa_water_to_bilayer_pagp_v2 = np.array([0.018,-2.025,0.116,2.716,-0.975,4.642,-1.431,-0.276,1.4,1.473,3.038,0.027,0.508,-0.297,0.566,-0.962])
    fa_water_to_bilayer_pagp_v2 = np.array(df_pagp['fa_water_to_bilayer:'])
    temp=np.array([fa_water_to_bilayer_pagp_v2.T])
    fa_water_to_bilayer_pagp_v2T = np.array(temp.T)
    
    
    # fa_imm_elec_pagp_v2 = np.array([-0.201,0.332,0.098,0.014,-0.026,0.142,0.415,-0.315,0.013,-0.675,-0.009,-0.423,-0.379,0.12,-0.096,0.106])
    fa_imm_elec_pagp_v2 = np.array(df_pagp['fa_imm_elec:'])
    temp=np.array([fa_imm_elec_pagp_v2.T])
    fa_imm_elec_pagp_v2T = np.array(temp.T)
    
    temp_pagp = np.hstack(( fa_water_to_bilayer_pagp_v2T,fa_imm_elec_pagp_v2T ))       
    A_pagp = np.hstack((temp_pagp, f_elec_lipidlayer_pagp_v2T))

    #-0.72,2.49, 1.18
    experimental_ddG_pagp = np.array([-0.72,-2.44,1.64,3.32,-2.17,3.54,-2.01,-1.15,2.95,2.54,3.22,1.83,0.95,-1.75,-2.21,-1.02])
    temp=np.array([experimental_ddG_pagp.T])
    experimental_ddG_pagp_v2T = np.array(temp.T)


    #-----------------------ompla as a function of height----------------------
    ompla_height_file = outdir + "C3_OmpLA_aro_ddGs/breakdown.sc" 
    df_ompla_height = pd.read_csv(ompla_height_file, delimiter=" ")
    # list_of_mutation = ['A93F','A135F','A144F','A181F','A185F','A194F','A210F','A214F']
    # list_of_mutation.append(['A91W','A93W','A107W','A181W','A185W','A194W','A210W','A214W','A91Y','A93Y','A144Y','A166Y','A181Y','A185Y','A194Y','A210Y','A214Y'])
    # print('list of mutation: ',list_of_mutation)
    
    df_ompla_height = df_ompla_height[(df_ompla_height['mutation']=='A91I') | (df_ompla_height['mutation']=='A91M') | (df_ompla_height['mutation']=='A183M') | (df_ompla_height['mutation']=='A185M') | (df_ompla_height['mutation']=='A91V') | (df_ompla_height['mutation']=='A135M')]#|(df_ompla_height['mutation']=='A135I')|(df_ompla_height['mutation']=='A183I')|(df_ompla_height['mutation']=='A185I')|(df_ompla_height['mutation']=='A135V')|(df_ompla_height['mutation']=='A183V')|(df_ompla_height['mutation']=='A185V')]
    df_ompla_height['total_except3'] = df_ompla_height.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    # print(df_ompla_height)
    # sys.exit()
    # fa_water_to_bilayer_ompla_height = np.array([-2.105,-1.374,	-1.308,	-1.432,-1.273,-1.368])
    fa_water_to_bilayer_ompla_height = np.array(df_ompla_height['fa_water_to_bilayer:'])
    temp=np.array([fa_water_to_bilayer_ompla_height.T])
    fa_water_to_bilayer_ompla_heightT = np.array(temp.T)
    
    # total_ompla_height = np.array([-2.517,-1.693,-1.467,3.696,-3.852,-1.109,0.416,4.258,-1.901,0.091,1.223,2.525,-2.037,-0.038,-1.221,6.208,7.023,-3.962,-3.685,-1.416,0.312,3.614,-4.236,-1.736,0.095,3.348,-3.704,0.07,0.745,-0.7,-0.966,-3.63])
    # total_ompla_height = np.array([-3.704,0.07,0.745,-0.7,-0.966,-3.63])
    total_ompla_height = np.array(df_ompla_height['total_except3'])
    temp=np.array([total_ompla_height.T])
    total_ompla_heightT = np.array(temp.T)
    
    # fa_imm_elec_ompla_height = np.array([-0.171,-0.299,-0.7,-0.318,0.57,0.336,0.643,-0.191,-0.281,-0.341,-0.684,-0.428,-0.012,0.516,0.223,-0.191,-0.138,0.494,-0.199,-0.243,-0.708,-0.291,0.05,0.221,0.091,-0.171,-0.376,-0.032,-0.339,-0.028,-0.224,0.099])
    # fa_imm_elec_ompla_height = np.array([-0.376,-0.032,-0.339,-0.028,-0.224,0.099])
    fa_imm_elec_ompla_height = np.array(df_ompla_height['fa_imm_elec:'])
    temp=np.array([fa_imm_elec_ompla_height.T])
    fa_imm_elec_ompla_heightT = np.array(temp.T)

    # f_elec_lipidlayer_ompla_height = np.array([-3.824,-2.606,-2.803,-3.203,-2.497,-2.199])
    f_elec_lipidlayer_ompla_height = np.array(df_ompla_height['f_elec_lipidlayer:'])
    temp=np.array([f_elec_lipidlayer_ompla_height.T])
    f_elec_lipidlayer_ompla_heightT = np.array(temp.T)
    
    # print(total_ompla_height)
    # print(fa_imm_elec_ompla_height)
    # print(f_elec_lipidlayer_ompla_height)
    # print(fa_water_to_bilayer_ompla_height)
    # sys.exit()
    
    # experimental_ddG_ompla_height = np.array([-3.25,-1.77,-2.51,-2.65,-2.32,-1.47,-2.99,-3.33,-2.16,-1.77,-1.29,-1.93,-0.13,-2.66,-3.1,-0.98,-3.4,-0.57,-1.09,-0.01,-0.23,-0.33,-2.66,-1.47,-0.4,-2.14,-2.72,-1.24,-0.4,-0.92,-1.15,-2.1])
    experimental_ddG_ompla_height = np.array([-2.72,-1.24,-0.40,-0.92,-1.15,-2.1])
    # experimental_ddG_ompla_height = np.array([-2.72,-1.27,-4.12,-1.61,-1.24,-0.40,-0.92,-1.15,-2.1,0.16,-3.45,-1.11])
    
    temp=np.array([experimental_ddG_ompla_height.T])
    experimental_ddG_ompla_heightT = np.array(temp.T)
    
    #-----------------------ompla aromatic as a function of height----------------------
    ompla_height_file = outdir + "C3_OmpLA_aro_ddGs/breakdown.sc" 
    df_ompla_aro_height = pd.read_csv(ompla_height_file, delimiter=" ")
    list_of_mutation = ['A93F','A135F','A144F','A181F','A185F','A194F','A210F','A214F','A91W','A93W','A107W','A181W','A185W','A194W','A210W','A214W','A91Y','A93Y','A144Y','A166Y','A181Y','A185Y','A194Y','A210Y','A214Y']
    print('list of mutation: ',list_of_mutation)
    
    df_ompla_aro_height['inlist'] = df_ompla_aro_height['mutation'].isin(list_of_mutation)
    # print(df_ompla_aro_height)
    df_ompla_aro_height = df_ompla_aro_height[(df_ompla_aro_height['inlist']==True)]
    # print(df_ompla_aro_height)
    df_ompla_aro_height['total_except3'] = df_ompla_aro_height.drop(['mutation','fa_water_to_bilayer:','fa_imm_elec:','f_elec_lipidlayer:'], axis=1).sum(axis=1)
    # print(df_ompla_aro_height)
    # sys.exit()
    # fa_water_to_bilayer_ompla_height = np.array([-2.105,-1.374,	-1.308,	-1.432,-1.273,-1.368])
    fa_water_to_bilayer_ompla_aro_height = np.array(df_ompla_aro_height['fa_water_to_bilayer:'])
    temp=np.array([fa_water_to_bilayer_ompla_aro_height.T])
    fa_water_to_bilayer_ompla_aro_heightT = np.array(temp.T)
    
    # total_ompla_height = np.array([-2.517,-1.693,-1.467,3.696,-3.852,-1.109,0.416,4.258,-1.901,0.091,1.223,2.525,-2.037,-0.038,-1.221,6.208,7.023,-3.962,-3.685,-1.416,0.312,3.614,-4.236,-1.736,0.095,3.348,-3.704,0.07,0.745,-0.7,-0.966,-3.63])
    # total_ompla_height = np.array([-3.704,0.07,0.745,-0.7,-0.966,-3.63])
    total_ompla_aro_height = np.array(df_ompla_aro_height['total_except3'])
    temp=np.array([total_ompla_aro_height.T])
    total_ompla_aro_heightT = np.array(temp.T)
    
    # fa_imm_elec_ompla_height = np.array([-0.171,-0.299,-0.7,-0.318,0.57,0.336,0.643,-0.191,-0.281,-0.341,-0.684,-0.428,-0.012,0.516,0.223,-0.191,-0.138,0.494,-0.199,-0.243,-0.708,-0.291,0.05,0.221,0.091,-0.171,-0.376,-0.032,-0.339,-0.028,-0.224,0.099])
    # fa_imm_elec_ompla_height = np.array([-0.376,-0.032,-0.339,-0.028,-0.224,0.099])
    fa_imm_elec_ompla_aro_height = np.array(df_ompla_aro_height['fa_imm_elec:'])
    temp=np.array([fa_imm_elec_ompla_aro_height.T])
    fa_imm_elec_ompla_aro_heightT = np.array(temp.T)

    # f_elec_lipidlayer_ompla_height = np.array([-3.824,-2.606,-2.803,-3.203,-2.497,-2.199])
    f_elec_lipidlayer_ompla_aro_height = np.array(df_ompla_aro_height['f_elec_lipidlayer:'])
    temp=np.array([f_elec_lipidlayer_ompla_aro_height.T])
    f_elec_lipidlayer_ompla_aro_heightT = np.array(temp.T)
    
    experimental_ddG_ompla_aro_height = np.array([-3.25,-1.75,-2.51,-2.65,-2.32,-1.47,-2.99,-3.33,-2.16,-2.62,-1.77,-0.13,-2.66,-3.1,-0.98,-3.4,-0.57,-1.09,-3.14,-0.23,-0.33,-2.66,-1.47,-0.4,-2.14])
    #=================================================================================
    
    temp_ompla_height = np.hstack(( fa_water_to_bilayer_ompla_heightT,fa_imm_elec_ompla_heightT ))       
    A_ompla_height = np.hstack((temp_ompla_height, f_elec_lipidlayer_ompla_heightT))
    
    B_ompla_height = np.array( experimental_ddG_ompla_heightT - total_ompla_heightT )

    #coefficients = np.array([0.0,0.0,0.0])
    for i in range(len(experimental_ddG_pagp)-2):
        if(experimental_ddG_ompla[i]==0):
            print(i)
            continue
        A = np.array([[fa_water_to_bilayer_ompla_v2[i], fa_imm_elec_ompla_v2[i], f_elec_lipidlayer_ompla_v2[i]],
                      [fa_water_to_bilayer_ompla_v2[i+1], fa_imm_elec_ompla_v2[i+1], f_elec_lipidlayer_ompla_v2[i+1]],
                      [fa_water_to_bilayer_ompla_v2[i+2], fa_imm_elec_ompla_v2[i+2], f_elec_lipidlayer_ompla_v2[i+2]]])
       #               [fa_water_to_bilayer[i+3], fa_imm_elec[i+3], f_elec_lipidlayer[i+3]]])
        B = np.array([experimental_ddG_ompla[i]-total_ompla_v2[i], experimental_ddG_ompla[i+1]-total_ompla_v2[i+1], experimental_ddG_ompla[i+2]-total_ompla_v2[i+2]])
    
        #print(A)
        #print(B)
        C = np.linalg.solve(A,B)
        if(i==0):
            coefficients = C 
        else:
             coefficients = np.vstack((coefficients, C))
    df = pd. DataFrame(coefficients, columns=['fa_water_bilayer', 'fa_imm_elec','fa_elec_bilayer'])
    print(df)
   
    fig, (axs1,axs2,axs3) = plt.subplots(3)
    axs1.hist(coefficients[:,0], bins=20, density=False)
    axs1.set_title('fa_water_to_bilayer')
    # axs1.xlabel('bins')
    # axs1.xlabel('#cases')
    counts, binEdges=np.histogram( coefficients[:,0],bins=20,density=False )
    temp=np.array([counts.T])
    counts = np.array(temp.T)
    temp=np.array([binEdges.T])
    binEdges = np.array(temp.T)
    
    hist = np.hstack((binEdges[0:len(binEdges)-1],counts))
    hist_A = pd.DataFrame(hist, columns=['Bins','fa_water_to_bilayer'])
    #print(hist_A)
    hist_A.to_csv('Ompla_coefficient1.csv',index=False)

    axs2.hist(coefficients[:,1], bins=50, density=False)
    # axs2.xlabel('bins')
    # axs2.xlabel('#cases')
    axs2.set_title('fa_imm_elec')
    counts,binEdges=np.histogram( coefficients[:,1],bins=50,density=False )
    temp=np.array([counts.T])
    counts = np.array(temp.T)
    temp=np.array([binEdges.T])
    binEdges = np.array(temp.T)

    hist = np.hstack((binEdges[0:len(binEdges)-1],counts))
    hist_B = pd.DataFrame(hist, columns=['Bins','fa_imm_elec'])
    # hist_B.to_csv('Ompla_coefficient2.csv',index=False)
    print(hist_B)
    #print(binEdges)

    axs3.hist(coefficients[:,2], bins=20, density=False)
    axs3.set_title('f_elec_lipidlayer')
    # axs3.xlabel('bins')
    # axs3.xlabel('#cases')
    counts, binEdges=np.histogram( coefficients[:,2],bins=20,density=False )
    temp=np.array([counts.T])
    counts = np.array(temp.T)
    temp=np.array([binEdges.T])
    binEdges = np.array(temp.T)
    
    hist = np.hstack((binEdges[0:len(binEdges)-1],counts))
    hist_C = pd.DataFrame(hist, columns=['Bins','f_elec_lipidlayer'])
    hist_C.to_csv('Ompla_coefficient3.csv',index=False)

    plt.savefig('PagP_coefficients.png')

    print('-----------------OmpLA----------------')
    #print(A_ompla)
    B_ompla = np.array( experimental_ddG_omplaT - total_ompla_v2T )
    #print(B_ompla)
    # reg_ompla = LinearRegression(fit_intercept=False).fit(A_ompla, B_ompla)
    # print(reg_ompla)
    # print(reg_ompla.coef_)
    # print(reg_ompla.intercept_)
    # print(reg_ompla.score(A_ompla, B_ompla, sample_weight=None))


    print('------------PagP--------------------')
    #print(A_pagp)
    B_pagp = np.array( experimental_ddG_pagp_v2T - total_pagp_v2T )
    #print(B_pagp)

    # reg_pagp = LinearRegression(fit_intercept=False).fit(A_pagp, B_pagp)
    # print(reg_pagp)
    # print(reg_pagp.coef_)
    # print(reg_pagp.intercept_)
    # print(reg_pagp.score(A_pagp, B_pagp, sample_weight=None))

    #coefficient_final = pd.DataFrame(coefficients, columns = ['fa_water_to_bilayer','fa_imm_elec','fa_elec_lipidlayer'])
    #print(coefficient_final)
    #coefficient_final.to_csv('Ompla_coefficient.csv',index=False)
    
    #hist = coefficient_final['fa_water_to_bilayer', 'fa_imm_elec', 'fa_elec_lipidlayer'].hist()
    #hist_faie = coefficient_final['fa_imm_elec'].hist()
    #hist_fael = coefficient_final['fa_elec_lipidlayer'].hist()
    #print(hist)
    #plt.show()
    #plt.savefig('hist_ompla_coefficients.png', bbox_inches = 'tight')
    A = np.vstack((A_ompla,A_pagp))
    A_combined = np.vstack((A, A_ompla_height))
    
    B = np.vstack((B_ompla,B_pagp))
    B_combined = np.vstack((B, B_ompla_height))

    del(A,B)
    # print("====================")
    # print(A_combined)
    # print(B_combined)
    # print("=====================")
    # reg_combined = LinearRegression(fit_intercept=False).fit(A_combined, B_combined)
    # print(reg_combined)
    # print(reg_combined.coef_)
    # print(reg_combined.intercept_)
    # print(reg_combined.score(A_pagp, B_pagp, sample_weight=None))


    #one crazy experiment
    #i am combining equations from both ompla, pagP
    # fa_water_to_bilayer_combined = fa_water_to_bilayer_ompla_v2
    # fa_water_to_bilayer_combined = fa_water_to_bilayer_ompla_height
    # fa_water_to_bilayer_combined = np.concatenate((fa_water_to_bilayer_ompla_v2, fa_water_to_bilayer_pagp_v2), axis =0)#np.array(fa_water_to_bilayer_pagp_v2)#
    # fa_water_to_bilayer_combined = np.concatenate((fa_water_to_bilayer_combined, fa_water_to_bilayer_ompla_height), axis =0)
    fa_water_to_bilayer_combined = np.concatenate((fa_water_to_bilayer_ompla_v2, fa_water_to_bilayer_ompla_aro_height), axis =0)#np.array(fa_water_to_bilayer_pagp_v2)#
    # fa_water_to_bilayer_combined = np.concatenate((fa_water_to_bilayer_combined, fa_water_to_bilayer_ompla_aro_height), axis =0)#np.array(fa_water_to_bilayer_pagp_v2)#
    fa_water_to_bilayer_ompla_center_combined = np.array(fa_water_to_bilayer_ompla_v2)
    fa_water_to_bilayer_pagp_center_combined = np.array(fa_water_to_bilayer_pagp_v2)
    fa_water_to_bilayer_ompla_height_combined = np.concatenate((fa_water_to_bilayer_ompla_height,fa_water_to_bilayer_ompla_aro_height), axis =0)
    
    
    # total_combined = total_pagp_v2
    # total_combined = total_ompla_height
    # total_combined = np.concatenate((total_ompla_v2, total_pagp_v2), axis = 0)#np.array(total_pagp_v2)#
    # total_combined = np.concatenate((total_combined, total_ompla_height), axis = 0)
    total_combined = np.concatenate((total_ompla_v2, total_ompla_aro_height), axis = 0)
    # total_combined = np.concatenate((total_combined, total_ompla_aro_height), axis = 0)
    total_ompla_center_combined = np.array(total_ompla_v2)
    total_pagp_center_combined = np.array(total_pagp_v2)
    total_ompla_height_combined = np.concatenate((total_ompla_height,total_ompla_aro_height), axis =0)
    
    # fa_imm_elec_combined = fa_imm_elec_pagp_v2
    # fa_imm_elec_combined = fa_imm_elec_ompla_height
    # fa_imm_elec_combined = np.concatenate((fa_imm_elec_ompla_v2, fa_imm_elec_pagp_v2), axis = 0)#np.array(fa_imm_elec_pagp_v2)#
    # fa_imm_elec_combined = np.concatenate((fa_imm_elec_combined, fa_imm_elec_ompla_height), axis = 0)#np.array(fa_imm_elec_pagp_v2)#
    fa_imm_elec_combined = np.concatenate((fa_imm_elec_ompla_v2, fa_imm_elec_ompla_aro_height), axis = 0)
    # fa_imm_elec_combined = np.concatenate((fa_imm_elec_combined, fa_imm_elec_ompla_aro_height), axis = 0)
    fa_imm_elec_ompla_center_combined = np.array(fa_imm_elec_ompla_v2)
    fa_imm_elec_pagp_center_combined = np.array(fa_imm_elec_pagp_v2)
    fa_imm_elec_ompla_height_combined = np.concatenate((fa_imm_elec_ompla_height,fa_imm_elec_ompla_aro_height), axis =0)
    
    
    # f_elec_lipidlayer_combined = f_elec_lipidlayer_pagp_v2
    # f_elec_lipidlayer_combined = f_elec_lipidlayer_ompla_height
    # f_elec_lipidlayer_combined = np.concatenate((f_elec_lipidlayer_ompla_v2, f_elec_lipidlayer_pagp_v2,), axis = 0)#np.array(f_elec_lipidlayer_pagp_v2)#
    # f_elec_lipidlayer_combined = np.concatenate((f_elec_lipidlayer_combined, f_elec_lipidlayer_ompla_height,), axis = 0)#np.array(f_elec_lipidlayer_pagp_v2)#
    f_elec_lipidlayer_combined = np.concatenate((f_elec_lipidlayer_ompla_v2, f_elec_lipidlayer_ompla_aro_height), axis = 0)
    # f_elec_lipidlayer_combined = np.concatenate((f_elec_lipidlayer_combined, f_elec_lipidlayer_ompla_aro_height), axis = 0)
    f_elec_lipidlayer_ompla_center_combined = np.array(f_elec_lipidlayer_ompla_v2)
    f_elec_lipidlayer_pagp_center_combined = np.array(f_elec_lipidlayer_pagp_v2)
    f_elec_lipidlayer_ompla_height_combined = np.concatenate((f_elec_lipidlayer_ompla_height,f_elec_lipidlayer_ompla_aro_height), axis =0)
    
    
    # experimental_combined = experimental_ddG_pagp
    # experimental_combined = experimental_ddG_ompla_height
    # experimental_combined = np.concatenate((experimental_ddG_ompla, experimental_ddG_pagp), axis =0)#np.array(experimental_ddG_pagp)#
    # experimental_combined = np.concatenate((experimental_combined, experimental_ddG_ompla_height), axis =0)#np.array(experimental_ddG_pagp)#
    experimental_combined = np.concatenate((experimental_ddG_ompla, experimental_ddG_ompla_aro_height), axis =0)
    # experimental_combined = np.concatenate((experimental_combined, experimental_ddG_ompla_aro_height), axis =0)
    experimental_ompla_center_combined = np.array(experimental_ddG_ompla)
    experimental_pagp_center_combined = np.array(experimental_ddG_pagp)
    experimental_ompla_height_combined = np.concatenate((experimental_ddG_ompla_height,experimental_ddG_ompla_aro_height), axis =0)
    
    # calculate_positive_lstsq_sol(fa_water_to_bilayer_combined,fa_imm_elec_combined,f_elec_lipidlayer_combined,total_combined,experimental_combined)
    # calculate_positive_lstsq_sol_with2var(fa_water_to_bilayer_combined,fa_imm_elec_combined,f_elec_lipidlayer_combined,total_combined,experimental_combined)
    
    print('----------ompla_height only---------------------')
    guess_afterfit = calculate_positivevar_maxr2(fa_water_to_bilayer_combined,fa_imm_elec_combined,f_elec_lipidlayer_combined,total_combined,experimental_combined)
    # guess_afterfit = calculate_positivevar_maxr2(fa_water_to_bilayer_ompla_height_combined,fa_imm_elec_ompla_height_combined,f_elec_lipidlayer_ompla_height_combined,total_ompla_height_combined,experimental_ompla_height_combined)
    # guess_afterfit = calculate_positivevar_maxr2(fa_water_to_bilayer_ompla_center_combined,fa_imm_elec_ompla_center_combined,f_elec_lipidlayer_ompla_center_combined,total_ompla_center_combined,experimental_ompla_center_combined)
    
    print(guess_afterfit)
    print('----------ompla_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_ompla_center_combined,fa_imm_elec_ompla_center_combined,f_elec_lipidlayer_ompla_center_combined,total_ompla_center_combined,experimental_ompla_center_combined,guess_afterfit)
    # calculate_positivevar_maxr2(fa_water_to_bilayer_ompla_center_combined,fa_imm_elec_ompla_center_combined,f_elec_lipidlayer_ompla_center_combined,total_ompla_center_combined,experimental_ompla_center_combined)
    
    print('----------pagp_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_pagp_center_combined,fa_imm_elec_pagp_center_combined,f_elec_lipidlayer_pagp_center_combined,total_pagp_center_combined,experimental_pagp_center_combined,guess_afterfit)
    # calculate_positivevar_maxr2(fa_water_to_bilayer_pagp_center_combined,fa_imm_elec_pagp_center_combined,f_elec_lipidlayer_pagp_center_combined,total_pagp_center_combined,experimental_pagp_center_combined)
    
    print('----------ompla_height_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_ompla_height_combined,fa_imm_elec_ompla_height_combined,f_elec_lipidlayer_ompla_height_combined,total_ompla_height_combined,experimental_ompla_height_combined,guess_afterfit)
    # calculate_positivevar_maxr2(fa_water_to_bilayer_ompla_height_combined,fa_imm_elec_ompla_height_combined,f_elec_lipidlayer_ompla_height_combined,total_ompla_height_combined,experimental_ompla_height_combined)
    
    guess_afterfit = np.array([1.5,0.0,0.0])
    print('----------ompla_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_ompla_center_combined,fa_imm_elec_ompla_center_combined,f_elec_lipidlayer_ompla_center_combined,total_ompla_center_combined,experimental_ompla_center_combined,guess_afterfit)
    
    print('----------pagp_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_pagp_center_combined,fa_imm_elec_pagp_center_combined,f_elec_lipidlayer_pagp_center_combined,total_pagp_center_combined,experimental_pagp_center_combined,guess_afterfit)
    
    print('----------ompla_height_only-----------------------------')
    calculate_rsq_with_guess(fa_water_to_bilayer_ompla_height_combined,fa_imm_elec_ompla_height_combined,f_elec_lipidlayer_ompla_height_combined,total_ompla_height_combined,experimental_ompla_height_combined,guess_afterfit)
    
    # print("===========")
    # result = pd.DataFrame()
    # print("fa_water_bilayer fa_imm_elec fa_elec_bilayer abs_error_sum corr_coeff corr_ompla corr_pagp\n")
    # for iter in range(1000):
    #         value = np.arange(len(experimental_combined))
    #         np.random.shuffle(value)
    #         # print(value)
    #         #forcefully input values
    #         #value=np.array([ 5,6,16,10,9,8,4,12,26,29,20,21,7,24,23,11,0,13,30,25,2,27,18,1,14,31,28,15,17,22,3,19])
    #         #print(value)
    #         for i in range(len(value)-2):
    #             if(experimental_combined[value[i]]==0):
    #                 print(i)
    #                 continue
    #             A = np.array([[fa_water_to_bilayer_combined[value[i]], fa_imm_elec_combined[value[i]], f_elec_lipidlayer_combined[value[i]]],
    #                         [fa_water_to_bilayer_combined[value[i+1]], fa_imm_elec_combined[value[i+1]], f_elec_lipidlayer_combined[value[i+1]]],
    #                         [fa_water_to_bilayer_combined[value[i+2]], fa_imm_elec_combined[value[i+2]], f_elec_lipidlayer_combined[value[i+2]]]])
    #         #               [fa_water_to_bilayer[i+3], fa_imm_elec[i+3], f_elec_lipidlayer[i+3]]])
    #             B = np.array([experimental_combined[value[i]]-total_combined[value[i]], experimental_combined[value[i+1]]-total_combined[value[i+1]], experimental_combined[value[i+2]]-total_combined[value[i+2]]])
            
    #             #print(A)
    #             #print(B)
    #             C = np.linalg.solve(A,B)

    #             if(i==0):
    #                 coefficients = C 
    #             else:
    #                 coefficients = np.vstack((coefficients, C))
    #         df = pd. DataFrame(coefficients, columns=['fa_water_bilayer', 'fa_imm_elec','fa_elec_bilayer'])
        
    #         # print('------------------printing the coefficient matrix------------------\n')
    #         A_min = np.stack((fa_water_to_bilayer_combined,fa_imm_elec_combined,f_elec_lipidlayer_combined), axis = 1)#np.array([fa_water_to_bilayer_combined[:],fa_imm_elec_combined[:],f_elec_lipidlayer_combined[:]])
    #         B_min = np.array(experimental_combined-total_combined)
            
    #         A_min_ompla = np.stack((fa_water_to_bilayer_ompla_v2,fa_imm_elec_ompla_v2,f_elec_lipidlayer_ompla_v2), axis = 1)#np.array([fa_water_to_bilayer_combined[:],fa_imm_elec_combined[:],f_elec_lipidlayer_combined[:]])
    #         B_min_ompla = np.array(experimental_ddG_ompla-total_ompla_v2)
            
    #         A_min_pagp = np.stack((fa_water_to_bilayer_pagp_v2,fa_imm_elec_pagp_v2,f_elec_lipidlayer_pagp_v2), axis = 1)#np.array([fa_water_to_bilayer_combined[:],fa_imm_elec_combined[:],f_elec_lipidlayer_combined[:]])
    #         B_min_pagp = np.array(experimental_ddG_pagp-total_pagp_v2)
            
    #         # guess = np.array([0.50,1.0,0.50])
    #         error = np.zeros(len(df))
    #         rsq_overall = np.zeros(len(df))
    #         rsq_ompla = np.zeros(len(df))
    #         rsq_pagp = np.zeros(len(df))
    #         for i in range(len(error)):
    #             guess = np.array([df.iat[i,0], df.iat[i,1], df.iat[i,2]])
    #             f = objective_func(A_min, B_min, guess)
    #             error[i] = f
    #             rsq_overall[i] = calculate_regression_coeff(A_min,B_min,guess)
    #             rsq_ompla[i] = calculate_regression_coeff(A_min_ompla,B_min_ompla,guess)
    #             rsq_pagp[i] = calculate_regression_coeff(A_min_pagp,B_min_pagp,guess)

    #         # bounds = Bounds([0.0,0.0,0.0], [np.inf, np.inf, np.inf])
    #         df['abs_error_sum'] = error
    #         df['corr_coeff'] = rsq_overall
    #         df['corr_ompla'] = rsq_ompla
    #         df['corr_pagp'] = rsq_pagp
    #         # print(df.iloc[np.where(error == np.amin(error))])
            
    #         for i in range(df.shape[0]):
    #             if(df.iat[i,0]>= 0 and df.iat[i,1]>= 0 and df.iat[i,2]>= 0 and df.iat[i,5].round(3)>0.50):
    #                 print(df.iat[i,0].round(3), df.iat[i,1].round(3), df.iat[i,2].round(3), df.iat[i,3].round(3), df.iat[i,4].round(3),df.iat[i,5].round(3),df.iat[i,6].round(3))
            
            
    #         df = df[(df['fa_water_bilayer']>=0) & (df['fa_imm_elec']>=0) & (df['fa_elec_bilayer']>=0)]
    #         labelout = outdir+'coefficient_f2021_wts_basedonompla.dat'
    #         df.to_csv(labelout, mode='a', sep=" ", index=False)       
            
    # TO Do:
    #read all files, and find coefficients which give min. error, max coeff. 
    #use code for 3 files. 
    #then submit the jobs ddg and tilt with those. 
def hist_csv(array, filename, arrayname, binsize):
    
    counts, binEdges=np.histogram( array,bins=binsize,density=False )
    temp=np.array([counts.T])
    counts = np.array(temp.T)
    temp=np.array([binEdges.T])
    binEdges = np.array(temp.T)
    
    hist = np.hstack((binEdges[0:len(binEdges)-1],counts))
    hist_A = pd.DataFrame(hist, columns=['Bins',arrayname])
    #print(hist_A)
    hist_A.to_csv(filename,index=False)
    color = ['#fafa6e','#23aa8f','#2a4858']
    fig, ax = plt.subplots(figsize=(10,10))
    
    if(arrayname =='fa_water_bilayer'):
            ax = plt.hist(array, density=True, bins=binsize, histtype ='bar', rwidth = 1.0, color = color[0])  # density=False would make counts
    elif(arrayname=='fa_imm_elec'):
            ax = plt.hist(array, density=True, bins=binsize, histtype ='bar', rwidth = 2.0, color = color[1])  # density=False would make counts
    else:
            ax = plt.hist(array, density=True, bins=binsize, histtype ='bar', rwidth = 1.0, color = color[2])  # density=False would make counts

    axis_font = {'fontname':'Arial', 'size':'50'}
    plt.xticks(rotation=0, **axis_font)
    plt.yticks(rotation=0, **axis_font)
    plt.tick_params(width=1)
    # ax.set_aspect("equal")
    # ax.axvline(x=0, color='k',linewidth=1)
    plt.ylabel('Probability', **axis_font)
    plt.xlabel(arrayname, **axis_font)
    # ax.get_frame().set_linewidth(1.0)
    plt.savefig( filename+ ".pdf", bbox_inches="tight", transparent=True, dpi=300, format="pdf" )

    
# def calculate_error(A, B):
#     error = 0.0

#     C = np.array([A.T])*np.array([p.T])
#     temp = np.array(C-B)
#     error = np.sum(temp[:])
#     print(error)
#     return(error)

def objective_func(A_min, B_min, x):
    f = np.sum(abs(np.matmul(A_min,x.T)-B_min.T))
    return(f)

def calculate_regression_coeff(A,B,x):
    #we are trying to calculate the regression coefficient for the 
    #solution x, where Ax=B
    from scipy import stats
    
    X = np.matmul(A,x.T)
    # print(X)
    # print(B)
    slope, intercept, r_value, p_value, std_err = stats.linregress(X,B)

    # print(reg.coef_)
    #coeffciient of regression
    return(r_value*r_value)

#     # for i in range(len(B)):
#     #     error = error + abs(A[:,i]*)

def calculate_positive_lstsq_sol(A1,A2,A3, B1, B2):
    from scipy.optimize import nnls
    from scipy.optimize import lsq_linear
    A=np.vstack([A1,A2,A3]).T 
    B=(B2-B1)
    
    # print(A)
    # print(B)
    x, rnorm = nnls(A,B)
    print('nnls sol: fa_water_to_bilayer:{} fa_imm_elec:{} f_elec_lipid:{} error:{}'.format(round(x[0],3),round(x[1],3),round(x[2],3),round(rnorm,3)))
    
    print(lsq_linear(A, B, bounds=(0.001, np.inf)))
    
def calculate_positivevar_maxr2(A1,A2,A3,B1,B2):
    from scipy import optimize
    A = np.vstack([A1,A2,A3]).T 
    B = (B2-B1)
    # print(A)
    # print(B)
    # sys.exit()
    guess = np.array([0.10,0.10,0.10])
    # b1 = [(0.01,np.inf),(1.0,1.0),(0.01,np.inf)]
    b1 = [(0.01,1.0),(0.01,1.0),(0.01,1.0)]
    
    def rsq_func(x0):
        r_sq = calculate_regression_coeff(A,B,x0)
        return(1.0 - r_sq)
    
    result = optimize.minimize(fun=rsq_func, bounds=b1,x0=guess)
    print("starting guess is:", guess)
    print("bound is:", b1)
    print("rsq is :{}".format(1.0 - rsq_func(guess)))
    print(result.x)
    print("rsq after minimization is :{}".format(1.0 - rsq_func(result.x)))
    return(result.x)

def calculate_rsq_with_guess(A1,A2,A3,B1,B2,guess):
    A = np.vstack([A1,A2,A3]).T 
    B = (B2-B1)
    r_sq = calculate_regression_coeff(A,B,guess)
    print('rsq for the guess:{} is {}'.format(guess,r_sq))
    
def calculate_positive2var_maxr2(A1,A2,A3,B1,B2):
    from scipy import optimize
    A = np.vstack([A2,A3]).T 
    print(A2)
    print(A3)
    B = (B2-B1-0.5*A1)
    guess = np.array([0.5,0.50])
    b1 = [(0.1,np.inf),(0.001,np.inf)]
    
    def rsq_func(x0):
        r_sq = calculate_regression_coeff(A,B,x0)
        return(1.0 - r_sq)
    
    result = optimize.minimize(fun=rsq_func, bounds=b1,x0=guess)
    print("rsq is :{}".format(rsq_func(guess)))
    print(result) 
       
def calculate_positive_lstsq_sol_with2var(A1,A2,A3, B1, B2):
    
    from scipy.optimize import lsq_linear
    A=np.vstack([A2,A3]).T 
    B=(B2-B1-0.5*A1)
    
    print(A)
    print(B)
    
    x, rnorm = nnls(A,B)
    print('nnls sol: fa_water_to_bilayer:{} fa_imm_elec:{} f_elec_lipid:{} error:{}'.format(round(0.500,3),round(x[0],3),round(x[1],3),round(rnorm,3)))
    x = lsq_linear(A, B, bounds=(0.001, np.inf))
    print(lsq_linear(A, B).x)#, bounds=(0.001, np.inf)))
    # print('unbound lstsq  sol: fa_water_to_bilayer:{} fa_imm_elec:{} f_elec_lipid:{} error:{}'.format(round(0.500,3),round(x[0],3),round(x[1],3)))

def calculate_dg_rosetta_atom_type():
    from scipy.optimize import lsq_linear
    from scipy import optimize
    
    atom_types = {'CH0':1,'CH2':2,'CH1':3,'CH3':4,
				'CNH2':5,'COO':6,'S':7,'SH1':8, 
				'aroC':9,'Nhis':10,'Ntrp':11,'Nlys':12,
				'NH2O':13,'Narg':14,'NtrR':15,
				'ONH2':16,'OOC':17,'OH':18}
    amino_acids = {'C':1,'D-1':2,'E-1':3,
				'F':4,'H':5,'I':6,'K':7, 
				'L':8,'M':9,'N':10,'Q':11,
				'R':12,'S':13,
				'T':14,'V':15,'W':16,'Y':17, 'D0':18,'E0':19}
    A = np.zeros((19,18))
    # Aalternate = np.zeros(17,18)
    B = np.zeros(19)
    #B is delG_w,l|sc from Moon and Fleming paper. This is not delg|x - delg|A. 
    for index in range(1,len(B)+1):
        print(index)
        if(index==amino_acids['C']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['SH1']-1] = 1.0
            B[index-1] = -1.08
            # Aalternate[index-1][amino_acids['CH2']-1] = 1.0
        elif(index==amino_acids['D-1']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['COO']-1] = 1.0
            A[index-1][atom_types['OOC']-1] = 2.0
            if(A[index-1][atom_types['OOC']-1]==1):
                A[index-1][atom_types['OH']-1] = 1.0
                B[index-1] = 1.38
            else:
                B[index-1] = 6.31
        elif(index==amino_acids['E-1']):
            A[index-1][atom_types['CH2']-1] = 2.0
            A[index-1][atom_types['COO']-1] = 1.0
            A[index-1][atom_types['OOC']-1] = 2.0
            # A[index-1][atom_types['OH']-1] = 1.0
            if(A[index-1][atom_types['OOC']-1]==1.0):
                A[index-1][atom_types['OH']-1] = 1.0
                B[index-1] = 0.07
            else:
                B[index-1] = 6.31
            # B[index-1] = 0.07
        elif(index==amino_acids['F']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH0']-1] = 1.0
            A[index-1][atom_types['aroC']-1] = 5.0
            B[index-1] = -3.77
        elif(index==amino_acids['H']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH0']-1] = 1.0
            A[index-1][atom_types['Nhis']-1] = 1.0
            A[index-1][atom_types['aroC']-1] = 2.0
            A[index-1][atom_types['Ntrp']-1] = 1.0
            B[index-1] = 3.19
        elif(index==amino_acids['I']):
            A[index-1][atom_types['CH3']-1] = 2.0
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH1']-1] = 1.0
            B[index-1] = -3.12 
        elif(index==amino_acids['K']):
            A[index-1][atom_types['CH2']-1] = 4.0
            A[index-1][atom_types['Nlys']-1] = 1.0
            B[index-1] = 3.82
        elif(index==amino_acids['L']):
            A[index-1][atom_types['CH3']-1] = 2.0
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH1']-1] = 1.0
            B[index-1] = -3.32
        elif(index==amino_acids['M']):
            A[index-1][atom_types['CH3']-1] = 1.0
            A[index-1][atom_types['CH2']-1] = 2.0
            A[index-1][atom_types['S']-1] = 1.0
            B[index-1] = -2.33
        elif(index==amino_acids['N']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['ONH2']-1] = 1.0
            A[index-1][atom_types['CNH2']-1] = 1.0
            A[index-1][atom_types['NH2O']-1] = 1.0
            B[index-1] = 1.91
        elif(index==amino_acids['Q']):
            A[index-1][atom_types['CH2']-1] = 2.0
            A[index-1][atom_types['ONH2']-1] = 1.0
            A[index-1][atom_types['CNH2']-1] = 1.0
            A[index-1][atom_types['NH2O']-1] = 1.0
            B[index-1] = 1.44
        elif(index==amino_acids['R']):
            A[index-1][atom_types['CH2']-1] = 3.0
            A[index-1][atom_types['aroC']-1] = 1.0
            A[index-1][atom_types['Narg']-1] = 2.0
            A[index-1][atom_types['NtrR']-1] = 1.0
            B[index-1] = 2.14
        elif(index==amino_acids['S']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['OH']-1] = 1.0
            B[index-1] = 0.26
        elif(index==amino_acids['T']):
            A[index-1][atom_types['CH3']-1] = 1.0
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['OH']-1] = 1.0
            B[index-1] = 0.21
        elif(index==amino_acids['V']):
            A[index-1][atom_types['CH3']-1] = 2.0
            A[index-1][atom_types['CH1']-1] = 1.0
            B[index-1] = -2.34
        elif(index==amino_acids['W']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH0']-1] = 3.0
            A[index-1][atom_types['aroC']-1] = 5.0
            A[index-1][atom_types['Ntrp']-1] = 1.0
            B[index-1] = -1.95
        elif(index==amino_acids['Y']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['CH0']-1] = 2.0
            A[index-1][atom_types['aroC']-1] = 4.0
            A[index-1][atom_types['OH']-1] = 1.0
            B[index-1] = -2.66
        elif(index==amino_acids['D0']):
            A[index-1][atom_types['CH2']-1] = 1.0
            A[index-1][atom_types['COO']-1] = 1.0
            A[index-1][atom_types['OOC']-1] = 1.0
            if(A[index-1][atom_types['OOC']-1]==1.0):
                A[index-1][atom_types['OH']-1] = 1.0
                B[index-1] = 1.38
            else:
                B[index-1] = 6.31
        elif(index==amino_acids['E0']):
            A[index-1][atom_types['CH2']-1] = 2.0
            A[index-1][atom_types['COO']-1] = 1.0
            A[index-1][atom_types['OOC']-1] = 1.0
            # A[index-1][atom_types['OH']-1] = 1.0
            if(A[index-1][atom_types['OOC']-1]==1.0):
                A[index-1][atom_types['OH']-1] = 1.0
                B[index-1] = 0.07
            else:
                B[index-1] = 6.31
            # B[index-1] = 0.07
        
    print(A)
    print(B)   
    sol = lsq_linear(A,B, bounds=(-1.50,7.0)).x
    print(lsq_linear(A,B).x)
     
    
    def rsq_func(x0):
        r_sq = calculate_regression_coeff(A,B,x0)
        return(1.0 - r_sq)
    
    guess = np.array([-0.76,-0.74,-0.69,-0.885,0.9,0.3,0.03,-0.25,-0.45,2.25,3.3,6.7,0.93,1.92,0.96,0.93,0.73,1.41])
    guess2 = np.array([-0.76,-0.74,-0.69,-0.885,0.9,-1.5,0.03,-0.25,-0.45,2.25,3.3,6.7,0.93,1.92,0.96,0.93,4.125,1.41])
    dg_imm1 = np.array([-0.3,-0.791,-0.3,-0.3,0,-0.211,0.612,0.612,-0.378,4.147,3.102,9.305,2.773,3.311,3.102,2.773,4.206,3.3])
    # guess = np.ones(18)
    result = optimize.minimize(fun=rsq_func, x0=guess)
    print(result.x)
    print("rsq with guess is :{}".format(1 - rsq_func(guess)))
    print("rsq with guess2 is :{}".format(1 - rsq_func(guess2)))
    print("rsq with r2 optimized result is: {}".format(1-rsq_func(result.x)))
    print("rsq with error optimized result is: {}".format(1-rsq_func(sol)))
    print('atom dg_reb dg_lsq dg_min dg_imm1 change_coo_ooc')
    for atoms in atom_types:
        # print(atoms)
        # print(sol[atom_types[atoms]-1])
        print('{} {} {} {} {} {}'.format(atoms, guess[atom_types[atoms]-1], round(sol[atom_types[atoms]-1],3), round(result.x[atom_types[atoms]-1],3), dg_imm1[atom_types[atoms]-1], guess2[atom_types[atoms]-1]))   

    guess_dg = np.matmul(A,guess.T)
    guess2_dg = np.matmul(A,guess2.T)
    imm1_dg = np.matmul(A,dg_imm1.T)
    sol_dg = np.matmul(A,sol.T)
    result_dg = np.matmul(A,result.x.T)
    print('AA dG_MF dG_reb dG_lsq dG_min dG_imm1 change_coo_ooc')
    for residue in amino_acids:
        print('{} {} {} {} {} {} {}'.format(residue, B[amino_acids[residue]-1], round(guess_dg[amino_acids[residue]-1],3), round(sol_dg[amino_acids[residue]-1],3), round(result_dg[amino_acids[residue]-1],3), round(imm1_dg[amino_acids[residue]-1],3), round(guess2_dg[amino_acids[residue]-1],3)))   

if __name__ == "__main__": main()
