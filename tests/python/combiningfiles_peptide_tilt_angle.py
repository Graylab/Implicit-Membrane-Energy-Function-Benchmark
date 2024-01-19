#!/usr/bin/env python
""" Master script for combininig all files into one file
and generating a energy landscape as a function of depth and tilt 
angle  minimized  over rotation (/azimuthal) angle. 
This module is used for tests 1,2,3,5 and 6. 

Authors: 
	Rituparna Samanta <rsamanta@utexas.edu> 

Example: 
	$ python3 combiningfiles.py --which_tests "tm-peptide-tilt-angle" --energy_fxn "franklin2019"

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
	- Rosetta release 246 or greater
	- PyRosetta4 for Python 3.6 or 3.7
"""

import pwd
import sys
import os
import read_config
import glob
import math
import fileinput
import pandas as pd
import numpy as np
from tqdm import tqdm
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)


def combiningallfiles(peptide_list, energy_fxn, pH=[], tilt_axis=False, testname=[]):

    for element_num in tqdm(range(len(peptide_list))):

        # start_index = ['0.000000', '20.000000','30.000000', '40.000000', '50.000000']
        # end_index = ['10.000000', '30.000000','40.000000', '50.000000', '60.000000']  # [-3, -2, -1, 0, 0, 0, 1, 2]
        if(tilt_axis):

            start_index = ['-50.000000', '-40.000000','-30.000000', '-20.000000', '-10.000000', '0.000000', '10.000000', '20.000000','30.000000', '40.000000']#, '50.000000']
            end_index = ['-40.000000', '-30.000000','-20.000000', '-10.000000', '1.000000', '10.000000', '20.000000', '30.000000','40.000000', '50.000000']#, '60.000000']  # [-3, -2, -1, 0, 0, 0, 1, 2]
        else:
            # start_index = ['-60.000000','-50.000000', '-40.000000','-30.000000', '-20.000000', '-10.000000']#, '10.000000', '20.000000','30.000000', '40.000000', '50.000000']
            # end_index = ['-50.000000','-40.000000', '-30.000000','-20.000000', '-10.000000','1.000000']
            start_index = ['-45.000000']#,'-50.000000', '-40.000000','-30.000000', '-20.000000', '-10.000000']#, '10.000000', '20.000000','30.000000', '40.000000', '50.000000']
            end_index = ['1.000000']#,'-40.000000', '-30.000000','-20.000000', '-10.000000','1.000000']
        
        # substring = "_after_relax"
        # substring = "_0001"
        substring=""
        
        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
                
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt" or testname=="amphiphilic-peptide-tilt-angle"):
            peptide_name = peptide_name2
            
            
            if(testname=="amphiphilic-peptide-tilt-angle"):
                if(peptide_name !='3tij_h1'):
                    substring = "_relaxed_002"
                else:
                    substring = "_relaxed_001"
        
        df=pd.DataFrame()
        if(pH == []):
                label_in = peptide_list[element_num] + '/' + \
                    peptide_list[element_num] + '_raw.dat'
                labelin1 = peptide_list[element_num] + \
                    '/' + peptide_list[element_num] + '_combined.dat'
                labelout1 = peptide_list[element_num] + '/' + peptide_list[element_num] + '_min.dat'
                print(labelin1)
                print(labelout1)

        else:
            label_in = peptide_list[element_num] + '_' + str(pH) + '/' + \
                peptide_list[element_num] + '_' + str(pH)+'_raw.dat'
            labelin1 = peptide_list[element_num] + '_' + str(pH)+ \
                '/' + peptide_list[element_num] + '_' + str(pH)+'_combined.dat'
            labelout1 = peptide_list[element_num] + '_' + str(pH)+ '/' + peptide_list[element_num] + '_' + str(pH)+'_min.dat'
            print(labelin1)
            print(labelout1)

        
        for index in range(len(start_index)):
                # print(index, start_index[index])
                if(pH==[]): 
                    # print(start_index[index])
                    label0 = peptide_list[element_num]+'/'+peptide_name + substring + '_'+energy_fxn+'_' + \
                        str(end_index[index])+'_'+str(start_index[index])+'_landscape.dat'
                    # print(label0)
                    if( index==0 and (not os.path.exists(label0))):
                        break
                    
                    df0 = pd.read_csv(label0)
                    # print(df0)
                    
                else:
                    label0 = peptide_list[element_num]+ '_' + str(pH) +'/'+peptide_name + substring + '_'+energy_fxn+'_' + \
                    str(end_index[index])+'_'+str(start_index[index])+'_landscape.dat'
                    df0 = pd.read_csv(label0)
                    
                df = pd.concat([df,df0])    
        
        # print(df)
        # sys.exit()
        # df = pd.concat([df0, df1])#, df2, df3, df4, df5])#, df6])#, df6])#, df7, df8, df9, df10])
        df.to_csv(labelin1, index=False)
        df_read = pd.read_csv(labelin1, delimiter=" ")
        X = unique(df_read['zcoord'])
        Y = unique(df_read['angle'])
        Z = unique(df_read['azimuthal'])
        
        #for ph-titration cases. 
        if(energy_fxn=='franklin2021'):
            df_read = df_read[df_read['loop']==0]
        
            
        #----extracting minimum scores-------------------
        mindata = []

        for i in range(len(X)):

#             newarr = df_read[df_read['zcoord'] == X[i]]
# #            print(newarr)
#             if(len(newarr) == 0):
#                 break

            for j in range(len(Y)):

                secondarr = df_read[(df_read['angle'] == Y[j]) & (df_read['zcoord'] == X[i])]
                if(len(secondarr) == 0):
                    break
                else:
                    arr = np.array(secondarr['total_score'])
#                    print(arr)
                    minpos = np.argmin(arr)
                    if(energy_fxn=='franklin2021'):
                        mindata.append(secondarr.iloc[minpos, :])
                    else:
                        mindata.append(secondarr.iloc[minpos, :])
        mindatadf = pd.DataFrame(mindata)
#        print(mindatadf)
        mindatadf.to_csv(labelout1, sep=" ", index=False)

        del df0#, df1#, df2
        # del df3, df4, df5#, df6
        # del df6, df7, df8
        # del df9, df10
        del df, df_read
        del X, Y

   # fclose(fileID);


def creating_1d_min_files(peptide_list, pH=[]):

    for element_num in range(len(peptide_list)):

        if(pH==[]):
            labelin = peptide_list[element_num] + '/' + \
                peptide_list[element_num] + '_min.dat'
            labelin1 = peptide_list[element_num] + '/' +\
                peptide_list[element_num] + '_combined.dat'
            
            labelout = peptide_list[element_num] + '/' + \
                peptide_list[element_num] + '_min1d.dat'
            
            labelout3 = peptide_list[element_num] + '/' + \
                peptide_list[element_num] + '_prev1d.dat'
            
            labelout1 = peptide_list[element_num] + '/' + \
                peptide_list[element_num] + '_prev.dat'
        else:
            labelin = peptide_list[element_num] + '_'+str(pH) + '/' + \
                peptide_list[element_num] + '_'+str(pH) + '_min.dat'
            labelout = peptide_list[element_num] + '_' + str(pH) + '/' + \
                peptide_list[element_num] + '_' + str(pH) + '_min1d.dat'
                
        print(labelin)
        print(labelout)

        labelout2 = peptide_list[element_num] + '/' + \
                peptide_list[element_num] + '_atd{}mindepth1d.dat'.format(0.0)
        df = pd.read_csv(labelin, delimiter=" ")
        df_zcoord = df[df['zcoord']==0]
        print(df_zcoord)
        df_zcoord.to_csv(labelout2, sep=" ", index=False)
        del df_zcoord
        
        X = unique(df['zcoord'])
        mindata = []
        mindata_prev=[]
        df0 = pd.read_csv(labelin1, delimiter=" ")
        df_prev = df0[(df0['azimuthal']==0)]
        # print(df_prev)
        df_prev.to_csv(labelout1, sep=" ", index=False)
        # df_prev_new = df_prev[df_prev['angle']==0]
        # df_prev_new.to_csv(labelout3, sep=" ", index=False)
        # del df_prev_new
        
        # arr_score = np.array(df['total_score'])
        # arr_depth = np.array(df['zcoord'])
        # minindex = np.where(arr_score == np.amin(df0['total_score']))
        # print(minindex[0])
        # zmin = arr_depth[minindex[0]] 
        # print(zmin[0])
        # df_angle1d = df[df['zcoord']==zmin[0].item()]
        
        # labelout2 = peptide_list[element_num] + '/' + \
        #         peptide_list[element_num] + '_atd{}mindepth1d.dat'.format(zmin)
        # print(df_angle1d)
        # df_angle1d.to_csv(labelout2, sep=" ", index=False)
        # #extract the angle dependence at min depth. 
        # del df_angle1d,arr_score,arr_depth,minindex,zmin
        
        
        for i in range(len(X)):
            newarr = df[df['zcoord'] == X[i]]
            arr = np.array(newarr['total_score'])

            # print(arr)
            minpos = np.argmin(arr)
            # print(minpos)
            mindata.append(newarr.iloc[minpos, :])
            
            newarr_prev = df_prev[df_prev['zcoord'] == X[i]]
            arr_prev = np.array(newarr_prev['total_score'])
            minpos_prev = np.argmin(arr_prev)
            mindata_prev.append(newarr_prev.iloc[minpos_prev, :])
            

        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout, sep=" ", index=False)
        
        mindatadf_prev = pd.DataFrame(mindata_prev)
        mindatadf_prev.to_csv(labelout3, sep=" ", index=False)

        

def combiningchargestatefiles(peptide_list, energy_fxn):

    for element_num in range(len(peptide_list)):

        start_index = ['0.000000', '10.000000', '20.000000','30.000000', '40.000000', '50.000000']
        end_index = ['10.000000', '20.000000', '30.000000','40.000000', '50.000000', '60.000000']  

        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA"):
            peptide_name = peptide_name2

        labelin1 = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_combined_chargestate.dat'
        label_in = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_raw.dat'

        print(labelin1)

        # print(str(end_index[2]))

        #fr = open(label1,'w')

        label0 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[0])+'_'+str(start_index[0])+'_chargestatevariation.dat'
        df0 = pd.read_csv(label0)
        label1 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[1])+'_'+str(start_index[1])+'_chargestatevariation.dat'
        df1 = pd.read_csv(label1)
        df1 = df1.fillna(0)
        label2 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[2])+'_'+str(start_index[2])+'_chargestatevariation.dat'
        df2 = pd.read_csv(label2)
        df2 = df2.fillna(0)
        label3 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[3])+'_'+str(start_index[3])+'_chargestatevariation.dat'
        df3 = pd.read_csv(label3)
        df3 = df3.fillna(0)
        label4 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[4])+'_'+str(start_index[4])+'_chargestatevariation.dat'
        df4 = pd.read_csv(label4)
        df4 = df4.fillna(0)
        label5 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
            end_index[5])+'_'+str(start_index[5])+'_chargestatevariation.dat'
        df5 = pd.read_csv(label5)
        df5 = df5.fillna(0)
       

        df = pd.concat([df0, df1, df2, df3, df4, df5])
        # df.fillna(0)
        df.to_csv(label_in, index=False)

        del df0
        del df1
        del df2
        del df3, df4, df5

        # del df,2,2,3,4,3]
#        df_read = pd.read_csv(label_in, delimiter=" ")
#        X = unique(df_read['zcoord'])
#        Y = unique(df_read['angle'])
#        Z = unique(df_read['azimuthal'])

        # print(X)
        # print(Y)
        # print(Z)

        #df = pd.DataFrame([],columns=list(df_read.columns))
#        final_arr = []
#        for i in range(len(X)):
        #newarr = df_read[df_read['zcoord']==X[i]]
        # if( len(newarr.index) == 0 ):
        #   break
#            for j in range(len(Y)):

        #   secondarr = newarr[newarr['angle'] == Y[j]]
        #   if( len(secondarr.index) == 0 ):
        #       break;
        #   for k in range(len(Z)):
#                 newarr = df_read[df_read['zcoord']==X[i]]
#                 secondarr = newarr[newarr['angle'] == Y[j]]
        # print(secondarr)
        # exit()
        # thirdarr = secondarr[secondarr['azimuthal']==0.000]#Z[k]]
        # print(thirdarr)
        # print(len(thirdarr.index))

#                 if( len(secondarr[secondarr['azimuthal']==0.000].index) == 0 ):
#                    break;
#                 else:
        #arr = np.array(thirdarr)
        # print(arr)
        # print(len(arr))
        # print(arr.shape[0])
        #avg_arr = np.mean(arr,axis=0)
        # print(avg_arr[0:4])
#                    final_arr.append(secondarr[secondarr['azimuthal']==0.0])
        # final_arr.append(arr.reshape(arr.shape[0],arr.shape[2])#avg_arr)
        # print(final_arr)
        # exit()

#        df=pd.concat(final_arr)
#        print(final_arr)
        # print(len(tfinal_arr)
        # print(len(final_arr[0][0]))
        # print(len(final_arr[0])
        # arr=np.reshape(final_arr,(len(final_arr),len(final_arr[0][0])))
        # print(final_arr.shape())
        # exit()
        #df = pd.DataFrame(arr, columns = list(df_read.columns))
        # print(df)
        # exit()

#        df.to_csv(labelin1, sep=" ", index=False, float_format='%.3f')
#        del final_arr,df, df_read


def main(args):

    all_sub_tests = ["ddG-of-insertion", "ddG-of-pH-insertion","amphiphilic-peptide-tilt-angle",
                     "tm-peptide-tilt-angle", "adsorbed-peptide-tilt-angle", "protein-tilt-angle","ph-titration"]

    # Read options from the command line
    parser = OptionParser(
        usage="usage %prog --energy_fxn franklin2019 --which_tests all")
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e', action="store",
                      help="Name of energy function weights file", )
    parser.add_option('--which_tests', '-w', action="store",
                      help="Specify tests run (comma separated list)", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options have been provided
    if (not Options.energy_fxn or not Options.which_tests):
        print("Missing required options --which_tests or --energy_fxn")
        sys.exit()

    # Check test categories
    test_names = []

    # Read path configuration file
    config = read_config.read_config()

    if (Options.which_tests == "all"):
        test_names = all_sub_tests
    else:
        test_names = Options.which_tests.split(",")
        # check that all names are valid
        for name in test_names:
            if name not in all_sub_tests:
                sys.exit("No such test " + name +
                         ". Exiting! or not among the tests 1-3 and 5-6")

            else:
                #fa_wb  f_elec_bilayer fa_imm_elec
                #1.179  -0.107  -0.015
                #1.375  0.286   -0.143
                #1.571  -0.179  -0.025
                #after fa_imm_modiefied to outside_memb fa_imm->0
                #1.629	    -1.561	    0.136
                #1.484	    0.571	    -0.462

                # wt_fa_water_to_bilayer = [1.0,0.863]#, 1.484] 
                # wt_fa_imm_elec = [0.01,0.001]#, 0.571]
                # wt_f_elec_bilayer = [0.128,0.152]#, -0.462]
                
                wt_fa_water_to_bilayer = [1.0]#, 1.484] 
                # wt_fa_imm_elec = [0.01]#, 0.571]
                # wt_f_elec_bilayer = [0.128]
                # wt_fa_water_to_bilayer = [0.863]#, 1.484] 
                # wt_fa_imm_elec = [0.001]#, 0.571]
                # wt_f_elec_bilayer = [0.152]#, -0.462]  
            
                # wt_fa_water_to_bilayer = [0.91]#, 1.484] 1.044,
                # wt_fa_imm_elec = [0.121]#, 0.571]0.107,
                # wt_f_elec_bilayer = [0.046]#, -0.462]  0.125,
        
                for i in range(0,len(wt_fa_water_to_bilayer)):
                    

                    if(Options.energy_fxn!="franklin2019"):
                        # datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                        #         Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                "franklin2021"+ "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"

                    else:
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        
                    # additional_folder='structure_fromA2ahelices_folder/'
                    additional_folder='test_f23_final/'
                    datadir = datadir + additional_folder
                    
                    # datadir = datadir + additional_folder + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                    #    wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                    # else:
                    # datadir = datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                    #          Options.energy_fxn + "/" + name + "/runs_for_paper1/"
                    #      # datadir = datadir + "fa_wb_" + str(0.966) + "_felecbilayer_" + str(
                    #     0.016) + "_fimm_" + str(0.379)   
                    print(datadir)
                    if (not os.path.isdir(datadir)):
                        sys.exit("No such test data available; test needs to finish before combining files")
                    os.chdir(datadir)

                    if(name == "tm-peptide-tilt-angle"):
                        # peptide_list = ['WALP23','1a11','1mp6','2nr1','polyA_cappedW', 'polyA_cappedY', 'AA20','AA25','AA28','AA30','AA35','AA40']#
                        # peptide_list = ['AA20','AA25','AA28','AA30','AA35','AA40']#,'WALP25','WALP31','WALP35','WALP39']#'1a11','1mp6','1pje',
                        # peptide_list = ['WALP31','WALP19','WALP25','WALP31','WALP35','WALP39']
                        peptide_list = ['WALP19','WALP25','WALP31','WALP35','WALP39']
                        #peptide_list = ['AA30','AA40']
                        # peptide_list = ['AA20','AA25','AA28','AA30','AA35','AA40']#, 'WALP19','WALP25','WALP31','WALP35','WALP39']
                        # peptide_list = ['WALP25','WALP35']#['WALP19','WALP25','WALP31','WALP35','WALP39']
                    elif(name == "adsorbed-peptide-tilt-angle"):
                        # peptide_list=['polylys_POPC','polylys_POPG']#,'polylys_DLPC']
                        peptide_list = ['1f0d','2mag']#,'1hu6','2mag']#,'LK_peptide_n6']
                        # peptide_list = ['LK_peptide_n6','2mag', '2mlt','1f0d','1f0d', '1f0g', '1hu5', '1hu6', '1hu7']#,'2khk','2mvm','5xng']
                    elif(name == "protein-tilt-angle"):
                        peptide_list = ['1fep', '1gzm', '1m0l', '1nqe', '1okc', '1p4t', '1qfg', '1qj8', '1qjp', '1r3j', '1yce', '2cfp', '2j8c', '2qom',
                                        '2x9k', '3aeh', '3dzm', '3pox', '3syb', '3wbn', '4afk', '4fqe', '4hyj', '4m48', '4n6h', '4rl8', '4uc2', '4x5n']
                    elif(name == "ddG-of-insertion"):
                        peptide_list = ['GL5','GL6', 'GL7', 'GL8', 'GWL6']#'GL5', 
                    elif(name == "ddG-of-pH-insertion"):
                        peptide_list = ['pHLIP-v1_4', 'pHLIP-v1_8', 'pHLIP-v2_4', 'pHLIP-v2_8', 'pHLIP-v3_4', 'pHLIP-v3_8', 'pHLIP-v4_4', 'pHLIP-v4_8', 'pHLIP-v5_4', 'pHLIP-v5_8', 'pHLIP-v6_4', 'pHLIP-v6_8', 'pHLIP-v7_4', 'pHLIP-v7_8', 'pHLIP-v8_4', 'pHLIP-v8_8',
                                        'pHLIP-v9_4', 'pHLIP-v9_8', 'pHLIP-v10_4', 'pHLIP-v10_8', 'pHLIP-v11_4', 'pHLIP-v11_8', 'pHLIP-v12_4', 'pHLIP-v12_8', 'pHLIP-v13_4', 'pHLIP-v13_8', 'pHLIP-v14_4', 'pHLIP-v14_8', 'pHLIP-v15_4', 'pHLIP-v15_8', 'pHLIP-v16_4', 'pHLIP-v16_8']
                    elif(name == "ph-titration"):
                        peptide_list = ['wt_W5','wt_W5_K12','wt_W5_K14','wt_Y5','wt_Y5_K12','wt_Y5_K14','polylys']#,'wt_W5_K12','wt_W5_K14_8','wt_Y5_8','wt_Y5_K12_8','wt_Y5_K14_8','polylys_8']
                    
                    elif(name == "amphiphilic-peptide-tilt-angle"):
                        peptide_list =['1h0a_h1','1q4g_h1','1q4g_h2','1q4g_h3','1q4g_h4','1rhz_h1','2hih_h1','2ziy_h1','3a7k_h1','3hyw_h1','3hyw_h2','3i9v_h1','3j5p_h1','3jw8_h1','3tij_h1','4hhr_h1','4hhr_h2','4hhr_h3','4m5e_h1','4nwz_h1','4qnd_h1',
                                        '4rp9_h3','4umw_h1','4ymk_h1','4ymk_h2','4ymk_h3','4zwn_h1','5ahv_h1','5dqq_h1','5ek8_h1','5f19_h3','5f19_h4','5lil_h1','5mlz_h2','5uz7_h1','5w7b_h1','5w7l_h1','5w7l_h2','5w7l_h3','6an7_h1','6d26_h1','6dvy_h1','6igk_h1']
                    

                    if(name == "amphiphilic-peptide-tilt-angle"):
                       combiningallfiles(peptide_list, Options.energy_fxn, tilt_axis=True, testname="amphiphilic-peptide-tilt-angle")
                    else:
                       combiningallfiles(peptide_list, Options.energy_fxn)
                    
                    creating_1d_min_files(peptide_list)
                    # combiningchargestatefiles(peptide_list, Options.energy_fxn)


if __name__ == "__main__":
    main(sys.argv)
