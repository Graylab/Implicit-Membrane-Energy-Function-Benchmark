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

import sys
import os
import read_config
import glob
import math
import fileinput
import pandas as pd
import numpy as np
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)


def combiningallfiles(datadir, peptide_list, energy_fxn):

    for element_num in range(len(peptide_list)):

        start_index = ['-60.000000','-50.000000','-40.000000', '-30.000000','-20.000000','-10.000000']#,'50.000000']#, '-30.000000',
                       #'-20.000000', '-10.000000', '0.000000', '10.000000', '20.000000']
        end_index = ['-50.000000','-40.000000','-30.000000', '-20.000000','-10.000000','0.000000']#,'60.000000']#, '-20.000000', '-10.000000',
                     #'0.000000', '10.000000', '20.000000', '30.000000']  # [-3, -2, -1, 0, 0, 0, 1, 2]

        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        pH = peptide_name2.split('_')[1]
        print("pH is: "+pH)

        if(peptide_name != "LK_peptide_n6"):
            peptide_name = peptide_name.split('_')[0]
        else:
            peptide_name = peptide_name2

        label_in = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_raw.dat'
        labelin1 = peptide_list[element_num] + \
            '/' + peptide_list[element_num] + '_combined.dat'
        labelout1 = peptide_list[element_num] + '/' + peptide_list[element_num] + '_min.dat'
        print(labelin1)
        print(labelout1)

        if(pH=='4'):
            
            label0 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                        str(end_index[0])+'_'+str(start_index[0])+'_landscape.dat'
            #label0 = peptide_list[element_num]+'/'+peptide_name+'_'+energy_fxn+'_' + \
                
            df0 = pd.read_csv(label0)
            label1 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                str(end_index[1])+'_'+str(start_index[1])+'_landscape.dat'
            df1 = pd.read_csv(label1)
            label2 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                str(end_index[2])+'_'+str(start_index[2])+'_landscape.dat'
            df2 = pd.read_csv(label2)
            label3 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                str(end_index[3])+'_'+str(start_index[3])+'_landscape.dat'
            df3 = pd.read_csv(label3)
            label4 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                str(end_index[4])+'_'+str(start_index[4])+'_landscape.dat'
            df4 = pd.read_csv(label4)
            label5 = peptide_list[element_num]+'/'+'ranked_0'+'_'+energy_fxn+'_' + \
                str(end_index[5])+'_'+str(start_index[5])+'_landscape.dat'
            df5 = pd.read_csv(label5)
        else:
            
            label0 = datadir + '/' + peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                        str(end_index[0])+'_'+str(start_index[0])+'_landscape.dat'
            #label0 = datadir + '/' + peptide_list[element_num]+'/'+peptide_name+'_'+energy_fxn+'_' + \
                
            df0 = pd.read_csv(label0)
            label1 = datadir + '/' + peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                str(end_index[1])+'_'+str(start_index[1])+'_landscape.dat'
            df1 = pd.read_csv(label1)

            label2 = peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                str(end_index[2])+'_'+str(start_index[2])+'_landscape.dat'
            df2 = pd.read_csv(label2)
            label3 = peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                str(end_index[3])+'_'+str(start_index[3])+'_landscape.dat'
            df3 = pd.read_csv(label3)
            label4 = peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                str(end_index[4])+'_'+str(start_index[4])+'_landscape.dat'
            df4 = pd.read_csv(label4)
            label5 = peptide_list[element_num]+'/'+peptide_name+'_disordered_'+energy_fxn+'_' + \
                str(end_index[5])+'_'+str(start_index[5])+'_landscape.dat'
            df5 = pd.read_csv(label5)

        df = pd.concat([df0, df1, df2, df3, df4, df5])#, df3, df4, df5])
        df.to_csv(labelin1, index=False)
        df_read = pd.read_csv(labelin1, delimiter=" ")
        X = unique(df_read['zcoord'])
        Y = unique(df_read['angle'])
        Z = unique(df_read['azimuthal'])

        final_arr = []
        for i in range(len(X)):
            newarr = df_read[df_read['zcoord'] == X[i]]
            if(len(newarr.index) == 0):
                break
            for j in range(len(Y)):

                secondarr = newarr[newarr['angle'] == Y[j]]
                if(len(secondarr.index) == 0):
                    break
                for k in range(len(Z)):
                    thirdarr = secondarr[secondarr['azimuthal'] == Z[k]]
                    if(len(thirdarr.index) == 0):
                        break
                    else:
                        arr = np.array(thirdarr)
                        avg_arr = np.mean(arr, axis=0)
                        # print(avg_arr[0:4])
                        if(energy_fxn=="franklin2021"):
                            final_arr.append(avg_arr[0:28])
                        else:
                            final_arr.append(avg_arr[0:25])

        # print(list(df_read.columns)[1:28])
        if(energy_fxn=="franklin2021"):
            df = pd.DataFrame(final_arr, columns=list(df_read.columns)[0:28])
        else:
            df = pd.DataFrame(final_arr, columns=list(df_read.columns)[0:25])
        
        df.to_csv(labelin1, sep=" ", index=False, float_format='%.3f')
        del df, df_read

        df_read = pd.read_csv(labelin1, delimiter=" ")

        mindata = []

        for i in range(len(X)):

            newarr = df_read[df_read['zcoord'] == X[i]]
            if(len(newarr) == 0):
                break

            for j in range(len(Y)):

                secondarr = newarr[newarr['angle'] == Y[j]]
                if(len(secondarr) == 0):
                    break
                else:
                    arr = np.array(secondarr['total_score'])
                    minpos = np.argmin(arr)
                    if(energy_fxn=="franklin2021"):
                        mindata.append(secondarr.iloc[minpos, 0:28])
                    else:
                        mindata.append(secondarr.iloc[minpos, 0:25])
        
        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout1, sep=" ", index=False, float_format='%.3f')

        del df0, df1, df2
        del df3, df4, df5
        #del df, df_read
        #del X, Y

   # fclose(fileID);


def creating_1d_min_files(peptide_list):

    for element_num in range(len(peptide_list)):

        labelin = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_min.dat'
        labelout = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_min1d.dat'

        print(labelin)
        print(labelout)

        df = pd.read_csv(labelin, delimiter=" ")
        X = unique(df['zcoord'])
        mindata = []

        for i in range(len(X)):
            newarr = df[df['zcoord'] == X[i]]
            arr = np.array(newarr['total_score'])

            # print(arr)
            minpos = np.argmin(arr)
            # print(minpos)
            mindata.append(newarr.iloc[minpos, 0:28])
            #mindata.append(newarr.iloc[minpos, 0:25])

        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout, sep=" ", index=False, float_format='%.3f')


def combiningchargestatefiles(datadir,peptide_list, energy_fxn):

    for element_num in range(len(peptide_list)):

        start_index = ['-60.000000','-50.000000','-40.000000', '-30.000000','-20.000000','-10.000000']#,'50.000000']#, '-30.000000',
                       #'-20.000000', '-10.000000', '0.000000', '10.000000', '20.000000']
        end_index = ['-50.000000','-40.000000','-30.000000', '-20.000000','-10.000000','0.000000']#,'60.000000']#, '-20.000000', '-10.000000',
           
        
        # start_index = ['0.000000','10.000000','20.000000', '30.000000','40.000000']#, '50.000000']#, '-30.000000',
        #               # '-20.000000', '-10.000000', '0.000000', '10.000000', '20.000000']

        # end_index = ['10.000000','20.000000','30.000000', '40.000000','50.000000']#,'60.000000']#, '-20.000000',
                     #'-10.000000', '0.000000', '10.000000', '20.000000', '30.000000']

        # start_index = ['-60.000000','-40.000000','-20.000000', '0.000000', '20.000000','40.000000']#[-1, -1, -1, -1, 0, 1, 1 , 1]
        # end_index = ['-40.000000','-20.000000', '0.000000', '20.000000', '40.000000','60.000000']#[-3, -2, -1, 0, 0, 0, 1, 2]

        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        pH = peptide_name2.split('_')[1]

        labelin1 = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_combined_chargestate.dat'
        label_in = peptide_list[element_num] + '/' + \
            peptide_list[element_num] + '_raw.dat'

        print(labelin1)

        # print(str(end_index[2]))

        #fr = open(label1,'w')
        if(pH=="8"):
            
            label0 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
            # label0 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
                end_index[0])+'_'+str(start_index[0])+'_chargestatevariation.dat'
            df0 = pd.read_csv(label0)
            #print(df0)
            df0 = df0.fillna(0)
            df0 = df0.replace("-nan", 0)
            print(df0)
            label1 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
                end_index[1])+'_'+str(start_index[1])+'_chargestatevariation.dat'
            df1 = pd.read_csv(label1)
            df1 = df1.fillna(0)
            
            label2 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
            end_index[2])+'_'+str(start_index[2])+'_chargestatevariation.dat'
            df2 = pd.read_csv(label2)
            df2 = df2.fillna(0)
            
            label3 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
                end_index[3])+'_'+str(start_index[3])+'_chargestatevariation.dat'
            df3 = pd.read_csv(label3)
            df3 = df3.fillna(0)
            
            label4 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
                end_index[4])+'_'+str(start_index[4])+'_chargestatevariation.dat'
            df4 = pd.read_csv(label4)
            df4 = df4.fillna(0)
            
            label5 = peptide_list[element_num] + '/' + 'ranked_0' + '_'+energy_fxn+'_'+str(
                end_index[5])+'_'+str(start_index[5])+'_chargestatevariation.dat'
            df5 = pd.read_csv(label5)
            df5 = df5.fillna(0)
        else:
            label0 = datadir + '/' + peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
            # label0 = datadir + '/' + peptide_list[element_num] + '/' + peptide_name + '_abinitio_relax_'+energy_fxn+'_'+str(
                end_index[0])+'_'+str(start_index[0])+'_chargestatevariation.dat'
            df0 = pd.read_csv(label0)
            # print(df0)
            df0 = df0.replace("-nan", 0)
            #print(df0)
            label1 = datadir + '/' + peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
                end_index[1])+'_'+str(start_index[1])+'_chargestatevariation.dat'
            df1 = pd.read_csv(label1)
            df1 = df1.fillna(0)
            label2 = peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
                end_index[2])+'_'+str(start_index[2])+'_chargestatevariation.dat'
            df2 = pd.read_csv(label2)
            df2 = df2.fillna(0)
            label3 = peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
                end_index[3])+'_'+str(start_index[3])+'_chargestatevariation.dat'
            df3 = pd.read_csv(label3)
            df3 = df3.fillna(0)
            label4 = peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
                end_index[4])+'_'+str(start_index[4])+'_chargestatevariation.dat'
            df4 = pd.read_csv(label4)
            df4 = df4.fillna(0)
            label5 = peptide_list[element_num] + '/' + peptide_name + '_disordered_'+energy_fxn+'_'+str(
                end_index[5])+'_'+str(start_index[5])+'_chargestatevariation.dat'
            df5 = pd.read_csv(label5)
            df5 = df5.fillna(0)
        #label6 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_'+str(
        #    end_index[6])+'_'+str(start_index[6])+'_chargestatevariation.dat'
        #df6 = pd.read_csv(label6)
        #df6 = df6.fillna(0)
        #label7 = peptide_list[element_num] + '/' + peptide_name + '_'+energy_fxn+'_' + str(
         #   end_index[7])+'_'+str(start_index[7])+'_chargestatevariation.dat'

        #df7 = pd.read_csv(label7)
        #df7 = df7.fillna(0)

        df = pd.concat([df0, df1, df2, df3, df4, df5])#, df3, df4, df5, df6, df7])
        # df.fillna(0)
        df.to_csv(label_in, index=False)

        del df0
        del df1
        del df2, df3, df4, df5
        #===========================      
       # # del df3, df4, df5
       # # del df6, df7

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

    all_sub_tests = ["ddG-of-insertion", "ddG-of-pH-insertion",
                     "tm-peptide-tilt-angle", "adsorbed-peptide-tilt-angle", "protein-tilt-angle"]

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
                wt_fa_water_to_bilayer = ["0.966"]#,"1.5"]#,"1.50"]  # ,"0.90","1.1","1.3","1.5"]
                wt_f_elec_bilayer = ["0.016"]#,"0.05","0.07","0.10"]
                wt_fa_imm_elec = ["0.379"]#,"0.10"]#["0.02","0.05","0.07","0.10"]
               
                for i in range(len(wt_fa_water_to_bilayer)):
                    for j in range(len(wt_f_elec_bilayer)):
                        for k in range(len(wt_fa_imm_elec)):

                            datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                    Options.energy_fxn + "/" + name + "/weights_from_test7/"
                            #datadir = config.benchmark_path + "data/" + \
                            #    Options.energy_fxn + "/" + name + "/"
                            datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                                wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                            print(datadir)
                            
                            #datadir1 = config.benchmark_path + "data/" + \
                            #    Options.energy_fxn + "/" + name + "/" + \
                            #    "weights_sweep_POPC_ordered_noPR/" + "fa_wb_" +\
                            #    str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(\
                            #    wt_f_elec_bilayer[j]) + "_fimm_" + str(wt_fa_imm_elec[k])

                            if (not os.path.isdir(datadir)):
                                sys.exit(
                                    "No such test data available; test needs to finish before combining files")
                            os.chdir(datadir)

                            if(name == "tm-peptide-tilt-angle"):
                                peptide_list = [
                                    '1a11', '1mp6', 'WALP23', '2nr1','polyA-cappedW', 'polyA-cappedY']
                            elif(name == "adsorbed-peptide-tilt-angle"):
                                peptide_list = [
                                    '1f0d', '1f0g', '1hu5', '1hu6', '1hu7', '2mag', 'LK_peptide_n6']
                            elif(name == "protein-tilt-angle"):
                                peptide_list = ['1fep', '1gzm', '1m0l', '1nqe', '1okc', '1p4t', '1qfg', '1qj8', '1qjp', '1r3j', '1yce', '2cfp', '2j8c', '2qom',
                                                '2x9k', '3aeh', '3dzm', '3pox', '3syb', '3wbn', '3wfd', '4afk', '4fqe', '4hyj', '4m48', '4n6h', '4rl8', '4uc2', '4x5n']
                            elif(name == "ddG-of-insertion"):
                                peptide_list = [
                                    'GL5', 'GL6', 'GL7', 'GL8', 'GWL6']
                            elif(name == "ddG-of-pH-insertion"):
    #                            peptide_list = ['pHLIP-v1_4', 'pHLIP-v1_8', 'pHLIP-v2_4', 'pHLIP-v2_8', 'pHLIP-v3_4', 'pHLIP-v3_8', 'pHLIP-v4_4', 'pHLIP-v4_8', 'pHLIP-v5_4', 'pHLIP-v5_8', 'pHLIP-v6_4', 'pHLIP-v6_8', 'pHLIP-v7_4', 'pHLIP-v7_8', 'pHLIP-v8_4', 'pHLIP-v8_8',
#                                                'pHLIP-v9_4', 'pHLIP-v9_8', 'pHLIP-v10_4', 'pHLIP-v10_8', 'pHLIP-v11_4', 'pHLIP-v11_8', 'pHLIP-v12_4', 'pHLIP-v12_8', 'pHLIP-v13_4', 'pHLIP-v13_8', 'pHLIP-v14_4', 'pHLIP-v14_8', 'pHLIP-v15_4', 'pHLIP-v15_8', 'pHLIP-v16_4', 'pHLIP-v16_8']
                                # peptide_list = ['pHLIP-v1_4','pHLIP-v2_4','pHLIP-v3_4','pHLIP-v4_4','pHLIP-v5_4','pHLIP-v6_4','pHLIP-v7_4','pHLIP-v8_4',\
                                #                 'pHLIP-v9_4','pHLIP-v10_4','pHLIP-v11_4','pHLIP-v12_4','pHLIP-v13_4','pHLIP-v14_4','pHLIP-v15_4','pHLIP-v16_4']

                                # peptide_list = ['pHLIP-v4_8','pHLIP-v15_8','pHLIP-v16_8']
                                peptide_list = ['pHLIP-v1_8','pHLIP-v2_8','pHLIP-v3_8','pHLIP-v4_8','pHLIP-v5_8','pHLIP-v6_8','pHLIP-v7_8','pHLIP-v8_8',\
                                                'pHLIP-v9_8','pHLIP-v10_8','pHLIP-v11_8','pHLIP-v12_8','pHLIP-v13_8','pHLIP-v14_8','pHLIP-v15_8','pHLIP-v16_8']
                               # end_index = ['10.000000','20.000000']
                               # start_index = ['0.000000','10.000000']
                                
                                
                              #  for l in range(len(peptide_list)):
                              #      for m in range(len(end_index)):
                              #          copy_cmdline ="scp" + " ../../weights_sweep_POPC/" + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                             #               wt_f_elec_bilayer[j]) + "_fimm_" + str(wt_fa_imm_elec[k]) + "/" + peptide_list[l] + "/" + (peptide_list[l]).split('_')[0]+"_franklin2021_"+end_index[m]+"_"+start_index[m]+"*" + " " + peptide_list[l]+"/."
                             #                                    copy_cmdline_sh = "scp" + " ../../weights_sweep_POPC/" + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                            #                wt_f_elec_bilayer[j]) + "_fimm_" + str(wt_fa_imm_elec[k]) + "/" + peptide_list[l] + "/" + (peptide_list[l]).split('_')[0]+"_"+(end_index[m]).split('.')[0]+"_"+(start_index[m]).split('.')[0]+"_*" + " " +peptide_list[l]+"/."
                                    
                                        #print(copy_cmdline)
                                        #print(copy_cmdline_sh)
                            #            os.system( copy_cmdline )
                            #            os.system( copy_cmdline_sh )

                            combiningallfiles(datadir,peptide_list,Options.energy_fxn)
                            # creating_1d_min_files(peptide_list)
                            # combiningchargestatefiles(datadir, peptide_list, Options.energy_fxn)


if __name__ == "__main__":
    main(sys.argv)
