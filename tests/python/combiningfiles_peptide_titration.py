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
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)


def combiningallfiles(peptide_list, energy_fxn):

    for element_num in range(len(peptide_list)):

        # start_index = ['0.000000', '20.000000','30.000000', '40.000000', '50.000000']
        # end_index = ['10.000000', '30.000000','40.000000', '50.000000', '60.000000']  # [-3, -2, -1, 0, 0, 0, 1, 2]
        for iter in range( 1,29 ):
            pH = (iter)*0.50
            
            start_index = ['0.000000','12.000000','55.000000']#, '0.000000', '10.000000', '20.000000','30.000000', '40.000000', '50.000000']
            end_index = ['12.000000','18.000000','60.000000']#, '10.000000', '20.000000', '30.000000','40.000000', '50.000000', '60.000000']  # [-3, -2, -1, 0, 0, 0, 1, 2]

            substring = "_after_relax"
            # substring = ""
            peptide_name2 = peptide_list[element_num]
            peptide_name = peptide_name2.split('_')[0]

            if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
                peptide_name = peptide_name2

            label_in = peptide_list[element_num] + '_' + str(pH) + '/' + \
                peptide_list[element_num] + '_raw.dat'
            labelin1 = peptide_list[element_num] + '_' + str(pH) + \
                '/' + peptide_list[element_num] + '_combined.dat'
            labelout1 = peptide_list[element_num] + '_' + str(pH) +  '/' + peptide_list[element_num] + '_min.dat'
            print(labelin1)
            print(labelout1)

            label0 = peptide_list[element_num]+ '_' + str(pH) +'/'+peptide_name + substring + '_'+energy_fxn+'_' + \
                str(end_index[0])+'_'+str(start_index[0])+'_landscape.dat'
            df0 = pd.read_csv(label0)
            label1 = peptide_list[element_num]+'_' + str(pH) + '/'+peptide_name + substring + '_'+energy_fxn+'_' + \
                str(end_index[1])+'_'+str(start_index[1])+'_landscape.dat'
            df1 = pd.read_csv(label1)
            label2 = peptide_list[element_num]+'_' + str(pH) + '/'+peptide_name + substring + '_'+energy_fxn+'_' + \
                str(end_index[2])+'_'+str(start_index[2])+'_landscape.dat'
            df2 = pd.read_csv(label2)
            
            

            df = pd.concat([df0, df1, df2])#, df3, df4, df5])#, df6])#, df6])#, df7, df8, df9, df10])
            # df = df0
            df.to_csv(labelin1, index=False)
            df_read = pd.read_csv(labelin1, delimiter=" ")
            X = unique(df_read['zcoord'])
            Y = unique(df_read['angle'])
            Z = unique(df_read['azimuthal'])
            
            #----extracting minimum scores-------------------
            mindata = []

            for i in range(len(X)):
                newarr1 = df_read[df_read['loop']==100]
                newarr = newarr1[newarr1['zcoord'] == X[i]]
                
    #            print(newarr)
                if(len(newarr) == 0):
                    break

                for j in range(len(Y)):

                    secondarr = newarr[newarr['angle'] == Y[j]]
                    if(len(secondarr) == 0):
                        break
                    else:
                        arr = np.array(secondarr['total_score'])
    #                    print(arr)
                        minpos = np.argmin(arr)
                        mindata.append(secondarr.iloc[minpos, 0:28])
            mindatadf = pd.DataFrame(mindata)
    #        print(mindatadf)
            mindatadf.to_csv(labelout1, sep=" ", index=False)

            del df0, df1, df2
            del df, df_read
            del X, Y

   # fclose(fileID);


def creating_1d_min_files(peptide_list):

    for element_num in range(len(peptide_list)):

        labelin = peptide_list[element_num] + '_' + str(pH) + '/' + \
            peptide_list[element_num] + '_min.dat'
        labelout = peptide_list[element_num] + '_' + str(pH) + '/' + \
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

        mindatadf = pd.DataFrame(mindata)
        mindatadf.to_csv(labelout, sep=" ", index=False)

def map_min_energy_pH(peptide_list):
    
    for element_num in range(len(peptide_list)):
        df_new = []
        for iter in range( 1,28 ):
            pH = (iter)*0.50
            
            peptide_name2 = peptide_list[element_num]
            peptide_name = peptide_name2.split('_')[0]
            residue_list=['Glu','His','Tyr','Lys','Asp'] 
            
            if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
                peptide_name = peptide_name2

            label_in = peptide_list[element_num] + '_' + str(pH) + '/' + \
                peptide_list[element_num] + '_min.dat'
            df = pd.read_csv(label_in, delimiter=" ")    
            
            df.loc[:,'pH'] = pH
            df_new.append(df)
            
        final_df = pd.concat(df_new)
        final_df = final_df.sort_values(['zcoord', 'pH'],ascending = [True, True])
        
        df_depth = unique(final_df['zcoord'])
        df_angle = unique(final_df['azimuthal'])
        df_ph = unique(final_df['pH'])

        depth_resize = np.resize(final_df['zcoord'],(len(df_depth), len(df_ph)))
        theta_resize = np.resize(final_df['azimuthal'],(len(df_depth), len(df_ph)))
        ph_resize = np.resize(final_df['pH'],(len(df_depth), len(df_ph)))
        
        outfile_name = peptide_name+'/'+peptide_name+'-'+'phdistribution'+'_heatmap.png'
        print(outfile_name)
        plot_heatmap(depth_resize,ph_resize,theta_resize, outfile_name, 'Depth (Å)', 'pH')
                          

def combiningchargestatefiles(peptide_list, energy_fxn, lastloop=100):

    for element_num in range(len(peptide_list)):

        start_index = ['0.000000','12.000000','55.000000']
        end_index = ['12.000000','18.000000','60.000000']  

        substring = "_after_relax"
        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
            peptide_name = peptide_name2

        label_out = peptide_list[element_num] + '_combined_titration.dat'
        label_in = peptide_list[element_num] + '_raw.dat' 

        # print(str(end_index[2]))
        df = []
        #fr = open(label1,'w')
        for iter in range( 1,28 ):
            pH = (iter)*0.50
            label0 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[0])+'_'+str(start_index[0])+'_chargestatevariation.dat'
            df0 = pd.read_csv(label0)
            label1 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[1])+'_'+str(start_index[1])+'_chargestatevariation.dat'
            df1 = pd.read_csv(label1)
            label2 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[2])+'_'+str(start_index[2])+'_chargestatevariation.dat'
            df2 = pd.read_csv(label2)
            
            df_full = pd.concat([df0, df1, df2])
            
            df_full.to_csv(label_in, index=False)
            
            df_read = pd.read_csv(label_in, delimiter=" ")
            df_new = df_read[df_read['loop']==lastloop]
            
            df_new.loc[:,'pH'] = pH
            df.append(df_new)#df3, df4, df5])
            # print(df)
            del df0, df_read, df_new, df1, df2, df_full
            # sys.exit()
            # df.fillna(0)
            print('pH: '+str(pH))
        final_df = pd.concat(df)
        final_df.to_csv(label_out,sep=" ",index=False, na_rep='NULL')

        # del df1
        # del df2
        # del df3, df4, df5

        # del df,2,2,3,4,3]
        
        

        # print(X)
        # print(Y)
        # print(Z)


#        df.to_csv(labelin1, sep=" ", index=False, float_format='%.3f')
#        del final_arr,df, df_read

def combiningchargestateresfiles(peptide_list, energy_fxn, lastloop=100):

    for element_num in range(len(peptide_list)):

        start_index = ['0.000000','12.000000','55.000000']
        end_index = ['12.000000','18.000000','60.000000']  

        substring = "_after_relax"
        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
            peptide_name = peptide_name2

        label_out = peptide_list[element_num] + '_combined_chargestateres.dat'
        label_in = peptide_list[element_num] + '_chargestateres_raw.dat'

        

        # print(str(end_index[2]))
        df = []
        #fr = open(label1,'w')
        for iter in range( 1,28 ):
            pH = (iter)*0.50
            label0 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[0])+'_'+str(start_index[0])+'_chargestateresidue.dat'
            df0 = pd.read_csv(label0)
            label1 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[1])+'_'+str(start_index[1])+'_chargestateresidue.dat'
            df1 = pd.read_csv(label1)
            label2 = peptide_list[element_num] + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
                end_index[2])+'_'+str(start_index[2])+'_chargestateresidue.dat'
            df2 = pd.read_csv(label2)
            
            df_full = pd.concat([df0, df1, df2])
            
            df_full.to_csv(label_in, index=False)
            
            df_read = pd.read_csv(label_in, delimiter=" ")
            # df_new = df_read[df_read['loop']==lastloop]
            df_new = df_read
            df_new.loc[:,'pH'] = pH
            df.append(df_new)#df3, df4, df5])
            # print(df)
            del df0, df_read, df_new, df1, df2, df_full
            # sys.exit()
            # df.fillna(0)
            print('pH: '+str(pH))
        final_df = pd.concat(df)
        final_df.to_csv(label_out,sep=" ",index=False, na_rep='NULL')
    # return(final_df)
        # del df1
        # del df2
        # del df3, df4, df5

        # del df,2,2,3,4,3]
        
        

        # print(X)
        # print(Y)
        # print(Z)


#        df.to_csv(labelin1, sep=" ", index=False, float_format='%.3f')
#        del final_arr,df, df_read

def combiningchargestateresfiles_pdb(peptide_list_name, energy_fxn):
    
    start_index = ['0.000000','12.000000','55.000000']
    end_index = ['12.000000','18.000000','60.000000']  

    substring = "_after_relax"
    peptide_name2 = peptide_list_name
    peptide_name = peptide_name2.split('_')[0]
    if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
        peptide_name = peptide_name2

    label_out = peptide_list_name + '_combined_chargestateres.dat'
    label_in = peptide_list_name + '_chargestateres_raw.dat'   

    # print(str(end_index[2]))
    df = []
    #fr = open(label1,'w')
    for iter in range( 1,28 ):
        pH = (iter)*0.50
        label0 = peptide_list_name + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
            end_index[0])+'_'+str(start_index[0])+'_chargestateresidue.dat'
        df0 = pd.read_csv(label0)
        label1 = peptide_list_name + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
            end_index[1])+'_'+str(start_index[1])+'_chargestateresidue.dat'
        df1 = pd.read_csv(label1)
        label2 = peptide_list_name + '_' + str(pH) + '/' + peptide_name + substring + '_'+energy_fxn+'_'+str(
            end_index[2])+'_'+str(start_index[2])+'_chargestateresidue.dat'
        df2 = pd.read_csv(label2)
        
        df_full = pd.concat([df0, df1, df2])
        
        df_full.to_csv(label_in, index=False)
        
        df_read = pd.read_csv(label_in, delimiter=" ")
        # df_new = df_read[df_read['loop']==lastloop]
        df_new = df_read
        df_new.loc[:,'pH'] = pH
        df.append(df_new)#df3, df4, df5])
        # print(df)
        del df0, df_read, df_new, df1, df2, df_full
        # sys.exit()
        # df.fillna(0)
        print('pH: '+str(pH))
    final_df = pd.concat(df)
    return(final_df)


def extractingchargestateresfiles(peptide_list, energy_fxn, lastloop=100, dist=0):
     
    for element_num in range(len(peptide_list)):
        substring = "_after_relax"
        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
            peptide_name = peptide_name2
            
        df_new = combiningchargestateresfiles_pdb(peptide_list[element_num], energy_fxn)
        df_read = df_new[(df_new['azimuthal']<100) & (df_new['depth']==dist)]
        print(df_read)
        # sys.exit()
        # df_read = df_new[(df_new['azimuthal']<100)]
        
        df_depth = unique(df_read['depth'])
        df_angle = unique(df_read['azimuthal'])
        df_ph = unique(df_read['pH'])
        df_loop = unique(df_read['loop'])
        df_sequence_length = unique(df_read['residue#'])
        max_seq_length = max(df_sequence_length)
        
        label_out = peptide_list[element_num] + '/' + peptide_list[element_num] + '_combined_chargestate_res_at'+str(dist)+'.dat'
        # df=[]
        df=pd.DataFrame([])
        for depth in range(len(df_depth)):
            
            # if((df_depth[depth]<14) or (df_depth[depth]>18)):
            #     if(df_depth[depth]<55):
            #         print("depth: {}".format(df_depth[depth]))
            #         continue
            # if(df_depth[depth]!=14):
            #     continue
            for angle in range(len(df_angle)):
                print("angle: {}".format(df_angle[angle]))
                if(df_angle[angle]>100):
                    continue
                for ph in range(len(df_ph)):
                    print("ph: {}".format(ph))
                    # df=[]
                    fraction_dict = {
                            "depth": [df_depth[depth]],
                            "tilt": [0.0],
                            "azimuthal": [df_angle[angle]],
                            "ph":[df_ph[ph]]
                            }
                    for res_num in range(len(df_sequence_length)):
                        res, fraction = calculate_protonation(df_read,df_sequence_length[res_num],\
                            max_seq_length,df_depth[depth],df_angle[angle], df_ph[ph], lastloop)
                        if(res==[]):
                            continue
                        
                        key_dict = '{}_{}'.format(res,res_num)
                        fraction_dict[key_dict] = fraction
                        
                    df_row=pd.DataFrame.from_dict(fraction_dict)
                    df = pd.concat([df,df_row], axis=0)
                    #df = pd.co.append(df_row)
                    # print(df)
                    # if(ph>2.0):
                    #     print(df)
                    #     sys.exit()
                    # del fraction_dict, df_row
                
            # print(df)
            final_data = pd.DataFrame(df)
            # print(final_data)
            final_data.to_csv(label_out,sep=" ",index=False, na_rep='NULL')
        # sys.exit()

def extract_pka_perres(peptide_list, dist=0):
    
    from scipy.optimize import curve_fit
        
    for element_num in range(len(peptide_list)):

        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
            peptide_name = peptide_name2

        labelout1 = peptide_name + '/' + peptide_name + '_combined_chargestate_res_at'+str(dist)+'.dat'
        if(not os.path.exists(labelout1)):
            sys.exit('generate the file {}'.format(labelout1))
            
        df_read = pd.read_csv(labelout1, delimiter=" ")
        
        depth = unique(df_read['depth'])
        theta = unique(df_read['azimuthal']) 
        print(theta)
        df=pd.DataFrame([])
        for element in df_read.columns:
            if(element not in ['depth','tilt','azimuthal','ph']):
                print(element)
                newarr = []
                
                
                # df = np.empty((4, 0))   
                for i in range(len(depth)):
                    for j in range(len(theta)):
                        newarr = df_read[(df_read['depth']==depth[i]) & (df_read['azimuthal']==theta[j])]
                        residue_element = element.split('_')[0]
                        # print(residue_element)
                        if(residue_element == 'TYR' or residue_element == 'LYS'):
                            # print('in loop 1')
                            popt_titration, pcov_titration = curve_fit(function_tyr_curve, newarr['ph'], 1.0 - newarr[element] )
                            # print(popt_titration[0])
                        elif(residue_element == 'ASP' or residue_element == 'GLU' or residue_element == 'HIS'):
                            # print('in loop 2')
                            popt_titration, pcov_titration = curve_fit(function_asp_curve, newarr['ph'], 1.0 - newarr[element])
                            # print(popt_titration[0])
                        # print(df)
                        pka_dict = {
                            "depth": [depth[i]],
                            "azimuthal": [theta[j]],
                            "index":[int(element.split('_')[1])],
                            "res":[element],
                            "pka":[popt_titration[0].round(3)]
                            }
                        
                        df_row = pd.DataFrame.from_dict(pka_dict)#, columns=['depth','azimuthal','element','pka'])    
                        # print(df_row)
                        
                        df = pd.concat([df,df_row], axis=0)
                        # df  = np.column_stack((df, np.array((depth[i], theta[j], element, popt_titration[0].round(3)))))

        print(df)
        df = df.sort_values(by=['azimuthal','index'])
        print(df)
        # depth_resize = np.resize(df[0],(len(depth), len(theta)))
        # theta_resize = np.resize(df[1],(len(depth), len(theta)))
        # pka_resize = np.resize(df[2],(len(depth), len(theta)))
        # outfile_name = peptide_name+'/'+peptide_name+'-'+element+'_heatmap.png'
        # plot_heatmap(depth_resize,theta_resize,pka_resize, outfile_name,'Depth (Å)','Rotation angle (degree)')
        label_out = peptide_name+'/'+peptide_name +'_finaldata_perres_at'+str(dist)+'.dat'
        df.to_csv(label_out,sep=" ",index=False, na_rep='NULL')
        # np.savetxt(peptide_name+'/'+peptide_name+'-'+ element +'_finaldata_perres_at15.dat', np.transpose(df), delimiter=' ')#, fmt='%1.3f')

def calculate_protonation(df,residue_number,max_seq_length,depth,angle, ph, num_records):
        # newarr = df[df['depth']==depth]
        # newarr1 = newarr[newarr['azimuthal'] == angle]
        # newarr2 = newarr2[newarr2['pH'] == ph]
        newarr = df[(df['depth']==depth) & (df['azimuthal'] == angle) & (df['pH'] == ph) & (df['residue#'] == residue_number)]                
        if(residue_number==1 or residue_number==max_seq_length):
            newarr[['res_name','end_name']] = newarr.residue.str.split(':',expand=True)
        else:
            newarr.loc[:,'res_name'] = newarr.loc[:,'residue']
            newarr.loc[:,'end_name'] = newarr.loc[:,'residue']
            
        # print('res:{} of {}'.format(residue_number,max_seq_length))
        #following the definition for membraneenergylandscapesampler
        # print(newarr['residue'].str.split(':')[0][0])
        if(newarr.iloc[0]['res_name'] == "LYS" or newarr.iloc[0]['res_name'] == "LYS_D"):
            # print('res:{} , #:{}'.format(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('LYS_D')).sum()/num_records))
            return(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('LYS_D')).sum()/num_records)
            # df['col'].eq('exact_string')).any()
        elif(newarr.iloc[0]['res_name'] == "TYR" or newarr.iloc[0]['res_name'] == "TYR_D"):
            # print('res:{} , #:{}'.format(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('TYR_D')).sum()/num_records))
            return(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('TYR_D')).sum()/num_records)
        elif(newarr.iloc[0]['res_name'] == "ASP" or newarr.iloc[0]['res_name'] == "ASP_P1"\
            or newarr.iloc[0]['res_name'] == "ASP_P2"):
            return(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('ASP_P1').sum() + newarr['res_name'].eq('ASP_P2').sum())/num_records)
        elif(newarr.iloc[0]['res_name'] == "GLU" or newarr.iloc[0]['res_name'] == "GLU_P1"\
            or newarr.iloc[0]['res_name'] == "GLU_P2"):
            return(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('GLU_P1').sum() + newarr['res_name'].eq('GLU_P2').sum())/num_records)
        elif(newarr.iloc[0]['res_name'] == "HIS" or newarr.iloc[0]['res_name'] == "HIS_P"\
            or newarr.iloc[0]['res_name'] == "HIS_D"):
            return(newarr.iloc[0]['res_name'][0:3], (newarr['res_name'].eq('HIS_P')).sum()/num_records)
        else:
            return([],0.0)
def function_tyr_curve(x, pKa):
    return 1.0/(1.0 + pow(10.0,(x-pKa)))


def function_asp_curve(x,pKa):
    return 1.0/(1.0 + pow(10.0,(-x+pKa)))

def plot_avg_protonation(peptide_name):
    import matplotlib.pyplot as plt
    residue_list=['Glu','His','Tyr','Lys','Asp']
    # residue_list=['Glu','His','Tyr','Lys','Asp']
    
    labelout1 = peptide_name + '_combined_titration.dat'
    df_read = pd.read_csv(labelout1, delimiter=" ")
    outfile = peptide_name + '/' + "plot_depth_dependent_protonation.png"
    
    outfile2 = peptide_name + '/' + "plot_theta_dependent_protonation.png"
    
    color_list = ['#fafa6e','#c4ec74','#92dc7e','#64c987','#39b48e','#089f8f','#00898a','#08737f','#215d6e','#2a4858']
    color_list_20 = ['#fafa6e','#e0f470','#c7ed73','#aee678','#97dd7d','#81d581','#6bcc86','#56c28a','#42b98d','#2eaf8f','#18a48f','#009a8f','#00908d','#008589','#007b84','#0c707d','#196676','#215c6d','#275263','#2a4858']
    color_list_30 = ['#fafa6e','#e9f66f','#d8f271','#c7ed73','#b7e876','#a8e379','#99de7c','#8ad87f','#7bd383','#6dcd85','#60c788','#52c08a','#45ba8c','#37b48e','#2aad8f','#1ca68f','#0a9f8f','#00998e','#00928d','#008b8b','#008489','#007d85','#027681','#0e6f7d','#166978','#1d6272','#225b6c','#255566','#284e5f','#2a4858']
    color_list_40 = ['#fafa6e','#edf76f','#e0f470','#d4f171','#c8ed73','#bcea75','#b0e678','#a5e27a','#99de7c','#8eda7f','#83d681','#79d283','#6ecd85','#64c987','#5ac489','#50bf8b','#46bb8c','#3cb68d','#32b18e','#28ac8f','#1ea78f','#12a28f','#039d8f','#00988e','#00938d','#008e8c','#00898a','#008488','#007e86','#007983','#057480','#0e6f7d','#156a79','#1a6575','#1e6071','#225b6c','#255667','#275163','#294d5d','#2a4858']
    
    depth_pH = unique(df_read['zcoord'])
    
    theta_pH = unique(df_read['azimuthal'])
    
    nrows = 3
    ncols = 2
    width = 7.5*ncols
    height = 6*nrows
    
    theme = {'axes.grid': True,
            'grid.linestyle': '',
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            "font.weight": 'regular',
            'xtick.color': 'black',
            'ytick.color': 'black',
            "axes.titlesize": 20,
            "axes.labelsize": 20
        }
    
    plt.rcParams['figure.figsize'] = width, height
    plt.rcParams.update(theme)
    plt.xlim([0,14])
    plt.ylim([0,1])
    plt.xlabel('pH')
    plt.ylabel('Dissociation')
    for element in range(len(residue_list)):
        col_name = 'avg_'+residue_list[element].upper()
        
        
        col_num = element//3 + 1
        row_num = element%3 + 1
        print('row {}, col {}'.format(row_num, col_num))
        
        plt.subplot(nrows, ncols, element+1)
        plt.title('residue: '+ residue_list[element])
        

        for i in range(len(depth_pH)):
            if(depth_pH[i]%2 == 0):
                newarr=df_read[df_read['zcoord']==depth_pH[i]]
                newarr2 = newarr[newarr['azimuthal'] == 0]
                plt.plot(newarr2['pH'],1.0 - newarr2[col_name],marker='D',color=color_list_40[i])
                 
    plt.tight_layout()
    plt.savefig( outfile, dpi=600, transparent=True  )
    plt.close()

def plot_avg_protonation_theta(peptide_name):
    import matplotlib.pyplot as plt
    residue_list=['Glu','His','Tyr','Lys','Asp']
    # residue_list=['Glu','His','Tyr','Lys','Asp']
    
    labelout1 = peptide_name + '_combined_titration.dat'
    df_read = pd.read_csv(labelout1, delimiter=" ")
    
    outfile2 = peptide_name + '/' + "plot_theta_dependent_protonation.png"
    
    color_list = ['#fafa6e','#c4ec74','#92dc7e','#64c987','#39b48e','#089f8f','#00898a','#08737f','#215d6e','#2a4858']
    color_list_20 = ['#fafa6e','#e0f470','#c7ed73','#aee678','#97dd7d','#81d581','#6bcc86','#56c28a','#42b98d','#2eaf8f','#18a48f','#009a8f','#00908d','#008589','#007b84','#0c707d','#196676','#215c6d','#275263','#2a4858']
    color_list_30 = ['#fafa6e','#e9f66f','#d8f271','#c7ed73','#b7e876','#a8e379','#99de7c','#8ad87f','#7bd383','#6dcd85','#60c788','#52c08a','#45ba8c','#37b48e','#2aad8f','#1ca68f','#0a9f8f','#00998e','#00928d','#008b8b','#008489','#007d85','#027681','#0e6f7d','#166978','#1d6272','#225b6c','#255566','#284e5f','#2a4858']
    color_list_40 = ['#fafa6e','#edf76f','#e0f470','#d4f171','#c8ed73','#bcea75','#b0e678','#a5e27a','#99de7c','#8eda7f','#83d681','#79d283','#6ecd85','#64c987','#5ac489','#50bf8b','#46bb8c','#3cb68d','#32b18e','#28ac8f','#1ea78f','#12a28f','#039d8f','#00988e','#00938d','#008e8c','#00898a','#008488','#007e86','#007983','#057480','#0e6f7d','#156a79','#1a6575','#1e6071','#225b6c','#255667','#275163','#294d5d','#2a4858']
    
    
    depth_pH = unique(df_read['zcoord'])
    
    theta_pH = unique(df_read['azimuthal'])
    
    nrows = 3
    ncols = 2
    width = 7.5*ncols
    height = 6*nrows
    
    theme = {'axes.grid': True,
            'grid.linestyle': '',
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            "font.weight": 'regular',
            'xtick.color': 'black',
            'ytick.color': 'black',
            "axes.titlesize": 20,
            "axes.labelsize": 20
        }
    #this part definitely needs to be better
    plt.rcParams['figure.figsize'] = width, height
    plt.rcParams.update(theme)
    
    for element in range(len(residue_list)):
        col_name = 'avg_'+residue_list[element].upper()
        
        row_num = element//3 + 1
        col_num = element%3 + 1
        print('row {}, col {}'.format(nrows, ncols))
        

        plt.subplot(nrows, ncols, element+1)
        plt.title('residue: '+ residue_list[element])
        plt.xlabel('pH')
        plt.ylabel('Dissociation')
        
        for i in range(len(theta_pH)):
            if(theta_pH[i]%3 == 0):
                newarr=df_read[df_read['azimuthal']==theta_pH[i]]
                newarr2 = newarr[newarr['zcoord'] == 15]
                plt.plot(newarr2['pH'],1.0 - newarr2[col_name],marker='D',color=color_list_40[i])
    
    plt.xlim([0,14])
    plt.ylim([0,1])             
    plt.tight_layout()
    plt.savefig( outfile2, dpi=600, transparent=True  )
    plt.close()
        
def plot_heatmap(x_axis,y_axis,z_map, outfile_name, xlabel, ylabel):
    import matplotlib.pyplot as plt
    
    theme = {'axes.grid': True,
            'grid.linestyle': '',
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            "font.weight": 'regular',
            'xtick.color': 'black',
            'ytick.color': 'black',
            "axes.titlesize": 20,
            "axes.labelsize": 20
        }
    plt.rcParams['figure.figsize'] = 10,10
    plt.rcParams.update(theme)
    fig, ax = plt.subplots()
    # ax.pcolormesh
    c = plt.contourf(x_axis, y_axis, z_map, cmap='viridis', vmin=z_map.min(), vmax=z_map.max())
    # set the limits of the plot to the limits of the data
    ax.axis([x_axis.min(), x_axis.max(), y_axis.min(), y_axis.max()])
    fig.colorbar(c, ax=ax)
    plt.ylabel(ylabel)#'Rotation angle (degree)')
    plt.xlabel(xlabel)#'Depth (Å)')
    plt.tight_layout()
    plt.savefig( outfile_name, dpi=600, transparent=True  )
    plt.close()
    

def extract_pka(peptide_list):
    from scipy.optimize import curve_fit
    
    for element_num in range(len(peptide_list)):

        peptide_name2 = peptide_list[element_num]
        peptide_name = peptide_name2.split('_')[0]
        if(peptide_name == "LK" or peptide_name=="polyA" or peptide_name=="wt"):
            peptide_name = peptide_name2

        
        residue_list=['Glu','His','Tyr','Lys','Asp']
        
        plot_avg_protonation(peptide_name)
        plot_avg_protonation_theta(peptide_name)

        labelout1 = peptide_name + '_combined_titration.dat'
        df_read = pd.read_csv(labelout1, delimiter=" ")
        
        depth_pH = unique(df_read['zcoord'])
        theta_pH = unique(df_read['azimuthal']) 
        
        for element in range(len(residue_list)):
            col_name = 'avg_'+residue_list[element].upper() 
            
            # print(df_read[col_name][0])
            # print(df_read[col_name][0]=='NaN')
            # print(df_read[col_name][0]=='nan')
            # print(math.isnan(df_read[col_name][0]))
            
            if(math.isnan(df_read[col_name][0])==True):
                continue
            # sys.exit()
            print(col_name)
            newarr = []
            newarr2 = []
            df = np.empty((3, 0))
            for i in range(len(depth_pH)):
                for j in range(len(theta_pH)):
                    newarr=df_read[df_read['azimuthal']==theta_pH[j]]
                    newarr2 = newarr[newarr['zcoord'] == depth_pH[i]]
                    
                    # print(newarr2[col_name][0])
                    # print(newarr2[col_name][0]=='NaN')
                    
                    
                    if(residue_list[element] == 'Tyr' or residue_list[element] == 'Lys'):
                        
                        popt_titration, pcov_titration = curve_fit(function_tyr_curve, newarr2['pH'], 1.0 - newarr2[col_name] )
                        # print(popt_titration[0])
                    elif(residue_list[element] == 'Asp' or residue_list[element] == 'Glu'):
                        
                        popt_titration, pcov_titration = curve_fit(function_asp_curve, newarr2['pH'], 1.0 - newarr2[col_name])
                        # print(popt_titration[0])
                        #print(pcov_titration)

                        # if(peptide_name == "veeks"):
                        #     print(newarr2['pH'])
                        #     print(newarr2[col_name])
                        #     print(popt_titration[0])
                            
                    df  = np.column_stack((df, np.array((depth_pH[i], theta_pH[j], popt_titration[0].round(3)))))
            # print(df)
            depth_resize = np.resize(df[0],(len(depth_pH), len(theta_pH)))
            theta_resize = np.resize(df[1],(len(depth_pH), len(theta_pH)))
            pka_resize = np.resize(df[2],(len(depth_pH), len(theta_pH)))
            outfile_name = peptide_name+'/'+peptide_name+'-'+residue_list[element]+'_heatmap.png'
            plot_heatmap(depth_resize,theta_resize,pka_resize, outfile_name,'Depth (Å)','Rotation angle (degree)')
            
    
            np.savetxt(peptide_name+'/'+peptide_name+'-'+residue_list[element]+'_finaldata.dat', np.transpose(df), delimiter=' ', fmt='%1.3f')
        del df_read
def main(args):

    all_sub_tests = ["ddG-of-insertion", "ddG-of-pH-insertion",
                     "tm-peptide-tilt-angle", "adsorbed-peptide-tilt-angle", "protein-tilt-angle", "ph-titration"]

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

                # wt_fa_water_to_bilayer = ["1.629"]  # ,"0.90","1.1","1.3","1.5"]
                # wt_fa_imm_elec = ["-1.561"]  # ,"0.02","0.05","0.07","0.10"]
                # wt_f_elec_bilayer = ["0.136"]  # ,"0.02","0.05","0.07","0.10"]
                wt_fa_water_to_bilayer = ["0.966"]#,"1.168" ]#, 1.484] 
                wt_fa_imm_elec = ["0.379"]#,"0.045"]#, 0.571]
                wt_f_elec_bilayer = ["0.016"]#,"0.126"]#, -0.462]  
        
        
                for i in range(len(wt_fa_water_to_bilayer)):
                    for j in range(len(wt_f_elec_bilayer)):
                        for k in range(len(wt_fa_imm_elec)):

                            datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                Options.energy_fxn + "/" + name + "/weights_from_test7_titration/"
                            # datadir = config.benchmark_path + "data/" + \
                            # Options.energy_fxn + "/" + name + "/weights_from_test7/"
                            datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                                wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                            
                            print(datadir)
                            if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                            os.chdir(datadir)

                            if(name == "tm-peptide-tilt-angle"):
                                # peptide_list = ['1a11', 'WALP23','polyA_cappedW', 'polyA_cappedY', 'AA20','AA25','AA28','AA30','AA35','AA40', 'WALP19','WALP25','WALP31','WALP35','WALP39']
                                #peptide_list = ['AA30','AA40']
                                # peptide_list = ['AA20','AA25','AA28','AA30','AA35','AA40']#, 'WALP19','WALP25','WALP31','WALP35','WALP39']
                                # peptide_list = ['WALP25','WALP35']#['WALP19','WALP25','WALP31','WALP35','WALP39']
                                peptide_list = ['pHLIP-L16H','pHLIP-wt']
                            elif(name == "adsorbed-peptide-tilt-angle"):
                                # peptide_list = [
                                #     '1f0d', '1f0g', '1hu5', '1hu6', '1hu7', '2mag', 'LK_peptide_n6']#,'2khk','2mvm','5xng']
                                peptide_list = ['AAsp','AHis','AGlu','ALys','ATyr']
                            elif(name == "protein-tilt-angle"):
                                peptide_list = ['1fep', '1gzm', '1m0l', '1nqe', '1okc', '1p4t', '1qfg', '1qj8', '1qjp', '1r3j', '1yce', '2cfp', '2j8c', '2qom',
                                                '2x9k', '3aeh', '3dzm', '3pox', '3syb', '3wbn', '4afk', '4fqe', '4hyj', '4m48', '4n6h', '4rl8', '4uc2', '4x5n']
                            elif(name == "ddG-of-insertion"):
                                peptide_list = [
                                    'GL5']#'GL6', 'GL7', 'GL8', 'GWL6']#'GL5', 
                            elif(name == "ddG-of-pH-insertion"):
                                peptide_list = ['pHLIP-v1_4', 'pHLIP-v1_8', 'pHLIP-v2_4', 'pHLIP-v2_8', 'pHLIP-v3_4', 'pHLIP-v3_8', 'pHLIP-v4_4', 'pHLIP-v4_8', 'pHLIP-v5_4', 'pHLIP-v5_8', 'pHLIP-v6_4', 'pHLIP-v6_8', 'pHLIP-v7_4', 'pHLIP-v7_8', 'pHLIP-v8_4', 'pHLIP-v8_8',
                                                'pHLIP-v9_4', 'pHLIP-v9_8', 'pHLIP-v10_4', 'pHLIP-v10_8', 'pHLIP-v11_4', 'pHLIP-v11_8', 'pHLIP-v12_4', 'pHLIP-v12_8', 'pHLIP-v13_4', 'pHLIP-v13_8', 'pHLIP-v14_4', 'pHLIP-v14_8', 'pHLIP-v15_4', 'pHLIP-v15_8', 'pHLIP-v16_4', 'pHLIP-v16_8']
                            elif(name == "ph-titration"):
                                peptide_list = ['veeks']#,'wt_W5_K12','wt_W5_K14','wt_Y5','wt_Y5_K12','wt_Y5_K14']#'polylys',
                            # combiningallfiles(peptide_list, Options.energy_fxn)
                            # creating_1d_min_files(peptide_list)
                            # combiningchargestatefiles(peptide_list, Options.energy_fxn,100)
                            # combiningchargestateresfiles(peptide_list, Options.energy_fxn,100)
                            # extract_pka(peptide_list)
                            # map_min_energy_pH(peptide_list)
                            
                            # extractingchargestateresfiles(peptide_list, Options.energy_fxn,100,17)
                            extract_pka_perres(peptide_list,17)
                            
                            # extractingchargestateresfiles(peptide_list, Options.energy_fxn,100,12)
                            extract_pka_perres(peptide_list,12)
                    
if __name__ == "__main__":
    main(sys.argv)
