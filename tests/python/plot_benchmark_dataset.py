#!/usr/bin/env python
""" Master script for plotting all benchmark data

This module plots for some of the tests.  
This is the last
step after generating all dataset.  

Authors: 
	Rituparna Samanta <rsamant2@jhu.edu>

Example: 
	

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
	- Rosetta release 246 or greater
	- PyRosetta4 for Python 3.6 or 3.7
"""

import sys, os, random
import hpc_util, read_config

import pandas as pd
import numpy as np
import plotly.graph_objects as go

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "sc-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "tm-peptide-tilt-angle", "helix-kinks", "hydrophobic-length", "adsorbed-peptide-tilt-angle", "protein-protein-docking", "protein-tilt-angle", "sequence-recovery", "amphiphilic-peptide-tilt-angle", "ph-titration"]


def main( args ): 
        import plotly.graph_objects as go
        # Read options from the command line
        parser = OptionParser(usage="usage %prog --energy_fxn franklin2023 --which_tests all --restore_talaris false" )
        parser.set_description(main.__doc__)

        parser.add_option( '--energy_fxn', '-e', action="store", help="Name of energy function weights file", )
        parser.add_option( '--which_tests', '-w', action="store", help="Specify tests run (comma separated list)", )

        (options, args) = parser.parse_args(args=args[1:])
        global Options
        Options = options

        # Check that required options have been provided
        if ( not Options.energy_fxn or not Options.which_tests ): 
                print("Missing required options --energy_fxn and/or --which_tests" )
                sys.exit()

        # Set restore variable based on energy function type
        restore = True
        if ( Options.energy_fxn == "franklin2019" or Options.energy_fxn == "ref2015" or Options.energy_fxn == "franklin2021" ): 
                restore = False 

        # Read path configuration file
        config = read_config.read_config()

        # Check test categories
        test_names = []
        if ( Options.which_tests == "all" ): 
                test_names = all_tests
        else: 
                test_names = Options.which_tests.split(",")
                print(test_names)
                # check that all names are valid
                for name in test_names: 
                        if name not in all_tests: 
                                sys.exit( "No such test " + name + ". Exiting!" )

        # wt_fa_water_to_bilayer = [1.0]#[1.0,1.0,0.5]#, 1.484] 
        # wt_fa_imm_elec = [1.0]#[0.01,1.0,1.0]#, 0.571]
        # wt_f_elec_bilayer = [0.172]#[0.128,0.172,1.0]#, -0.462]
        
        wt_fa_water_to_bilayer = [1.0]#[0.863]#, 1.484] 
        wt_fa_imm_elec = [0.01]#[0.001]#, 0.571]
        wt_f_elec_bilayer = [0.128]#[0.152]#, -0.462] 
        # Test #1: Tilt angles for transmembrane peptides
        for name in test_names: 
                if ( "adsorbed-peptide-tilt-angle" in test_names ):
                        #test 1a: A1+A3
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_min.dat'
                        comparefile = '../../franklin2019_summary.dat'
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        print(df_final)
                        plotly_scatter(df_final,"adsorbed_petide_p1_label_absolutetilts.png", datadir, labels=True)
                        plotly_scatter(df_final,"adsorbed_petide_p1_absolutetilts.png", datadir, labels=False) 
                elif ( "tm-peptide-tilt-angle" in test_names ):
                        #test 1a: A2+A4
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_min.dat'
                        comparefile = '../../franklin2019_summary.dat'
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        #     print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        #     print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        #     print(df_final)  
                        df_final = normalize_angle(df_final,'tilt_f21')  
                        # plotly_scatter(df_final,"tm_petide_p1_label_2022.png", datadir, labels=True)
                        # plotly_scatter(df_final,"tm_petide_p1_2022.png", datadir)
                        
                        #test1b: polyalanines: polyalanines are blocked with amidateation in ctrem and methylated in Nterm. 
                        #also like the paper, I restrict it at z=0. 
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_AA_min.dat'
                        comparefile = '../../franklin2019_summary.dat'
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        #     print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        #     print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        
                        # plotly_scatter(df_final,"tm_petide_AA_mp.png", datadir,compare_mp=True,labels=False)
                        # plotly_scatter(df_final,"tm_petide_AA_label.png", datadir, labels=True)
                        # plotly_scatter(df_final,"tm_petide_AA.png", datadir)
                        
                        #test1c: WALP: leaving teh charged ends. 
                        #I tried acetylation and amidation. acetylated end is increasing the del G at the end weirdly. 
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG_WALP_dmpc/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_WALP_min.dat'
                        comparefile = '../../franklin2019_summary_WALPdmpc.dat'
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        #     print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        #     print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        #     print(df_final)   
                        plotly_scatter(df_final,"tm_petide_WALP_mp.png", datadir,compare_mp=True, labels=False)
                        # plotly_scatter(df_final,"tm_petide_WALP_labels.png", datadir, compare_mp=False, labels=True)
                        # plotly_scatter(df_final,"tm_petide_WALP.png", datadir, compare_mp=False)
                        
                elif("ddG-of-mutation" in test_names):
                        #test 1a: A1+A3
                        datadir = "/home/rsamant2/scratch16-jgray21/rsamant2/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/" + "data/"+ \
                                Options.energy_fxn + "/" + name + "/test8_afternheavyatomcorrection_changedddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'C1_OmpLA_canonical_ddGs/ddG_franklin2021_allruns_grouped.dat'
                        comparefile = 'ompla_franklin2019_allruns_grouped.dat'
                        #to group the data
                        #remove the headers
                        #df_ddgmutation_ompla=pd.read_csv('ddG_franklin2021.dat',sep=" ")
                        #df_new=df_ddgmutation_ompla.groupby(['Nat','Pos','Mut','experimental_ddG','class','depth'])['predicted_ddG'].agg([('mean_predictied_ddG','mean'),('std_dev','std')]).reset_index()
                        # df_new.to_csv('ddG_franklin2021_allruns_grouped.dat',sep=" ",index=False)

                        df_data = pd.read_csv(datafile, delimiter=' ')
                        print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on=['Nat','Pos','Mut','experimental_ddG','class','depth'])
                        print(df_final)
                        plotly_scatter_colorbyclass(df_final, "ompla_dd_werror.png", datadir,err=True)

                        datafile = 'C2_PagP_canonical_ddGs/ddG_franklin2021_allruns_grouped.dat'
                        comparefile = 'pagp_franklin2019_allruns_grouped.dat'
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on=['Nat','Pos','Mut','experimental_ddG','class','depth'])
                        print(df_final)
                        plotly_scatter_colorbyclass(df_final, "pagp_ddG_werror.png", datadir,err=True)
            
                elif("ddG-of-insertion" in test_names):
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changedddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_ddG_min.dat'
                        comparefile = '../../franklin2019_ddg.dat'
                        experimentfile = '../../experiment_test5.dat'
                        
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        df_f21 = df_data[['seq','exp.ddG','predicted.ddG']]
                        # df_f21['sfxn']='Franklin2021'
                        df_f21['sfxn']='F23'
                        
                        print(df_f21)
                        
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        df_experiment = pd.read_csv(experimentfile,delimiter=' ')
                        df_f19 = pd.merge(df_experiment,df_compare, on='seq')
                        df_f19.columns=['seq','exp.ddG','predicted.ddG']
                        # df_f19['sfxn']='Franklin2019'
                        df_f19['sfxn']='F19'
                        print(df_f19)
                        
                        df_final=pd.concat([df_f19,df_f21])
                        print(df_final)        
                        plotly_scatter_withbestfit(df_final, "ddG_of_insertion.png", datadir)
                        
                elif("protein-tilt-angle" in test_names):
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        Options.energy_fxn + "/" + name + "/weights_from_test8_heavyatomcorrection_notilt_changeddG/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                        wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'summary_min.dat'
                        comparefile = '../../franklin2019_predictions.dat'
                        experimentfile = '../../opm_predictions.dat'
                        
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        #     print(df_data)
                        df_compare = pd.read_csv(comparefile,delimiter=' ')
                        #     print(df_compare)
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        df_exp = pd.read_csv(experimentfile, delimiter=' ')
                        df_final = pd.merge(df_final,df_exp[['pdb','delcenter']])
                        #     print(df_final)
                        
                        plotly_scatter_multipassprotein(df_final, "f21_tiltprediction.png", datadir, 'exp tilt', 'tilt_f21', 'tilt_f19', 'circle')
                        plotly_scatter_multipassprotein(df_final, "f19_tiltprediction.png", datadir, 'exp tilt','tilt_f19', 'tilt_f21', 'square')
                                        
                        plotly_scatter_multipassprotein(df_final, "f21_depthprediction.png", datadir, 'exp depth','depth_f21', 'depth_f19', 'circle')
                        plotly_scatter_multipassprotein(df_final, "f19_depthprediction.png", datadir, 'exp depth','depth_f19', 'depth_f21', 'square')
                        
                        #generate residual plot
                        df_final = normalize_angle(df_final, 'exp tilt')
                        df_final = normalize_angle(df_final, 'tilt_f19')
                        df_final = normalize_angle(df_final, 'tilt_f21')
                        
                        #     print(df_final)
                        #calculate_residuals
                        # have taken absolute value of the depths, since the system is symmetric about the 
                        #center of the membrane and 180 degrees. a protein at x degrees and depth = y is same as at 180-x degrees and depth = -y. 
                        #Rosetta calculated center of the protein for rotation whereas OPM of the TM section. 
                        df_final['depth.res_f21']=  abs(df_final['exp depth']) - abs(df_final['depth_f21'] + df_final['delcenter'] )
                        df_final['depth.res_f19']=  abs(df_final['exp depth']) - abs(df_final['depth_f19'] + df_final['delcenter'] )
                
                        df_final['tilt.res_f21'] =  ( abs(df_final['exp tilt']) - abs(df_final['tilt_f21']) )
                        df_final['tilt.res_f19'] =  ( abs(df_final['exp tilt']) - abs(df_final['tilt_f19']) )
                        
                        print(df_final)
                        
                        plot_residue_ecdf(df_final, "cdf_compare_depth.png", datadir, 'depth.res_f21','depth.res_f19')
                        plot_residue_ecdf(df_final, "cdf_compare_tilt.png", datadir, 'tilt.res_f21','tilt.res_f19')
                        
                elif("sequence-recovery" in test_names):
                        datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                Options.energy_fxn + "/" + name + "/weights_from_test7_heavyatomcorrection_disallowC/"
                        datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                                wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        print(datadir)
                        if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                        os.chdir(datadir)
                        datafile = 'sequence_recovery_per_pdb.csv'
                        comparefile =  '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                        'franklin2019/' + name + '/sequence_recovery_per_pdb.csv'
                        
                        df_data = pd.read_csv(datafile, delimiter=' ')
                        df_data.columns=['pdb','seq_rec_f21']
                        df_compare = pd.read_csv(comparefile, delimiter=' ')
                        df_compare.columns=['pdb','seq_rec_f19']
                        df_final = pd.merge(df_data,df_compare, on='pdb')
                        print(df_final)
                        alpha_list = config.benchmark_path + "targets/" +"design/alpha_monomer_chains.list"
                        beta_list = config.benchmark_path + "targets/" +"design/beta_monomer_chains.list"
                        
                        df_alpha= pd.read_csv(alpha_list)
                        df_alpha.columns=['pdb']
                        df_alpha['type'] = 'alpha'
                        
                        df_beta = pd.read_csv(beta_list)
                        df_beta.columns=['pdb']
                        df_beta['type'] = 'beta'

                        df_type=pd.concat([df_alpha,df_beta])
                        df_final = pd.merge(df_final,df_type, on='pdb')
                        print(df_final)
                        alpha = df_final[df_final['type']=='alpha'].sort_values(by='seq_rec_f21',ascending=False)
                        print(alpha)
                        alpha.to_csv('alpha_sorted_list.csv', index=False, sep=' ')
                        beta = df_final[df_final['type']=='beta'].sort_values(by='seq_rec_f21',ascending=False)
                        beta.to_csv('beta_sorted_list.csv', index=False, sep=' ')
                
                elif("sc-distribution" in test_names):
                        # datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                        #         Options.energy_fxn + "/" + "sequence-recovery"+ "/weights_from_test7_heavyatomcorrection_disallowC/"
                        # datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[0]) + "_felecbilayer_" + str(
                        #         wt_f_elec_bilayer[0]) + "_fimm_" + str(wt_fa_imm_elec[0])
                                
                        # print(datadir)
                        # if (not os.path.isdir(datadir)):
                        #         sys.exit("No such test data available; test needs to finish before combining files")
                        # os.chdir(datadir)
                        
                        # alpha_datafile = 'alpha_sorted_list.csv'
                        # beta_datafile = 'beta_sorted_list.csv'
                        # df_alpha = pd.read_csv(alpha_datafile, delimiter=' ')
                        # df_beta = pd.read_csv(beta_datafile, delimiter=' ')
                        # df_final = pd.concat([df_alpha,df_beta])
                        
                        # compare_file = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/mpnn_seq_recov.dat'
                        # df_compare = pd.read_csv(compare_file,delimiter=' ')
                        # df_final = pd.merge(df_final,df_compare, on='pdb')
                        # print(df_final)
                        # df_final.to_csv('/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/full_seq_recov.dat', sep=' ', index=False)
                        # plot_distribution_persfxn(df_final,'seq_rec_f21','franklin2021_seqrecovery_plot.png','Franklin2021 Sequence Recovery')
                        # plot_distribution_persfxn(df_final,'seq_rec_f19','franklin2019_seqrecovery_plot.png','Franklin2021 Sequence Recovery')
                        
                        # print("==================================")
                        # print(df_beta)
                        # df_beta = pd.merge(df_beta,df_compare, on='pdb')
                        # print(df_beta)
                        # df_beta.to_csv('/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/beta_full_seq_recov.dat', sep=' ', index=False)
                        # plot_distribution_persfxn(df_beta,'seq_rec_f21','beta_franklin2021_seqrecovery_plot.png','Franklin2021 Sequence Recovery')
                        # plot_distribution_persfxn(df_beta,'seq_rec_f19','beta_franklin2019_seqrecovery_plot.png','Franklin2021 Sequence Recovery')
                        
                        print("============hydropathy plot===========")
                        franklin_file = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/weights_from_test7_heavyatomcorrection_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001/franklin_hydropathybeta.dat'
                        compare_file = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/beta_hydropathy.dat'
                        df_compare = pd.read_csv(compare_file,delimiter=' ')
                        df_franklin = pd.read_csv(franklin_file,delimiter=' ')
                        df_franklin['model'] = 'franklin21'
                        df_mpnn = df_compare[['pdb','avg_aqueous_hydropathy','avg_lipid_hydropathy','avg_interface_hydropathy']]
                        df_mpnn['model'] = 'mpnn'
                        df_native = df_compare[['pdb','native_aqueous_hp','native_lipid_hp','native_interface_hp']]
                        df_native.columns = ['pdb','avg_aqueous_hydropathy','avg_lipid_hydropathy','avg_interface_hydropathy']
                        df_native['model'] = 'native'
                        df_final = pd.concat([df_franklin,df_mpnn,df_native])
                        print(df_final)
                        
                        df_final2=pd.DataFrame()
                        for col_name in ['avg_aqueous_hydropathy','avg_lipid_hydropathy','avg_interface_hydropathy']:
                                df_final_temp = df_final[['pdb',col_name,'model']]
                                df_final_temp.columns = ['pdb','hydropathy','model']
                                df_final_temp['type'] = col_name.split('_')[1]
                                
                                df_final2 = pd.concat([df_final2,df_final_temp])
                        print(df_final2)
                        barplot_hydropathy_persfxn(df_final2, '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/sequence-recovery/protein_mpnn/no_design_C/beta_franklin2019_hydropathy_plot.png')
                        
def normalize_angle(df,column_name):
        df['diff']=df[column_name]
        print(df['diff'])
        df[(df[column_name]>90)]['diff'] = 180-df[(df[column_name]>90)][column_name]
        df[(df[column_name]<-90)]['diff'] = -180-df[(df[column_name]>90)][column_name]
        print(df['diff'])
        df[column_name] = df['diff']
        print(df)
        df.drop(['diff'], axis=1)
        return(df)
        
def plotly_scatter(df_final, filename, datadir, compare_mp=False, labels=False, adsorbed=False):
        fig = go.Figure()
        if(compare_mp):
                comparefile = '../../mp15_summary.dat'
                df_compare = pd.read_csv(comparefile,delimiter=' ')
                df_final = pd.merge(df_final,df_compare, on='pdb')
        
        if(adsorbed):
                #here the angles vary from -90 to +90 about the surface. -90 meaning along z-axis. 
                #+90 meaning opposite to z-axis. they have all been transformed to lower membrane surface. 
                max_lim = float(max(max(df_final['tilt_f19']),max(df_final['tilt_f21']),max(df_final['exp tilt']))) + 2
                min_lim = float(min(min((df_final['tilt_f19'])),min((df_final['tilt_f21'])),min((df_final['exp tilt'])))) - 2 
       
        else:
                if(compare_mp):
                        max_lim = float(max(max(df_final['tilt_f19']),max(df_final['tilt_f21']),max(df_final['exp tilt']),max(df_final['tilt_mp15']))) + 2
                        min_lim = float(min(min(abs(df_final['tilt_f19'])),min(abs(df_final['tilt_f21'])),min(abs(df_final['exp tilt'])),min(abs(df_final['tilt_mp15'])))) 
                else:
                        max_lim = float(max(max(df_final['tilt_f19']),max(df_final['tilt_f21']),max(df_final['exp tilt']))) + 2
                        min_lim = float(min(min(abs(df_final['tilt_f19'])),min(abs(df_final['tilt_f21'])),min(abs(df_final['exp tilt'])))) 

        print(max_lim)
        print(min_lim)
        x=np.linspace(min_lim, max_lim, num=5)
        x_10 = np.linspace(min_lim, max_lim-10, num=5)
        x_20 = np.linspace(min_lim, max_lim-20, num=5)
        x_m10 = np.linspace(min_lim+10, max_lim, num=5)
        x_m20 = np.linspace(min_lim+20, max_lim, num=5)
        
        y = x
        y_10 = x_10+10
        y_m10 = x_m10-10
        y_20 = x_20+20
        y_m20 = x_m20-20
        
        # Add traces
        if(labels):
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=abs(df_final['tilt_f21'][:-2]),text=list(df_final['pdb'][:-2]),
                        mode='markers+text',marker_symbol='circle',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", marker_line_width=1.0,marker_size=10,
                        textposition='top left',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),
                        name='F23'))
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=abs(df_final['tilt_f19'][:-2]),text=list(df_final['pdb'][:-2]),
                        mode='markers+text',marker_symbol='square',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=1.0, marker_size=10,
                        textposition='bottom left',textfont=dict(family="Helvetica",size=15,color="midnightblue"),
                        name='F19'))
                #for WALP
                # fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=abs(df_final['tilt_f21'][-2:]),text=list(df_final['pdb'][-2:]),
                #         mode='markers+text',marker_symbol='circle-open',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", marker_line_width=2.0,marker_size=10,
                #         textposition='top right',textfont=dict(family="Helvetica",size=15,color="rgb(240,96,96)"),name='F23',showlegend=False))
                # fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=abs(df_final['tilt_f19'][-2:]),text=list(df_final['pdb'][-2:]),
                #         mode='markers+text',marker_symbol='square-open',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=2.0, marker_size=10,
                #         textposition='bottom right',textfont=dict(family="Helvetica",size=15,color="midnightblue"),
                #         name='F19', showlegend=False))
                
        
                if(compare_mp):
                        fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=df_final['tilt_mp15'][:-2],
                                mode='markers+text',marker_symbol='star',marker_color="rgb(123, 110, 250)",marker_line_color="black", marker_line_width=1.0,marker_size=10,
                                name='M12(IMM1)'))
                        #for WALP
                        # fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=df_final['tilt_mp15'][-2:],
                        #         mode='markers+text',marker_symbol='star-open',marker_color="rgb(123, 110, 250)",marker_line_color="black", marker_line_width=2.0,marker_size=10,
                        #         name='M12(IMM1)',showlegend=False))
        else:
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=df_final['tilt_f21'][:-2],
                        mode='markers',marker_symbol='circle',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", marker_line_width=1.0,marker_size=10,
                        name='F23'))
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=df_final['tilt_f19'][:-2],
                        mode='markers',marker_symbol='square',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=1.0, marker_size=10,
                        name='F19'))
                #the later commented part must ne removed for WALP
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=abs(df_final['tilt_f21'][-2:]),text=list(df_final['pdb'][-2:]),
                        mode='markers',marker_symbol='circle-open',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", marker_line_width=2.0,marker_size=10,
                        name='F23', showlegend=False))
                fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=abs(df_final['tilt_f19'][-2:]),text=list(df_final['pdb'][-2:]),
                        mode='markers',marker_symbol='square-open',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=2.0, marker_size=10,
                        name='F19', showlegend=False))
        
                if(compare_mp):
                        fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][:-2]), y=df_final['tilt_mp15'][:-2],
                                mode='markers',marker_symbol='star',marker_color="rgb(123, 110, 250)",marker_line_color="black", marker_line_width=1.0,marker_size=10,
                                name='M12(IMM1)'))
                        #for WALP
                        fig.add_trace(go.Scatter(x=abs(df_final['exp tilt'][-2:]), y=df_final['tilt_mp15'][-2:],
                                mode='markers',marker_symbol='star-open',marker_color="rgb(123, 110, 250)",marker_line_color="black", marker_line_width=1.0,marker_size=10,
                                name='M12(IMM1)', showlegend=False))
                
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=4), showlegend=False))
        fig.add_trace(go.Scatter(x=x_10, y=y_10,mode='lines',
                        line = dict(color='black', width=4,dash='dash'),showlegend=False))
        fig.add_trace(go.Scatter(x=x_20, y=y_20,mode='lines',
                        line = dict(color='black', width=4,dash='dot'),showlegend=False)) 
        fig.add_trace(go.Scatter(x=x_m10, y=y_m10,mode='lines',
                        line = dict(color='black', width=4,dash='dash'),showlegend=False)) 
        fig.add_trace(go.Scatter(x=x_m20, y=y_m20,mode='lines',
                        line = dict(color='black', width=4,dash='dot'),showlegend=False))           

        #courtesy: Ameyas plots for docking .
        fig.update_xaxes(range=[min_lim, max_lim], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Helvetica", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_yaxes(range=[min_lim, max_lim],tickfont=dict(family="Helvetica", color='black',size=20),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_layout(
                font=dict(family="Helvetica", color='black',size=20))
        
        fig.update_layout(xaxis=dict(tick0 = 0, dtick=10))
        fig.update_layout(yaxis=dict(tick0 = 0, dtick=10))
        
        fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        # fig.update_yaxes(tickmode="auto", nticks=7)
        # fig.update_traces(marker=dict(size=20,line=dict(width=0.4,color='Black')))
        fig.update_layout(showlegend=True,legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
                ))
        fig.update(layout_coloraxis_showscale=True)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        # fig.update_layout( yaxis2=dict( title = pdbid, titlefont=dict( family="Helvetica", size=24), side="right" ) )
        
        fig.write_image( os.path.join(datadir,filename ), format="png", width=550, height=550, scale=2)
        # plotly.io.write_image(fig, file, format=None, scale=None, width=None, height=None, validate=True, engine='auto')

def plotly_scatter_multipassprotein(df_final, filename, datadir, exp_column_name, column_name, compare_column_name, symbol):
        fig = go.Figure()
        max_lim = float(max(max(df_final[column_name]),max(df_final[compare_column_name]),max(df_final[exp_column_name])))
        min_lim = float(min(min(abs(df_final[column_name])),min(abs(df_final[compare_column_name])),min(abs(df_final[exp_column_name]))))
        print(max_lim)
        print(min_lim)
        
        x=np.linspace(min_lim, max_lim, num=5)
        x_10 = np.linspace(min_lim, max_lim-10, num=5)
        x_20 = np.linspace(min_lim, max_lim-20, num=5)
        x_m10 = np.linspace(min_lim+10, max_lim, num=5)
        x_m20 = np.linspace(min_lim+20, max_lim, num=5)
        
        y = x
        y_10 = x_10+10
        y_m10 = x_m10-10
        y_20 = x_20+20
        y_m20 = x_m20-20
        

        # Add traces
        fig.add_trace(go.Scatter(x=abs(df_final[exp_column_name]), y=abs(df_final[df_final['alpha_or_beta']=='α'][column_name]),
                        mode='markers',marker_symbol=symbol,marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)", marker_line_width=1.0,marker_size=10,
                        name='α'))
        fig.add_trace(go.Scatter(x=abs(df_final[exp_column_name]), y=abs(df_final[df_final['alpha_or_beta']=='β'][column_name]),
                        mode='markers',marker_symbol=symbol,marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=1.0, marker_size=10,
                        name='β'))
        
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=4), showlegend=False))
        fig.add_trace(go.Scatter(x=x_10, y=y_10,mode='lines',
                        line = dict(color='gray', width=4,dash='dash'),showlegend=False))
        fig.add_trace(go.Scatter(x=x_20, y=y_20,mode='lines',
                        line = dict(color='gray', width=4,dash='dot'),showlegend=False)) 
        fig.add_trace(go.Scatter(x=x_m10, y=y_m10,mode='lines',
                        line = dict(color='gray', width=4,dash='dash'),showlegend=False)) 
        fig.add_trace(go.Scatter(x=x_m20, y=y_m20,mode='lines',
                        line = dict(color='gray', width=4,dash='dot'),showlegend=False))           

        #courtesy: Ameyas plots for docking .
        fig.update_xaxes(range=[min_lim, max_lim], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Helvetica", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_yaxes(range=[min_lim, max_lim],tickfont=dict(family="Helvetica", color='black',size=20),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_layout(
                font=dict(family="Helvetica", color='black',size=20))
        
        fig.update_layout(xaxis=dict(tick0 = 0, dtick=10))
        fig.update_layout(yaxis=dict(tick0 = 0, dtick=10))
        
        fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        # fig.update_yaxes(tickmode="auto", nticks=7)
        # fig.update_traces(marker=dict(size=20,line=dict(width=0.4,color='Black')))
        fig.update_layout(showlegend=True,legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
                ))
        fig.update(layout_coloraxis_showscale=True)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        # fig.update_layout( yaxis2=dict( title = pdbid, titlefont=dict( family="Helvetica", size=24), side="right" ) )
        
        fig.write_image( os.path.join(datadir,filename ), format="png", width=550, height=550, scale=2)
        # plotly.io.write_image(fig, file, format=None, scale=None, width=None, height=None, validate=True, engine='auto')

def plotly_scatter_colorbyclass(df_final, filename, datadir,err=False):
       
        print(list(df_final['Mut']))
        list_res=df_final['Mut'].unique()
        list_filters=list([x for x in list_res if x not in ['P','A','D','E']])
        df_final_temp=df_final[df_final['Mut'].isin(list_filters)]
        df_final=df_final_temp

        fig = go.Figure()
        if(err):
                max_lim = float(max(max(df_final['mean_predicted_ddG'])+1,max(df_final['mean_ddG_f19'])+1,max(df_final['experimental_ddG'])))
                min_lim = float(min(min(df_final['mean_predicted_ddG'])+1,min(df_final['mean_ddG_f19'])+1,min(df_final['experimental_ddG'])))
        else:
                max_lim = float(max(max(df_final['predicted_ddG']),max(df_final['ddG_f19']),max(df_final['experimental_ddG'])))
                min_lim = float(min(min(df_final['predicted_ddG']),min(df_final['ddG_f19']),min(df_final['experimental_ddG'])))
        print(max_lim)
        print(min_lim)
        x=np.linspace(min_lim, max_lim, num=5)
        # x_10 = np.linspace(min_lim, max_lim-10, num=5)
        # x_20 = np.linspace(min_lim, max_lim-20, num=5)
        # x_m10 = np.linspace(min_lim+10, max_lim, num=5)
        # x_m20 = np.linspace(min_lim+20, max_lim, num=5)
        
        y = x
        # y_10 = x_10+10
        # y_m10 = x_m10-10
        # y_20 = x_20+20
        # y_m20 = x_m20-20
        
        # Add traces
        import plotly.express as px
        
        if(err):
                fig.add_trace(go.Scatter(x=df_final['experimental_ddG'], y=df_final['mean_predicted_ddG'],error_y=dict(type='data',array=df_final['std_dev'],visible=True),text=list(df_final['Mut']),
                        mode='markers+text',marker_symbol='circle',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)",marker_line_width=1.0,marker_size=10, 
                        textposition='bottom center',textfont=dict(family="Helvetica",size=15,color="rgb(130, 3, 51)"),
                        name='F23'))
        
                fig.add_trace(go.Scatter(x=df_final['experimental_ddG'], y=df_final['mean_ddG_f19'],error_y=dict(type='data',array=df_final['std_dev_f19'],visible=True),text=list(df_final['Mut']),
                        mode='markers+text',marker_symbol='square',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=1.0,marker_size=10, 
                        textposition='top center',textfont=dict(family="Helvetica",size=15,color="rgb(2, 8, 115)"),
                        name='F19'))
        else:
                fig.add_trace(go.Scatter(x=df_final['experimental_ddG'], y=df_final['predicted_ddG'],text=list(df_final['Mut']),
                        mode='markers+text',marker_symbol='circle',marker_color="rgb(240,96,96)",marker_line_color="rgb(191,191,191)",marker_line_width=1.0,marker_size=10, 
                        textposition='bottom center',textfont=dict(family="Helvetica",size=15,color="rgb(130, 3, 51)"),
                        name='F23'))
        
                fig.add_trace(go.Scatter(x=df_final['experimental_ddG'], y=df_final['ddG_f19'],text=list(df_final['Mut']),
                        mode='markers+text',marker_symbol='square',marker_line_color="midnightblue", marker_color="lightskyblue",marker_line_width=1.0,marker_size=10, 
                        textposition='top center',textfont=dict(family="Helvetica",size=15,color="rgb(2, 8, 115)"),
                        name='F19'))
        
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=4), showlegend=False))
        # fig.add_trace(go.Scatter(x=x_10, y=y_10,mode='lines',
        #                 line = dict(color='Gray', width=4,dash='dash'),showlegend=False))
        # fig.add_trace(go.Scatter(x=x_20, y=y_20,mode='lines',
        #                 line = dict(color='Gray', width=4,dash='dot'),showlegend=False)) 
        # fig.add_trace(go.Scatter(x=x_m10, y=y_m10,mode='lines',
        #                 line = dict(color='Gray', width=4,dash='dash'),showlegend=False)) 
        # fig.add_trace(go.Scatter(x=x_m20, y=y_m20,mode='lines',
        #                 line = dict(color='Gray', width=4,dash='dot'),showlegend=False))           

        #courtesy: Ameyas plots for docking .
        fig.update_xaxes(range=[min_lim, max_lim], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Helvetica", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_yaxes(range=[min_lim, max_lim],tickfont=dict(family="Helvetica", color='black',size=20),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        # fig.update_traces(textposition='top right',textfont_size=14)
        fig.update_layout(
                font=dict(family="Helvetica", color='black',size=20))
        
        fig.update_layout(xaxis=dict(tick0 = round(min_lim), dtick=2.0))
        fig.update_layout(yaxis=dict(tick0 = round(min_lim), dtick=2.0))
        
        fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        # fig.update_yaxes(tickmode="auto", nticks=7)
        # fig.update_traces(marker=dict(size=20,line=dict(width=0.4,color='Black')))
        fig.update_layout(showlegend=True,legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
                ))
        fig.update(layout_coloraxis_showscale=True)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        # fig.update_layout( yaxis2=dict( title = pdbid, titlefont=dict( family="Helvetica", size=24), side="right" ) )
        
        fig.write_image( os.path.join(datadir,filename ), format="png", width=550, height=550, scale=2)
        # plotly.io.write_image(fig, file, format=None, scale=None, width=None, height=None, validate=True, engine='auto')

def plotly_scatter_withbestfit(df_final, filename, datadir, labels=False):
        fig = go.Figure()
        y_max_lim = float(max(df_final['predicted.ddG']))
        y_min_lim = float(min(df_final['predicted.ddG']))
        x_max_lim = float(max(df_final['exp.ddG']))
        x_min_lim = float(min(df_final['exp.ddG']))
        
        # Add traces
        import plotly.express as px
        
        color_discrete_map = {'F19': "lightskyblue", 'F23': "rgb(240,96,96)"}
        symbol_map = {'F19': "square", 'F23': "circle"}
        
        fig = px.scatter(df_final, x="exp.ddG", y="predicted.ddG", symbol="sfxn", color="sfxn", trendline="ols", color_discrete_map=color_discrete_map,symbol_map=symbol_map )
        results = px.get_trendline_results(fig)
        print(results)
        results_f19 = results.iloc[0]["px_fit_results"]
        results_f21 = results.iloc[1]["px_fit_results"]
        # results = results.iloc[1]["px_fit_results"].params
        print("f19 results c,m {} rsqd {}".format(results_f19.params, results_f19.rsquared))
        print("f19 results c,m {} rsqd {}".format(results_f21.params, results_f21.rsquared))
        
        # results.query("sfxn =='Franklin2019'").px_fit_results.iloc[0].summary()
        
        #courtesy: Ameyas plots for docking .
        fig.update_xaxes(range=[x_min_lim, x_max_lim], mirror=True, linewidth=2, linecolor='black',
                     tickfont=dict(family="Helvetica", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        fig.update_yaxes(range=[y_min_lim, y_max_lim],tickfont=dict(family="Helvetica", color='black',size=20),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)
        
        # fig.update_traces(textposition='top right',textfont_size=14)
        fig.update_layout(
                font=dict(family="Helvetica", color='black',size=20))
        
        fig.update_layout(xaxis=dict(title="",tick0 = round(x_min_lim), dtick=(round(x_max_lim)-round(x_min_lim))/5))
        fig.update_layout(yaxis=dict(title="",tick0 = round(y_min_lim), dtick=2.0))
        
        fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        # fig.update_yaxes(tickmode="auto", nticks=7)
        
        fig.update_traces(marker=dict(size=10,line=dict(width=1.0,color='Black')))
        fig.update_layout(legend_title="",showlegend=True,legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
                ))
        fig.update(layout_coloraxis_showscale=True)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        # fig.update_layout( yaxis2=dict( title = pdbid, titlefont=dict( family="Helvetica", size=24), side="right" ) )
        
        fig.write_image( os.path.join(datadir,filename ), format="png", width=550, height=550, scale=2)
        # plotly.io.write_image(fig, file, format=None, scale=None, width=None, height=None, validate=True, engine='auto')

def calculate_ecdf(data):
        # https://stackoverflow.com/questions/24788200/calculate-the-cumulative-distribution-function-cdf-in-python
        
        data_sorted = np.sort(data)
        # calculate the proportional values of samples
        p = 1. * np.arange(len(data)) / (len(data) - 1)

        return(data_sorted,p)


def plot_residue_ecdf(df, filename, datadir, data_colname,compare_colname):
        import plotly.express as px

        x_max_lim = float(max(max(abs(df[data_colname])),max(abs(df[compare_colname]))))
        x_min_lim = float(min(min(abs(df[data_colname])),min(abs(df['depth.res_f19']))))
        y_max_lim = 1.0
        y_min_lim = 0.0
        
        fig = go.Figure()
        
        x,y = calculate_ecdf(abs(df[data_colname]))
        print('F23')
        print(x[23:25],y[23:25])
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=2), name='F23'))
        del x,y 
        
        x,y = calculate_ecdf(abs(df[df['alpha_or_beta']=='α'][data_colname]))
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='red', width=2), showlegend=False))
        del x,y
        
        x,y = calculate_ecdf(abs(df[df['alpha_or_beta']=='β'][data_colname]))
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='blue', width=2), showlegend=False))
        
        del x,y
        
        x,y = calculate_ecdf(abs(df[compare_colname]))
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='black', width=2, dash='dash'), name='F19'))
        
        print('F19')
        print(x[15:25],y[15:25])
        del x,y 
        
        x,y = calculate_ecdf(abs(df[df['alpha_or_beta']=='α'][compare_colname]))
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='red', width=2, dash='dash'), showlegend=False))
        del x,y
        
        x,y = calculate_ecdf(abs(df[df['alpha_or_beta']=='β'][compare_colname]))
        fig.add_trace(go.Scatter(x=x, y=y,mode='lines',
                        line = dict(color='blue', width=2, dash='dash'), showlegend=False))
        
        del x,y
        
        
        
        #courtesy: Ameyas plots for docking .
        fig.update_xaxes(range=[x_min_lim, x_max_lim], mirror=True, linewidth=2, linecolor='black',
                        tickfont=dict(family="Helvetica", color='black',size=20), ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)

        fig.update_yaxes(range=[y_min_lim, y_max_lim],tickfont=dict(family="Helvetica", color='black',size=20),
                        mirror=True, linewidth=2, linecolor='black',showgrid=True, ticks="outside",tickson="boundaries",ticklen=5,tickwidth=2)

        # fig.update_traces(textposition='top right',textfont_size=14)
        fig.update_layout(
                font=dict(family="Helvetica", color='black',size=20))

        fig.update_layout(xaxis=dict(tick0 = round(x_min_lim), dtick=(round(x_max_lim)-round(x_min_lim))/4))
        fig.update_layout(yaxis=dict(tick0 = round(y_min_lim), dtick=0.2))

        fig.update_xaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        fig.update_yaxes(showgrid=True, gridwidth=0.75, gridcolor='Gray')
        # fig.update_yaxes(tickmode="auto", nticks=7)
        # fig.update_traces(marker=dict(size=20,line=dict(width=0.4,color='Black')))
        fig.update_layout(showlegend=True,legend=dict(
                yanchor="bottom",
                y=0.01,
                xanchor="right",
                x=0.99
                ))
        fig.update(layout_coloraxis_showscale=True)
        fig.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        # fig.update_layout( yaxis2=dict( title = pdbid, titlefont=dict( family="Helvetica", size=24), side="right" ) )

        fig.write_image( os.path.join(datadir,filename ), format="png", width=550, height=550, scale=2)
        # plotly.io.write_image(fig, file, format=None, scale=None, width=None, height=None, validate=True, engine='auto')

def plot_distribution_persfxn(df,column_name,outfile,xlabel_input):
    import seaborn as sns
    import matplotlib
    import matplotlib.pyplot as plt
    
    theme = {'axes.grid': True,
        'grid.linestyle': '',
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        "font.weight": 'regular',
        'xtick.color': 'black',
        'ytick.color': 'black',
        "axes.titlesize": 20,
        "axes.labelsize": 18
    }
    
    matplotlib.rcParams.update(theme)
    #to have a distribution you can change kde to True
#     sns.histplot(data=df, x=column_name, bins=10, kde=True)
    sns.kdeplot(data=df, x='seq_rec_f19',bw_adjust=.5, color='blue')
    sns.kdeplot(data=df, x='seq_rec_f21',bw_adjust=.5, color='red')
    sns.kdeplot(data=df, x='avg_seq_recov_mpnn',bw_adjust=.5, color='green')
    
#     plt.xlabel(xlabel_input)
    plt.xlabel('Sequence Recovery per protein')
    plt.ylabel("Count")
    
    delta = (max(df[column_name])-min(df[column_name]))/(4.0*5.0)
    #mean of design
    plt.axvline(x=df['seq_rec_f21'].mean(), color='red', linewidth=3)
    plt.text(df['seq_rec_f21'].mean()-delta,3,'Franklin2021',color='red',rotation='vertical')
    
    plt.axvline(x=df['seq_rec_f19'].mean(), color='blue', linewidth=3, ls='--')
    plt.text(df['seq_rec_f19'].mean()+delta,3,'Franklin2019',color='blue',rotation='vertical')
    
    plt.axvline(x=df['avg_seq_recov_mpnn'].mean(), color='green', linewidth=3, ls='-.')
    plt.text(df['avg_seq_recov_mpnn'].mean()+delta,3,'<PMPNN>',color='green',rotation='vertical')
    
#     plt.axvline(x=df[df['type']=='alpha']['avg_seq_recov_mpnn'].mean(), color='green', linewidth=3,ls='-.')
#     plt.text(df[df['type']=='alpha']['avg_seq_recov_mpnn'].mean()+delta,5,'α-helices',color='green',rotation='vertical')
    
#     plt.axvline(x=df[df['type']=='beta']['avg_seq_recov_mpnn'].mean(), color='green', linewidth=3, ls='-.')
#     plt.text(df[df['type']=='beta']['avg_seq_recov_mpnn'].mean()+delta,5,'β-barrels',color='green',rotation='vertical')
    
    plt.tight_layout()
    plt.savefig(outfile, transparent=True, dpi=600)
    plt.close()

def barplot_hydropathy_persfxn(df, outfile):
    import seaborn as sns
    import matplotlib
    import matplotlib.pyplot as plt
    
    theme = {'axes.grid': True,
        'grid.linestyle': '',
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        "font.weight": 'regular',
        'xtick.color': 'black',
        'ytick.color': 'black',
        "axes.titlesize": 20,
        "axes.labelsize": 18
    }
    
    ax = sns.boxplot(x = df['type'],
                y = df['hydropathy'],
                hue = df['model'],
                palette = 'Set2',
                order=['aqueous','interface','lipid'])
    ax.set_xlabel('')
    ax.set_ylabel('Hydropathy')
    plt.legend(bbox_to_anchor=(0.5, 1.01), loc='lower center', borderaxespad=0)
    plt.tight_layout(rect=[0, 0.03, 1, 0.9])
    plt.yticks(np.arange(-4.0, 4.0, 0.5))
    plt.savefig(outfile, dpi=600, transparent=True)
    plt.close()
    
    
if __name__ == "__main__":
    main(sys.argv)
