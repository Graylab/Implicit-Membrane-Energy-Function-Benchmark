#!/usr/bin/env python
""" Predict orientation of multi-pass transmembrane proteins

This script runs the MembraneEnergyLandscapeSampler on each target
in the dataset. The result is a complete mapping of energies to 
peptide orientation, given as a function of depth and either tilt 
or rotation angle. This script generates data for Taest 3.

Authors: 
	Rebecca Alford <ralford3@jhu.edu> 
	Rituparna Samanta <rituparna@utexas.edu>
Example: 
	$ import make_peptide_energy_landscape
	$ make_protein_energy_landscape.run_peptide_energy_landscape_calc( 
	  energy_fxn, config, test_name, targets_dir, xml )

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- config: Container with path to benchmark and rosetta files
	- test_name: Name of benchmark test
	- targets_dir: Location of test cases
	- xml: Path to MembraneEnergyLandscapeSampler XML application

Requirements: 
	- Rosetta release 246 or greater
"""

import random
import sys, os
import hpc_util, read_config
from string import Template

def run_protein_energy_landscape_calc( energy_fxn, restore, config, test_name ): 
    """
	A function to generate orientation-dependent peptide energy landscapes

	Arguments: 
		energy_fxn = energy function to use for calculations (typically, name of the weights file)
		config = path to benchmark, Rosetta executables
		test_name = name of test directory
		xml = Path to RosettaScript defining the protein energy landscape search protocol
    """

    print( "Generating data for protein energy landscape test..." ) 

    # Read list of energy landscape test cases
    targets_path = config.benchmark_path + "targets/orientation/B1_multispan_proteins/"
    list_of_targets = targets_path + "target_master.list"
    #print(list_of_targets)

    with open( list_of_targets, 'rt' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate path to executable
    executable = config.rosetta_path + "/rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
    xml_script =  config.benchmark_path + "tests/xml/make_depth_vs_tilt_protein_energy.xml"
    print(xml_script)

    # Change directories to a data analysis dir
    #outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name +"/weights_from_test7"
    outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/" 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )
    
    outdir = outdir + energy_fxn + "/" 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )

    outdir = outdir + "/" + test_name 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )
   
    outdir = outdir + "/weights_from_test8_heavyatomcorrection_notilt_changeddG"
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )
    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = targets_path + case + "/" + case + "_after_relax.pdb"
        spanfile = "from_structure"

        if ( "protein-tilt-angle" in test_name ):
            interface = "0"
			
        #flag_axis=2:eigenvector; 0:vector joining top and bottom centers; 1:average of all tm-axis
        flag_axis = ["1.0"]
        print(flag_axis[0])


        if (interface=="0"):
            start_z = ["0","10","20","30","40","50"]
            end_z = ["10","20","30","40","50","60"]
		#fa_w_b    fa_imm_elec f_elec_lipidlayer 	
        #1.375		-0.143		0.286
        ##0.933		0.733		0.094
        #1.179		-0.015		-0.107
        #1.571		-0.025		-0.179
        #after fa_imm_modiefied to outside_memb fa_imm->0
        #1.629	    -1.561	    0.136
        #1.484	    0.571	    -0.462
        # wt_fa_water_to_bilayer = [0.863]#, 1.484] 
        # wt_fa_imm_elec = [0.001]#, 0.571]
        # wt_f_elec_bilayer = [0.152]#, -0.462]       
        
        wt_fa_water_to_bilayer = [0.863, 1.0]#, 1.484] 
        wt_fa_imm_elec = [0.001, 0.01]#, 0.571]
        wt_f_elec_bilayer = [0.152, 0.128]#, -0.462]        
        
        for i in range(len(wt_fa_water_to_bilayer)):
            # for k in range(len(wt_f_elec_bilayer)):
            #     for l in range(len(wt_fa_imm_elec)):
                              
            for j in range(len(start_z)):
    
                # Default composition is DLPC
                arguments = " -out:nooutput -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:script_vars interface=" + interface + " -parser::script_vars start_z=" + start_z[j] + " -parser::script_vars end_z=" + end_z[j] + " -parser::script_vars flag_axis=" + flag_axis[0] +  " -parser:protocol " + xml_script 
                
                #print(arguments)

                if ( restore ): 
                    arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior "
                else:
                    arguments = arguments + " -mp:lipids:composition DOPC"

                
                casedir = outdir + "/fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(\
                            wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])

                if ( not os.path.isdir( casedir ) ): 
                    os.system( "mkdir " + casedir )
                os.chdir( casedir )

                casedir = casedir + "/" + case
                if ( not os.path.isdir( casedir ) ): 
                    os.system( "mkdir " + casedir )
                os.chdir( casedir )

                if(energy_fxn=='franklin2021'):
                    if( wt_fa_imm_elec[i]<0 or wt_f_elec_bilayer[i]<0 ):
                        filename = 'flag_file'
                        #the line below should be "w" not "a" or append mode. 
                        with open( filename, 'w' ) as f: 
                            f.write(  "-set_weights fa_water_to_bilayer "+ str(wt_fa_water_to_bilayer[i]) + "\n" )
                            f.write(  "-set_weights f_elec_lipidlayer "+ str(wt_f_elec_bilayer[i]) + "\n" )
                            f.write(  "-set_weights fa_imm_elec "+ str(wt_fa_imm_elec[i]) + "\n" )

                        f.close()
                        os.system( "chmod +x " + filename )
                        arguments = arguments + " @flag_file"

                    else:
                        arguments = arguments + " -set_weights fa_water_to_bilayer " + \
                                    str(wt_fa_water_to_bilayer[i]) + " -set_weights f_elec_lipidlayer " + str(\
                                    wt_f_elec_bilayer[i]) + " -set_weights fa_imm_elec " + str(wt_fa_imm_elec[i])
                    
                

                # Write jobfile and submit to the HPC
                print("Submitting protein orientation sampling calculations for case:", case ) 
                jobname = case + "_" + end_z[j] + "_" + start_z[j] + "_energy_landscape"
                print(casedir)
                jobfile = hpc_util.make_jobfile( casedir, jobname, executable, arguments )
                
                hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile )
                #hpc_util.submit_condor_job( casedir, jobname, jobfile, "", 1 )
