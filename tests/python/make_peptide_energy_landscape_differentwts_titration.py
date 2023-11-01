#!/usr/bin/env python
""" Predict orientation of helical peptides

The main objective of the test is to do a paramter sweep on the
different weights of the score function. We then chose the 
weights combination which produces the minimum error. This script runs the MembraneEnergyLandscapeSampler on each target
in the dataset. The result is a complete mapping of energies to 
peptide orientation, given as a function of depth and tilt 
either at rotation angle=0 or minimized over all rotation 
angles. This script generates data for Tests 1,2,5,6.

Authors: 
	Rituparna Samanta <rituparna@utexas.edu>
Example: 
	$ import make_peptide_energy_landscape
	$ make_peptide_energy_landscape.run_peptide_energy_landscape_calc( 
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
import os
import sys
import hpc_util
import read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname(os.path.realpath(__file__))


def run_peptide_energy_landscape_calc(energy_fxn, config, test_name, targets, xml, pH_mode=False, pH=7):
    """
    A function to generate orientation-dependent peptide energy landscapes

    Arguments: 
            energy_fxn = energy function to use for calculations (typically, name of the weights file)
            config = path to benchmark, Rosetta executables
            test_name = name of test directory
            targets = Name of targets subdirectory (e.g. A4_designed_surface_ahelices)
            xml = Path to RosettaScript defining the protein energy landscape search protocol
    """

    print("Generating data for peptide tilt and rotation angle test...")

    # Read list of energy landscape test cases
    #targets_path = config.benchmark_path + "targets/" + targets + "/"
    #targets_path = config.benchmark_path + targets + "/"
    
    # if(pH==8):
    #         targets_path = config.benchmark_path + targets + "/"
    # elif(pH==4):
    #         targets_path = config.benchmark_path + targets + "/"
    # else:
    targets_path = config.benchmark_path + "targets/" + targets + "/"
    
    
    list_of_targets = targets_path + "targets.list"
    
    if(not os.path.exists(list_of_targets)):
        print("the list of targets do not exist")
        sys.exist()
        
    with open(list_of_targets, 'rt') as f:
        test_cases = f.readlines()
    test_cases = [x.strip() for x in test_cases]

    if(energy_fxn == "franklin2019"):
        executable = "/home/rsamant2/scratch16-jgray21/rsamant2/Rosetta/main/source/bin"
        executable = executable + "/rosetta_scripts" + \
            "." + config.platform + config.compiler + config.buildenv
    else:
        executable = config.rosetta_path + "/rosetta_scripts" + \
            "." + config.platform + config.compiler + config.buildenv
    xml_script = config.benchmark_path + xml
    print("xml_scrip:~", xml_script)

    # Change directories to a data analysis dir
    # outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name
    # print("outdir:~", outdir)

    outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/" 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )
    
    outdir = outdir + energy_fxn + "/" 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )

    outdir = outdir + test_name 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )
    outdir = outdir + "/weights_from_test7_titration"
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        # if(pH==8):
        #     pdbfile = targets_path + case + "/" + case + "_disordered.pdb"
            
        # elif(pH==4):
        #     pdbfile = targets_path + case + "/" + "ranked_0.pdb"
        # else:
        if("tm-peptide-tilt-angle" in test_name):
            pdbfile = targets_path + case + "/" + case + "_relaxed.pdb"
        elif("ph-titration" in test_name):
            pdbfile = targets_path + case + "/" + case + "_after_relax.pdb"
        else:
            pdbfile = targets_path + case + "/" + case + ".pdb"
        
        
        spanfile = "single_TM_mode"
        print("pdbfile:~", pdbfile)

        # flag_axis=2:eigenvector; 0:vector joining top and bottom centers; 1:average of all tm-axis
        flag_axis = ["0.0"]

        # Default composition is DLPC
        if ("tm-peptide-tilt-angle" in test_name):
            interface = "0"
        if ("adsorbed-peptide-tilt-angle" in test_name):
            interface = "1"
            if("tilt_angle/A4_design_surface_ahelices" not in targets):
                flag_axis = ["2.0"]

        if ("ddG-of-insertion" in test_name):
            interface = "0"
        if ("ddG-of-pH-insertion" in test_name):
            interface = "0"
        #similar to the adsorbed case. We want it to lie on the surface. 
        if("ph-titration" in test_name):
            interface = "1"
            if(case == "polylys"):
                flag_axis = ["2.0"]
        #fa_w_b    fa_imm_elec f_elec_lipidlayer 	
        #1.375		-0.143		0.286
        #1.571		-0.025		-0.179
        #after fa_imm_modiefied to outside_memb fa_imm->0
        #1.629	    -1.561	    0.136
        #1.484	    0.571	    -0.462

        wt_fa_water_to_bilayer = [0.966]#, 1.484] 
        wt_fa_imm_elec = [0.379]#, 0.571]
        wt_f_elec_bilayer = [0.016]#, -0.462]  
          

        start_z = ["0","55"]#["0","10","20","30","40","50"]
        end_z = ["12","60"]#["10","20","30","40","50","60"]  

        for i in range(len(wt_fa_water_to_bilayer)):
            for k in range(len(wt_f_elec_bilayer)):
                for l in range(len(wt_fa_imm_elec)):

                    for j in range(len(start_z)):

                        # the default lipid composition is DLPC, based on the test, the lipid composition and temperature needs to be changed.
                        if (pH_mode):
                                arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:script_vars interface=" + interface + " -parser::script_vars start_z=" + start_z[j] + " -parser::script_vars end_z=" + end_z[j] + " -parser::script_vars flag_axis=" + flag_axis[0] + "  -parser:protocol " + xml_script + " -mp:lipids:temperature 30.0"
                                arguments = arguments + \
                                    " -pH_mode true -value_pH " + str(pH)
                                
                        else:
                            print("pH mode is a must in titration")
                            sys.exit()    
                        if("ph-titration" in test_name): 
                            arguments = arguments + " -mp:lipids:composition DOPC"
                        else:
                            arguments = arguments + " -mp:lipids:composition POPC" 
                        casedir = outdir + "/fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                            wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                        if (not os.path.isdir(casedir)):
                            os.system("mkdir " + casedir)
                        os.chdir(casedir)

                        casedir = casedir + "/" + case
                        # if (pH_mode and ("ddG-of-pH-insertion" in test_name)):
                        #     casedir = casedir + "_" + str(pH)
                        
                        # if (pH_mode and ("adsorbed-peptide-tilt-angle" in test_name)):
                        #     casedir = casedir + "_" + str(pH)
                        if(pH_mode):
                            casedir = casedir + "_" + str(pH)
                            
                                
                        if (not os.path.isdir(casedir)):
                            os.system("mkdir " + casedir)
                        os.chdir(casedir)

                        if(energy_fxn=="franklin2021"):
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
                                    str(wt_fa_water_to_bilayer[i]) + " -set_weights f_elec_lipidlayer " + str(
                                        wt_f_elec_bilayer[i]) + " -set_weights fa_imm_elec " + str(wt_fa_imm_elec[i])
                        
                        # Write jobfile and submit to the HPC
                        print(
                            "Submitting orientation sampling calculations for energy landscape case:", case)
                        jobname = case + "_" + \
                            end_z[j] + "_" + start_z[j] + \
                            "_protein_energy_landscape"
                        print(jobname)
                        jobfile = hpc_util.make_jobfile(
                            casedir, jobname, executable, arguments)
                        #print( "outdir:~", outdir )
                        print("casedir:~", casedir)
                        print("case:~", case)
                        print("executable:~", executable)
                        print("arguments:~", arguments)

                        if(pH_mode and ("adsorbed-peptide-tilt-angle" in test_name) and pH == 14):
                            #hpc_util.submit_condor_job(
                            #    casedir, jobname, jobfile, "")
                            hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile )
                            # elif(("adsorbed-peptide-tilt-angle" not in test_name)):
                        else:
                            #hpc_util.submit_serial_rockfish_job(
                            #    casedir, jobname, jobfile)
                            hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile )
                            #hpc_util.submit_condor_job( casedir, jobname, jobfile, "" )
