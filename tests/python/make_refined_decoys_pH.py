#!/usr/bin/env python
###################################################################
#@file:         make_refined_models.py                                  
#@description:  Generate relaxed candidate structures                                                       
#@author:     	Rebecca F. Alford                   
#@email:    		rfalford12@gmail.com                                          
###################################################################

## TODO: Rename sampled_candidates files

import random, os, sys
import hpc_util, read_config
from string import Template

def run_refine_decoys_calc( energy_fxn, restore, config, test_name, resolution, xml ): 
    """
    A function for refining canddiate structures for the decoy discrimination test

    Arguments: 
        energy_fxn = energy function to use for calculations (typically, name of the weights file)
        config = path to benchmark, Rosetta executables
        targets = list of targets
        test_name = Name of test
    """

    print( "Generating data for decoy discrimination test...")

    # Read list of test case IDs and PDBs
    targets_path = config.benchmark_path + "targets/structure/D6_decoy_discrimination/" + resolution + "/"
    #targets_path = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/native_structure_recovery/'
    list_of_targets = targets_path + "targets.list"

    with open( list_of_targets, 'rt' ) as f: 
        targets = f.readlines()
        targets = [ x.strip() for x in targets ]

    # Generate path to executable
    #executable = config.rosetta_path + "/rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
    executable = config.rosetta_path + "/relax" + "." + config.platform + config.compiler + config.buildenv
    xml_script =  config.benchmark_path + "tests/xml/" + xml

    # Change directories to a data analysis dir
    #outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name + "/highdGforDE/"   
    outdir = '/home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/native_structure_recovery/'
    #outdir = targets_path
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )

    

    # Iterate through each target
    for target in targets:

        print("Submitting refinement calculations for case:", target)

        # Read the list of decoy lists
        list_of_decoy_lists = []
        if ( resolution == "highres" ):
        	#list_of_decoy_lists = targets_path + target + "/list.of.decoy.lists"
            list_of_decoys = targets_path + target + "/decoy_test.list"
        else: 
        	list_of_decoy_lists = targets_path + target + "/" + target + "_model_subset_lists"
        #list_of_decoy_lists = targets_path + target + "/output/" + "decoy.list"
        #list_of_decoys = targets_path + target + "/output/" + "br_structure.list"
        
        with open( list_of_decoys, 'rt' ) as f: 
            list_of_lists = f.readlines()
            list_of_lists = [ x.strip() for x in list_of_lists ]

        # Read general target variables
        spanfile = targets_path + target + "/" + target + ".span" 
        scorefile = target + "_refined.sc"
        casedir = outdir + target
        if ( not os.path.isdir( casedir ) ): 
            os.system( "mkdir " + casedir )
        os.chdir( casedir )

        # Iterate through each decoy list
        i = 0
        for dlist in list_of_lists: 

            i = i + 1

            # Setup case-specific variables
            #models_list = targets_path + target + "/output/" + dlist
            structure = dlist
            # Generate a string of arguments from the case-specific variables
            #s = Template( " -relax:constrain_relax_to_start_coords -in:file:l $models_list -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct $nmodels -run:multiple_processes_writing_to_one_directory ")
            s = Template( "-ignore_unrecognized_res -relax:constrain_relax_to_start_coords -in:file:s $structure -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct $nmodels -run:multiple_processes_writing_to_one_directory")
            #arguments = s.substitute( models_list=models_list, span=spanfile, xml=xml_script, sfxn=energy_fxn, outdir=casedir, nmodels=5, scorefile=scorefile )
            arguments = s.substitute( structure=structure, span=spanfile, xml=xml_script, sfxn=energy_fxn, outdir=casedir, nmodels=5, scorefile=scorefile )

            # Restore/lipid composition flags
            if ( restore ): 
                arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
            else:
                arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:has_pore false"
            
            if( energy_fxn=="franklin2021"):
                wt_fa_water_to_bilayer = [1.629]#, 1.484] 
                wt_fa_imm_elec = [-1.561]#, 0.571]
                wt_f_elec_bilayer = [0.136]#, -0.462]
                if( wt_fa_imm_elec[0]<0 or wt_f_elec_bilayer[0]<0 ):
                        filename = 'flag_file'
                        #the line below should be "w" not "a" or append mode. 
                        with open( filename, 'w' ) as f: 
                            f.write(  "-set_weights fa_water_to_bilayer "+ str(wt_fa_water_to_bilayer[0]) + "\n" )
                            f.write(  "-set_weights f_elec_lipidlayer "+ str(wt_f_elec_bilayer[0]) + "\n" )
                            f.write(  "-set_weights fa_imm_elec "+ str(wt_fa_imm_elec[0]) + "\n" )

                        f.close()
                        os.system( "chmod +x " + filename )
                        arguments = arguments + " @flag_file"
                
            if(target == "7sqg"):
                arguments = arguments + " -pH_mode true -value_pH 8.0"
            elif( target == "7sqh"):
                arguments = arguments + " -pH_mode true -value_pH 4.50"
            elif( target == "7sqf"):
                arguments = arguments + " -pH_mode true -value_pH 4.50"
            else:
                arguments = arguments + " -pH_mode true -value_pH 7.0"
            
            # Write jobfile and submit to the HPC
            jobname = target + "_refine_" + str(i)
            jobfile = hpc_util.make_jobfile(
                            casedir, jobname, executable, arguments)

            #hpc_util.submit_condor_job( casedir, jobname, executable, arguments, 5 )
            hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile )
                        


def run_refine_kined_structures( energy_fxn, restore, config, test_name, xml ): 
    """
    A function for refining canddiate structures for the decoy discrimination test

    Arguments: 
        energy_fxn = energy function to use for calculations (typically, name of the weights file)
        config = path to benchmark, Rosetta executables
        targets = list of targets
        test_name = Name of test
    """

    print( "Generating data for decoy discrimination test...")

    # Read list of test case IDs and PDBs
    targets_path = config.benchmark_path + "targets/structure/D5_helix_kinks/"
    list_of_targets = targets_path + "targets.list"

    with open( list_of_targets, 'rt' ) as f: 
        targets = f.readlines()
        targets = [ x.strip() for x in targets ]

    # Generate path to executable
    executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
    xml_script =  config.benchmark_path + "tests/xml/" + xml

    # Change directories to a data analysis dir
    outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name 
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Iterate through each target
    for target in targets:

        print("Submitting refinement calculations for helix kink case:", target)

        # Read the list of decoy lists
        decoy_list = targets_path + "/" + target + "/models.tr.clean.list"
        with open ( decoy_list ) as f: 
            pdbs = f.readlines()
            pdbs = [ x.strip() for x in pdbs ]

        i = 0
        for pdb in pdbs: 

            i = i + 1

            # Read general target variables
            spanfile = targets_path + "/" + target + "/" + target + ".span" 
            scorefile = target + "_refined_" + str(i) + ".sc"
            casedir = outdir + "/" + target
            if ( not os.path.isdir( casedir ) ): 
                os.system( "mkdir " + casedir )
                os.chdir( casedir )

            # Generate a string of arguments from the case-specific variables
            s = Template( " -relax:constrain_relax_to_start_coords -in:file:s $models_list -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct $nmodels -run:multiple_processes_writing_to_one_directory ")
            arguments = s.substitute( models_list=pdb, span=spanfile, xml=xml_script, sfxn=energy_fxn, outdir=casedir, nmodels=5, scorefile=scorefile )

            # Restore/lipid composition flags
            if ( restore ): 
                arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior "
            else:
                arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:has_pore false "

            # Write jobfile and submit to the HPC
            jobname = target + "_refine_" + str(i)
            hpc_util.submit_condor_job( casedir, jobname, executable, arguments, 5 )

