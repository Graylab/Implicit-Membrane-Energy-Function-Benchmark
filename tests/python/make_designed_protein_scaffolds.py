#!/usr/bin/env python
""" Generate redesigned protein scaffolds

This script runs fixed backbone design on native membrane
protein structures, required data for the sequence recovery
and side chain distribution tests (Test 8 and 9). 

Authors: 
  Rebecca Alford <ralford3@jhu.edu> 

Example: 
  $ import make_designed_protein_scaffolds
  $ make_designed_protein_scaffolds.run_fixed_backbone_design_calc( 
    config, energy_fxn,  targets, test_name, restore )

Arguments: 
  - config: Container with path to benchmark and rosetta files
  - energy_fxn: Weights file for energy function of interest
  - targets: Location of design targets
  - test_name: Name of test
  - restore: Use talaris2014 energy function parameters

Requirements: 
  - Rosetta release 246 or greater
"""

import random, os, sys
import hpc_util, read_config
from string import Template

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )
from pyrosetta import * 
	
init()
  
def run_fixed_backbone_design_calc( config, energy_fxn, targets, test_name, restore, res_file_name=[],analysis = False ): 
  """
  A function to perform fixed backbone design on a set of targets
  
  Arguments: 
    energy_fxn = energy function to use for calculations (typically, name of the weights file)
    config = path to benchmark, Rosetta executables
    targets = list of targets
    xml = Path to RosettaScript defining the peptide landscape search protocol
  """

  print( "Performing fixed backbone design on set", targets ) 

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/design/monomer_chains_design.list"
  # list_of_targets = config.benchmark_path + "targets/design/monomer_chains_design_failed.list"
  with open( list_of_targets, 'rt' ) as f: 
    test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

  list_of_targets_highmem = config.benchmark_path + "targets/design/monomer_chains_design_failed.list"
  with open( list_of_targets_highmem, 'rt' ) as f: 
    test_cases_highmem = f.readlines()
    test_cases_highmem = [ x.strip() for x in test_cases_highmem ]
    
  
  # Generate path to executable
  executable = config.rosetta_path + "/fixbb" + "." + config.platform + config.compiler + config.buildenv

  # Make a directory for the subset and lipid composition
  #outdir = config.benchmark_path + "data/" + energy_fxn + "/sequence-recovery/highdGforDE" 
  #os.system( "mkdir " + outdir )
  #os.chdir( outdir )

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
  outdir = outdir + "/weights_from_test8_heavyatomcorrection_changedddG_disallowC"
  if ( not os.path.isdir( outdir ) ): 
      os.system( "mkdir " + outdir )
  os.chdir( outdir )
  
  wt_fa_water_to_bilayer = [0.863] 
  wt_fa_imm_elec = [0.001]
  wt_f_elec_bilayer = [0.152] 
  i=0
  
  outdir = outdir + "/fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                  wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                       
  if ( not os.path.isdir( outdir ) ): 
      os.system( "mkdir " + outdir )
  os.chdir( outdir )
  
  # resfileinfo = open('reslevelinfofile', 'w')
          
  # For each test case, generate specific arguments, condor files, and then run
  for case in test_cases:

    print(test_cases)
    # Make one directory per case
    casedir = outdir + "/" + case
    if ( not os.path.isdir( casedir ) ): 
      os.system( "mkdir " + casedir )
    
    print(casedir)
    
    os.chdir( casedir )

    # Setup arguments by substitution
    pdbfile = config.benchmark_path + "targets/design/" + case + "/" + case + "_tr_ignorechain.pdb"
    spanfile = config.benchmark_path + "targets/design/" + case + "/" + case + "_tr.span"
    s = Template( "-in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
    arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir )
    
    
    if(energy_fxn =="franklin2021"):
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
              wt_f_elec_bilayer[i]) + " -set_weights fa_imm_elec " + str(wt_fa_imm_elec[i]) + \
              " -set_weights menv_pH 0.0 "  
          #adding extra commands
                  
    else:
      arguments = arguments 
      #+ " -set_weights menv_pH 0.0 " 

    
    #default options
    #By default, fixbb will attempt to redesign all residues using all amino acids.
    #Interaction Graph (Default is to precompute all rotamer pair energies)
    #-use_input_sc: Include the side chain from the input pdb.  False by default.
    #Annealer: (Default is the original annealer used in [Kuhlman et al. 2003])
    #-score:weights: score12 is the default
    # analysis = True
    
    if(analysis):
        if(energy_fxn =="franklin2021"):  
          analysis_executable = "sort -n -k2 score.sc | head -n 100 | awk " + "'{print $25 " + '"\t"' + " $2}'" + " > score_rmsd"
          os.system(analysis_executable)
        elif(energy_fxn =="franklin2019"):
          analysis_executable = "sort -n -k2 score.sc | head -n 100 | awk " + "'{print $23 " + '"\t"' + " $2}'" + " > score_rmsd"
          os.system(analysis_executable)
        else:
          sys.exit('write the column which has description')  
    else:
      
        if(not res_file_name==[]):
          resfile = casedir + "/" + case + "_resfile"
          with open( resfile, 'w' ) as f: 
              f.write( "ALLAA \n" )
              f.write( "start \n" )
              
              native = pose_from_pdb(pdbfile)
              # resfileinfo.write("{} has {} res \n".format(case,native.total_residue()))
            
              for index in range(1,native.total_residue()): 
                  if(native.residue(index).name1() in ["C"]):
                      f.write( str(index) + " " + native.pdb_info().chain(index) + " NATAA" + "\n" )
              f.close()
          os.system( "chmod +x " + resfile )
            
        # res_file_name = config.benchmark_path + "tests/python/" + res_file_name
        arguments = arguments + " -resfile " + resfile
          
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
        else: 
          arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:temperature 37"

        arguments = arguments + " -nstruct 50 -ex1 -ex2 -ex1aro -ex2aro -pH_mode false -out:pdb_gz -run:multiple_processes_writing_to_one_directory"
        #based on the methods section of following paper 
        #A large scale test of computational protein 
        # design: folding and stability of nine completely redesigned globular proteins. Dantas G, et al. J Mol Biol. 2003 Sep 12;332(2):449-60
        # Write arguments and executable to a separate file
        jobfile = casedir + "/" + case + "_seqrecov.sh"
        with open( jobfile, 'w' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
        print("Submitting fixed backbone design calculation for sequence recovery case:", case) 
        queue_no = 1
        if(case in test_cases_highmem):
          high_mem = True
          print(case + 'is a high mem case')
        else:
          high_mem = False
          
        jobname = case+"_seqrec"
        if(high_mem):
          print('in highmem loop')
          hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile, num_tasks=24, high_mem=True ) 
        else:
          hpc_util.submit_serial_rockfish_job( casedir, jobname, jobfile) 
    #    hpc_util.submit_condor_job( casedir, case, executable, arguments, queue_no, high_mem )
  
  resfileinfo.close()
