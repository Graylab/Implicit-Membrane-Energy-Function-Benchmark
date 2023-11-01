#!/usr/bin/env python
""" relax structures from either ab-initio or crystal structures
Authors:
	Rituparna Samanta <rituparna@utexas.edu>
Example:

Arguments:


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

def run_abinitio_structures(nstruct, multithread=True):
    # Read path configuration file
    config = read_config.read_config()
    print("generating the sh file for relax")

    # Read list of energy landscape test cases
    targets_path = config.benchmark_path + "../../folding_sequence/disordered_states/"
    list_of_targets = targets_path + "testcases1.list"
    os.chdir(targets_path)
    
    with open(list_of_targets, 'rt') as f:
        testcases = f.readlines()
    testcases = [ x.strip() for x in testcases ]
    print(testcases)
    
    if(multithread):
      executable = config.rosetta_mpi_path +"/AbinitioRelax" + ".mpi." + config.platform + config.compiler + config.buildenv
    else:
      executable = config.rosetta_path + "/AbinitioRelax" + "." + config.platform + config.compiler + config.buildenv
    
    database = config.rosetta_path + "/../../database"

    for case in testcases:
    # Change directories to a data analysis dir
            outdir = targets_path + case
            outdir = outdir + "/output/"
            print("outdir:~",outdir)

            if ( not os.path.isdir( outdir ) ):
               os.system( "mkdir " + outdir )

            inputdir = targets_path + case

            pdb_input = inputdir + "/"+ case + "after_relax.pdb"
            print("pdbfile:~",pdb_input)

           # fasta_input = config.benchmark_path + "/targets/stability/C5_pHLIP_helical_peptides" + "/" + case + ".fasta"
            fasta_input = inputdir + "/" + case + ".fasta"
            input_3mer = inputdir + "/" +  case + ".robetta.3mer"
            input_9mer = inputdir + "/" + case + ".robetta.9mer"  

            argument = "-database " + database + " -in:file:native " + pdb_input + " -in:file:fasta " + fasta_input
            argument = argument + " -in:file:frag3 " + input_3mer + " -in:file:frag9 " + input_9mer 
            argument = argument + " -out:file:silent " + case + "_silent.out -out:path " + outdir 
            argument = argument + " -abinitio:recover_low_in_stages 1 2 3 4 -abinitio::increase_cycles 2 -overwrite " + "-nstruct " + str(nstruct)
          
            if(multithread):
                argument = argument+" -multiple_processes_writing_to_one_directory"

            os.chdir( inputdir )
            print("Submitting abinitio of case:", case )

            jobname = case + "_abinitio"
            print( jobname)

            jobfile = hpc_util.make_jobfile( inputdir, jobname, executable, argument )

            hpc_util.submit_mpi_rockfish_job( inputdir, jobname, jobfile, 2 )
            
            #all poses have to be converted to all_atom
            # can we do the last relax part in jupyter
            #Plot the score files. if they are mostly withing the range of 40-50
            #we got the same thing before as well.
            #or another thing we can select the best structure, do the last set of protocols 
            #in pyrosetta and match their backbone with the current set up. then relax to see. 

def run_relax_on_crystal(targets, targetlist = [], nstruct = 6, ramp_constraints=False):
    import hpc_util
    import read_config
    
    os.chdir('../../../tests/python/')
    # Read path configuration file
    config = read_config.read_config()
    print("generating the sh file for relax")

    # Read list of energy landscape test cases
    targets_path = config.benchmark_path + "targets/tilt_angle/" + targets + "/"
    if(targetlist==[]):
      list_of_targets = targets_path + targetlist
    else:
      list_of_targets = targets_path + "targets.list"
    os.chdir(targets_path)
    with open(list_of_targets, 'rt') as f:
        testcases = f.readlines()
    testcases = [ x.strip() for x in testcases ]
    print(testcases)
    executable = config.rosetta_path + "/relax" + "." + config.platform + config.compiler + config.buildenv
    database = config.rosetta_path + "/../../database"
    # For each test case, generate specific arguments, a condor file, and then run
    for case in testcases:
    # Change directories to a data analysis dir
            outdir = targets_path + case
            if ( not os.path.isdir( outdir ) ):
               os.system( "mkdir " + outdir )
    
            # inputdir = targets_path + case
            inputdir = targets_path 

            
            #pdb_input = inputdir + "/"+ case + "_ignorechain.pdb"
            pdb_input = inputdir + case + "/" + case +".pdb"
            #pdb_input = inputdir + case + "_renum.pdb"
            # pdb_input = inputdir + case + "_relaxed_001.pdb"

            print("pdbfile:~",pdb_input)

            if(ramp_constraints):
               argument = "-database " + database + " -relax:constrain_relax_to_start_coords " + "-in:file:s " + pdb_input
               argument = argument + " -ex1 -ex2 -ex1aro -overwrite " + "-nstruct " + str(nstruct) + " -out:path:pdb " + outdir + " -out:path:score " + outdir       
               argument = argument + " -relax:ramp_constraints true"
            else:
               #the vanilla code by AJ Vincelli: https://docs.google.com/document/d/1pekk-bYLQPX_qNVZjqppM-vheqFsw8Mwq8yDbI1U7GU/edit
               argument = "-database " + database + " -relax:constrain_relax_to_start_coords " + "-in:file:s " + pdb_input
               argument = argument + " -ex1 -ex2 -ex1aro -ex2aro -overwrite " + "-nstruct " + str(nstruct) + " -out:path:pdb " + outdir + " -out:path:score " + outdir
               argument = argument + " -relax:ramp_constraints false"
            
            os.chdir( inputdir )
            print("Submitting relax of case:", case )
    
            jobname = case + "_relaxedv3"
            print( jobname) 
            jobfile = hpc_util.make_jobfile( inputdir, jobname, executable, argument )
            print( "casedir:~", inputdir )
            print( "case:~", case )
            # print( "executable:~", executable )
            # print( "arguments:~", arguments )
            # hpc_util.submit_serial_rockfish_job( inputdir, jobname, jobfile )

def sort_scores(targets):
    # Read path configuration file
    config = read_config.read_config()
    print("generating the sh file for relax")

    # Read list of energy landscape test cases
    targets_path = config.benchmark_path + "targets/" + targets + "/"
    list_of_targets = targets_path + "targets1.list"
    os.chdir(targets_path)
    with open(list_of_targets, 'rt') as f:
        testcases = f.readlines()
    testcases = [ x.strip() for x in testcases ]
   #  print(testcases)
    
    # For each test case, generate specific arguments, a condor file, and then run
    for case in testcases:
    # Change directories to a data analysis dir
            outdir = targets_path + case
            # outdir = outdir + "/relaxed_structure_v2/"
            outdir = outdir + "_relaxed_structure_v3/"
            print("case:~",case)

            if ( not os.path.isdir( outdir ) ):
               os.system( "mkdir " + outdir )
    
            # inputdir = targets_path + case
            inputdir = targets_path 
            os.chdir( outdir )
            cmdline = "sort -n -k2 score.sc | awk " + "\'{print $2 \"\t\" $23}\'" + "> sorted_file"
            os.system(cmdline)
            sorted_score_file = outdir + 'sorted_file'
            with open(sorted_score_file, 'rt') as f:
               scores = f.readlines()
            scores = [ x.strip() for x in scores ]
            print(scores[0].split("\t")[:])
            print(len(scores[0].split("\t")))
            if(len(scores[0].split("\t"))==1):
               cmdline = "scp "+ outdir + scores[3].split("\t")[1] +".pdb " + outdir + "../" + case + "_relaxed_002.pdb"
               print(cmdline)
               print(scores[3].split("\t")[0])
               os.system(cmdline)
            elif(scores[0].split("\t")[1]=="description"):
               cmdline = "scp "+ outdir + scores[3].split("\t")[1] +".pdb " + outdir + "../" + case + "_relaxed_002.pdb"
               print(cmdline)
               print(scores[3].split("\t")[0])
               os.system(cmdline)
            else:   
               cmdline = "scp "+ outdir + scores[0].split("\t")[1] +".pdb " + outdir + "../" + case + "_relaxed_002.pdb"
               print(cmdline)
               print(scores[0].split("\t")[0])
               os.system(cmdline)
            # print(cmdline)
            

# def main( ):
#     # Read path configuration file
#     #config = read_config.read_config()
#     #targets = "stability/C5_pHLIP_helical_peptides/disordered_states"
#     #nstruct = 100 
#     #ramp_constraints = True
#     #run_relax_on_crystal(targets, nstruct, ramp_constraints)

#    #  targets = "stability/C5_pHLIP_helical_peptides"
#     targets = "tilt_angle/A6_AmphiScan_pdbs"
#     #targets = "orientation/B1_multispan_proteins"
#     nstruct = 500 
#    #  run_relax_on_crystal(targets, nstruct)
#     sort_scores(targets)

#     #run_abinitio_structures(50000,"TRUE")
#     #to analyze the structures:
#     #sort -n -k2 example_score_file | awk '{print $2 "\t" $3}'
# if __name__ == "__main__":
#     main()    

