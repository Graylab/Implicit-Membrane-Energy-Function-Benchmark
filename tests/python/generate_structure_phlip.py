#!/usr/bin/env python

import sys, os, random
import hpc_util
import numpy as np
import read_config

def analyzing_dataset():
    executable = "sort -n -k2 score.sc | head -n 100 | awk " + "'{print $26" + '"\t" $2}' + " > score_rmsd.dat"
    return(executable)

def main():
    
    config = read_config.read_config()
    
    # Read list of energy landscape test cases
    list_of_targets = config.benchmark_path + "targets/stability/C5_pHLIP_helical_peptides/target1.list"
    with open( list_of_targets, 'rt' ) as f: 
        test_cases = f.readlines()
        test_cases = [ x.strip() for x in test_cases ]
    
    energy_fxn = "franklin2021"
    
    for name in test_cases:
        
        print("Submitting structure prediction for:", name )
        jobname = name + "_abinitio"
        print("jobname: "+ jobname)

        # Make output directories
        outdir = config.benchmark_path + "data/" 
        
        if ( not os.path.isdir( outdir ) ): 
            os.system( "mkdir " + outdir )
        os.chdir( outdir )
        
        outdir = outdir + energy_fxn 
        if ( not os.path.isdir( outdir ) ): 
            os.system( "mkdir " + outdir )
        os.chdir( outdir )

        outdir = outdir + "/abinitio_structure_of_pHLIP"
        if ( not os.path.isdir( outdir ) ): 
            os.system( "mkdir " + outdir )
        os.chdir( outdir ) 
        
        
        casedir = outdir + "/" + name
        if ( not os.path.isdir( casedir ) ): 
            os.system( "mkdir " + casedir ) 
        os.chdir( casedir )
        
        executable = "python3"
        arguments = config.benchmark_path + "tests/python/Predict_abinitio_structure.py -f " + energy_fxn +" -n "+ name\
            + " -p " + config.benchmark_path
        analysis = False
        if(analysis): 
            jobname = "analysis"  
            analysis_executable = analyzing_dataset()
            arguments = ""
            jobfile = hpc_util.make_jobfile( casedir, jobname, analysis_executable, arguments )
        else:
            jobfile = hpc_util.make_jobfile( casedir, jobname, executable, arguments )
            print( "casedir:~", casedir )
            hpc_util.submit_serial_rockfish_job_py36( casedir, jobname, jobfile )

if __name__ == "__main__": main()
                                               