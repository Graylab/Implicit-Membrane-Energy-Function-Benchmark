#!/usr/bin/env python
###################################################################
#@file:         hpc_util.py                                                                                  
#@description:  Methods for submitting Rosetta jobs to SLURM and CONDOR clusters                                             
#@args:			jobname, executable, arguments, queue_no, high_mem                                                
#@author: 		Rebecca F. Alford                   
#modified by Rituparna Samanta (rituparna@utexas.edu)
#@email: 		rfalford12@gmail.com                                          
###################################################################

import random, sys, os

def make_jobfile( casedir, jobname, executable, arguments, suffix=".sh" ): 

    jobfile = casedir + "/" + jobname + suffix
    #the line below should be "w" not "a" or append mode. 
    with open( jobfile, 'w' ) as f: 
        f.write( "#!/bin/bash\n" )
        f.write( executable + " " + arguments + "\n" )
        f.close()
    os.system( "chmod +x " + jobfile )
    print(jobfile)
    return jobfile


def submit_condor_job( path, jobname, executable, arguments, queue_no=1, high_mem=False ):

    # Given a test path and
    filename = path + "/" + jobname + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + jobname + ".out\n" )
        f.write( "error = " + path + "/" + jobname + ".err\n" )
        if ( high_mem == True ): 
            f.write( "request_memory = 15000\n")
        f.write( "notify_user = rsamanta89@gmail.com\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    # os.system( "condor_submit " + filename )

def submit_stampede_job( path, jobname, jobfile, num_nodes=1 ): 
  
    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Stampede2 with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH -J " + jobname + "\n" )
        f.write( "#SBATCH -p normal\n" )
        f.write( "#SBATCH -N " + str(num_nodes) + "\n" )
        f.write( "#SBATCH -n 16\n" )
        f.write( "#SBATCH -t 24:0:0\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH -o " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH -e " + path + "/" + jobname + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "#SBATCH -A TG-MCB180056\n" )
        f.write( "\n" )

        # Specify required modules
        f.write( "module load intel/18.0.0\n" )
        f.write( "export MKL_MIC_ENABLE=1\n" )

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpiexec bash " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    #sbatch_command = "sbatch " + filename
    #os.system( sbatch_command )


def submit_marcc_job( path, jobname, jobfile, num_nodes=1 ): 

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on MARCC with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + jobname + "\n" )
        f.write( "#SBATCH --partition=parallel\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=60:0:0\n" )
        f.write( "#SBATCH --mem=120GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + jobname + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        f.write( "module unload openmpi gcc\n" )
        f.write( "module load intel-mpi git\n" )
        f.write( "ml --gcc\n")
        f.write( "ml gcc/4.8.2\n")

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpirun ./" + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    #sbatch_command = "sbatch " + filename
    #os.system( sbatch_command )

def submit_serial_rockfish_job( path, jobname, jobfile, num_nodes=1, num_tasks=1, high_mem=False):

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch"
    with open( filename, 'w' ) as f:

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Rockfish for serial applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + jobname + "\n" )
        f.write( "#SBATCH --partition=defq\n" )
        f.write( "#SBATCH --account=jgray21\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        
        if(high_mem):
            f.write( "#SBATCH --ntasks-per-node=" + str(num_tasks) + "\n" )
        f.write( "#SBATCH --time=2:0:0\n" )
        # if(high_mem):
        #     f.write( "#SBATCH --mem-per-cpu=15GB\n" )
        # else:
        #     f.write( "#SBATCH --mem-per-cpu=2GB\n" )
        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + jobname + ".%j.err\n" )
        #f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        #f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        #f.write( "module unload openmpi gcc\n" )
        #f.write( "module load intel-mpi git\n" )
        f.write( "ml gcc\n")
        #f.write( "ml gcc/4.8.2\n")

        # Provide a description of the job
        f.write(  "echo Starting serial job" + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "srun " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    # os.system( sbatch_command )

def submit_serial_rockfish_job_py36( path, jobname, jobfile, num_nodes=1 ):

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch"
    with open( filename, 'w' ) as f:

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Rockfish for serial applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + jobname + "\n" )
        f.write( "#SBATCH --partition=defq\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=36:0:0\n" )
        f.write( "#SBATCH --mem-per-cpu=2GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + jobname + ".%j.err\n" )
        #f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        #f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        #f.write( "module unload openmpi gcc\n" )
        #f.write( "module load intel-mpi git\n" )
        f.write( "source /home/rsamant2/scratch16-jgray21/rsamant2/venv_py3.6/bin/activate\n" )
       

        # Provide a description of the job
        f.write(  "echo Starting serial job" + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "srun " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    os.system( sbatch_command )


def submit_array_rockfish_job( outdir, case, pH_value, num_nodes=1 ):

    # Create a new sbatch file named for the job type and test
    casedir = outdir + "/" + case + "_" + str(pH_value)
    filename = casedir + "/" + case + ".sbatch"
    with open( filename, 'w' ) as f:

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Rockfish for serial applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )
        

        # Write the job information
        f.write( "#SBATCH --job-name=" + case + "\n" )
        f.write( "#SBATCH --partition=defq\n" )
        f.write( "#SBATCH --account=jgray21\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --ntasks-per-node=" + str(1) + "-" + str(12) + "\n" )
        f.write( "#SBATCH --time=10:0:0\n" )
        f.write( "#SBATCH --mem-per-cpu=2GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + casedir + "/" + case + "_%a.out\n" )
        f.write( "#SBATCH --error " + casedir + "/" + case + "_%a.err\n" )
        #f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        #f.write( "module unload openmpi gcc\n" )
        #f.write( "module load intel-mpi git\n" )
        f.write( "ml gcc\n")
        #f.write( "ml gcc/4.8.2\n")

        # Provide a description of the job
        f.write(  "echo Starting serial array job" + case + "\n" )

        # Run the job
        f.write( "time\n" )
        start_z = ["-60","-50","-20","-10","0","-40","-30","10","20"]
        end_z = ["-50","-40","-10","0","10","-30","-20","20","30"]
        #pH_value = ["4","8"]

        
        for j in range(len(start_z)):
            casedir = outdir + "/" + case + "_" + str(pH_value)
            jobname = case + "_" + end_z[j] + "_" + start_z[j] + "_protein_energy_landscape"
            jobfile = casedir + "/" + jobname + ".sh"
            # Write job specific output and reporting information
            #f.write( "#SBATCH --output " + casedir + "/" + jobname + ".out\n" )
            #f.write( "#SBATCH --error " + casedir + "/" + jobname + ".err\n" )

            f.write( "srun " + jobfile + " &" + "\n" )
        #f.write( "time\n" )
        f.write("wait")
        f.close()

    # Run the slurm job file
    #sbatch_command = "sbatch " + filename
    #os.system( sbatch_command )

def submit_mpi_rockfish_job( path, jobname, jobfile, num_nodes=1 ):
    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + "_mpi_"+".sbatch"
    with open( filename, 'w' ) as f:

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Rockfish for serial applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + jobname + "\n" )
        f.write( "#SBATCH --account=jgray21\n" )
        f.write( "#SBATCH --partition=defq \n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=6:0:0 \n" )
        f.write( "#SBATCH --mem-per-cpu=2GB \n" )
        f.write( "#SBATCH --ntasks-per-node=" + str(24) + "\n" )
 
        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + jobname + ".%j.err\n" )
        #f.write( "#SBATCH --mail-user=rsamanta89@gmail.com\n" )
        #f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )
    
        #f.write("ml purge \n")
        f.write("ml gcc \n")
        #f.write("ml intel intel-mpi intel-mkl \n")
        #f.write("ml python \n")
        # ml gcc openmpi 
        # ml python

        # Provide a description of the job
        f.write(  "echo Starting mpi job " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpirun -np 24 " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()
    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    os.system( sbatch_command )
