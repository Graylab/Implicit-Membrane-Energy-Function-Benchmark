#!/bin/bash -l

#----------------------------------------------------
# SLURM job script for membrane force field benchmarking applications
# Runs on Rockfish for serial applications
#----------------------------------------------------

#SBATCH --job-name=1fep_relaxedv2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --time=2:0:0
#SBATCH --mem-per-cpu=2GB
#SBATCH --output /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/1fep/1fep_relaxedv2.%j.out
#SBATCH --error /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/1fep/1fep_relaxedv2.%j.err

ml gcc
echo Starting serial job/home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/1fep/1fep_relaxedv2.sh
time
srun /home/rsamant2/Softwares/Implicit-Membrane-Energy-Function-Benchmark-Electrostatics/Implicit-Membrane-Energy-Function-Benchmark-master/targets/orientation/B1_multispan_proteins/1fep/1fep_relaxedv2.sh
time
