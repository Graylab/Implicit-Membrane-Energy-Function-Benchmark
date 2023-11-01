#!/bin/bash -l

#----------------------------------------------------
# SLURM job script for membrane force field benchmarking applications
# Runs on Rockfish for serial applications
#----------------------------------------------------

#SBATCH --job-name=ph-titration
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --time=10:0:0
#SBATCH --mem-per-cpu=12GB
#SBATCH --output /home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/ph-titration/weights_from_test7_titration/fa_wb_0.966_felecbilayer_0.016_fimm_0.379/titration.%j.out
#SBATCH --error /home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2021/ph-titration/weights_from_test7_titration/fa_wb_0.966_felecbilayer_0.016_fimm_0.379/titration.%j.err
conda deactivate
source /home/rsamant2/Softwares/venv/bin/activate

time
# python3 combiningfiles_peptide_titration.py --energy_fxn "franklin2021" --which_tests "ph-titration"
python3 combiningfiles_peptide_titration.py --energy_fxn "franklin2021" --which_tests "adsorbed-peptide-tilt-angle"
python3 combiningfiles_peptide_titration.py --energy_fxn "franklin2021" --which_tests "tm-peptide-tilt-angle"
time
