#!/bin/bash -l

#----------------------------------------------------
# SLURM job script for membrane force field benchmarking applications
# Runs on Rockfish for serial applications
#----------------------------------------------------

#SBATCH --job-name=sc-distribution
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --time=10:0:0
#SBATCH --mem-per-cpu=12GB
#SBATCH --output /home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2019/sequence-recovery/weights_from_test8_heavyatomcorrection_changedddG_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001/sc-distribution.%j.out

#SBATCH --error /home/rsamant2/scratch16-jgray21/rsamant2/data/franklin2019/sequence-recovery/weights_from_test8_heavyatomcorrection_changedddG_disallowC/fa_wb_0.863_felecbilayer_0.152_fimm_0.001/sc-distribution.%j.err

conda deactivate
source /home/rsamant2/Softwares/py3.9_pyrosetta_gitcommit_7a030b9/venv_py3.9_pyR7a030b9/bin/activate

time
# python3 combiningfiles_peptide_titration.py --energy_fxn "franklin2021" --which_tests "ph-titration"
#python3 process_benchmark_data_diffwts.py --energy_fxn 'franklin2021' --which_tests 'sequence-recovery'

#python3 process_benchmark_data_diffwts.py --energy_fxn 'proteinmpnn' --which_tests 'sc-distribution'
#python3 process_benchmark_data_diffwts.py --energy_fxn 'franklin2021' --which_tests 'sc-distribution'
python3 process_benchmark_data_diffwts.py --energy_fxn 'franklin2019' --which_tests 'sc-distribution'
# python3 combiningfiles_peptide_titration.py --energy_fxn "franklin2021" --which_tests "tm-peptide-tilt-angle"
# python3 generate_benchmark_data.py --energy_fxn "franklin2021" --which_tests "sequence-recovery" 


time
