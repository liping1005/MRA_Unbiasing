#!/bin/bash -login
#SBATCH --job-name=sim_f5_oracle
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00 --mem=50gb
module load MATLAB/2022a
cd /mnt/ufs18/home-109/yinlipi1/bispectrumV2/NumericalExperiments/Oracle
matlab -nodisplay -r "sim_f5"