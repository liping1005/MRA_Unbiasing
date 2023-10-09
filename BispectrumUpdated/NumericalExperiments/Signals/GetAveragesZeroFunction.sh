#!/bin/bash -login
#SBATCH --job-name=ZeroFunction
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2 
#SBATCH --time=168:00:00  
#SBATCH --mem=100G
#SBATCH -A cmse
module load MATLAB/2022a
/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Signals/
matlab -nodisplay -r "GetAveragesZeroFunction"