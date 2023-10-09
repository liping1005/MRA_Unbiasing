#!/bin/bash -login
#SBATCH --job-name=sinc
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=12:00:00 --mem=50gb
#SBATCH --constraint=[intel16|intel18|amd20]
module load MATLAB/2022a
cd /mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical
matlab -nodisplay -r "sim_sinc_learn_eta"