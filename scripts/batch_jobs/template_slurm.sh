#!/bin/bash

#SBATCH -C haswell
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem 5000
#SBATCH --time=6:00:00
##SBATCH --mail-user=sri@physics.ucla.edu
##SBATCH --mail-type=ALL

##module load Python/2.7.9-GCC-4.9.2-bare
##module load xcb-proto/1.11-intel-2016.u3-Python-2.7.9
##module use /home/sri/modulefiles/
##module load anaconda
module load python

#python downloa_dr14_CasJobs.py 430000 500000
###python s2_fit_kappa_crossmaps_stacksample_Mlamb_0.2.py data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3_lgt_5/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/ sptpol_des 200 20 1000 0
### python s2_fit_kappa_crossmaps_stacksample_Mlamb_0.2.py data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3_lgt_5/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/ sptpol_des 200 20 40 0
##python s2_fit_kappa_crossmaps_stacksample_Mlamb_0.2.py data/sims/sptpol_des/crossmaps/sptpol_map_like_sims/data/year3_lgt_5/SNR_20.0_1000.0/0.5_no_tSZ_sptpol/ sptpol_des 200 20 30 0
