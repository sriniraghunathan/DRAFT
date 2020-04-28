#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o batch_jobs/joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=8:00:00,h_data=10G

#modules

# Your script content goes here...
##python get_pspec_galactic_sims.py -dust_or_sync dust -which_mask 2
python get_pspec_galactic_sims.py -dust_or_sync dust -which_mask 2 -nside 2048 -lmax 3500
