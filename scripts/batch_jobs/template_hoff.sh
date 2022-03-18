#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o batch_jobs/joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=6:00:00,h_data=6G

#modules

# Your script content goes here...
