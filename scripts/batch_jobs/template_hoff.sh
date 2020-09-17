#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o batch_jobs/joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=4:00:00,h_data=10G

#modules

# Your script content goes here...
