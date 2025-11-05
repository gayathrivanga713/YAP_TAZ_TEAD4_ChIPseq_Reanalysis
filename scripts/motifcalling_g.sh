#!/bin/bash

#SBATCH --chdir=/home/g/gayathripriyav/midterm/
#SBATCH --job-name=MotifCalling
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH -t 02:00:00
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --mail-user=gayathripriyav@usf.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bash_profile
module purge
conda activate homer

findMotifsGenome.pl /home/g/gayathripriyav/midterm/YAP_peak/YAP_peaks.narrowPeak hg38 /home/g/gayathripriyav/midterm/YAP_peak/motif_results/ -size 100 -mask
findMotifsGenome.pl /home/g/gayathripriyav/midterm/TAZ_peak/TAZ_peaks.narrowPeak hg38 /home/g/gayathripriyav/midterm/TAZ_peak/motif_results_TAZ/ -size 100 -mask
findMotifsGenome.pl /home/g/gayathripriyav/midterm/TEAD4_peak/TEAD4_peaks.narrowPeak hg38 /home/g/gayathripriyav/midterm/TEAD4_peak/motif_results_TEAD4/ -size 100 -mask