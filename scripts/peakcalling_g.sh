#!/bin/bash

#BATCH --chdir=/home/g/gayathripriyav/midterm
#SBATCH --job-name=MACS2_PeakCalling
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --mem=10000
#SBATCH -t 48:00:00
#SBATCH -o output=macs2.out
#SBATCH -o output=macs2.err
#SBATCH -e run.err
#SBATCH --mail-user=gayathripriyav@usf.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bash_profile
module purge
conda activate macs3

macs3 callpeak -t SRR1810900.bam -c SRR1810912.bam -f BAM -n YAP -g hs --outdir ~/midterm/YAP_peak
macs3 callpeak -t SRR1810907.bam -c SRR1810912.bam -f BAM -n TAZ -g hs --outdir ~/midterm/TAZ_peak
macs3 callpeak -t SRR1810918.bam -c SRR1810912.bam -f BAM -n TEAD4 -g hs --outdir ~/midterm/TEAD4_peak
