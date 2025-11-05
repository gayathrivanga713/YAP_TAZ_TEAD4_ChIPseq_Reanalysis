#!/bin/bash

#BATCH --chdir=/home/v/gayathripriyav
#SBATCH --job-name=qualityreads
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=2
#SBATCH --mem=10000
#SBATCH -t 03:00:00
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --mail-user=gayathripriyav@usf.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

conda create -n mid_term
conda activate mid_term

mkdir -p midterm
cd midterm

#wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/000/SRR1810900/SRR1810900.fastq.gz"
#wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/007/SRR1810907/SRR1810907.fastq.gz"
#wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/008/SRR1810918/SRR1810918.fastq.gz"
#wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR181/002/SRR1810912/SRR1810912.fastq.gz"

module purge
module add apps/fastqc/0.11.5

fastqc ~/midterm/SRR1810900.fastq.gz
fastqc ~/midterm/SRR1810907.fastq.gz
fastqc ~/midterm/SRR1810918.fastq.gz
fastqc ~/midterm/SRR1810912.fastq.gz
fastqc ~/midterm/h3k4me1.fastq.gz
fastqc ~/midterm/h3k4me3.fastq.gz
fastqc ~/midterm/h3k27ac.fastq.gz
