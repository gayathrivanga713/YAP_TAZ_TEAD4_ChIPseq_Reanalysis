#!/bin/bash

#SBATCH --chdir=/home/g/gayathripriyav/midterm
#SBATCH --job-name=ChIPseqMapping
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --mail-user=gayathripriyav@usf.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module purge
module add apps/bowtie/2.3.2
module add apps/samtools/1.3.1

bowtie2 -x ~/midterm/refgenome/GRCh38_noalt_as -U ~/midterm/SRR1810900.fastq.gz -S ~/midterm/SRR1810900.sam --threads 6 -k 1 2> ~/midterm/SRR1810900.out
samtools view -u ~/midterm/SRR1810900.sam | samtools sort -o ~/midterm/SRR1810900.bam
samtools index ~/midterm/SRR1810900.bam

bowtie2 -x ~/midterm/refgenome/GRCh38_noalt_as -U ~/midterm/SRR1810907.fastq.gz -S ~/midterm/SRR1810907.sam --threads 6 -k 1 2> ~/midterm/SRR1810907.out 
samtools view -u ~/midterm/SRR1810907.sam | samtools sort -o ~/midterm/SRR1810907.bam
samtools index ~/midterm/SRR1810907.bam

bowtie2 -x ~/midterm/refgenome/GRCh38_noalt_as -U ~/midterm/SRR1810912.fastq.gz -S ~/midterm/SRR1810912.sam --threads 6 -k 1 2> ~/midterm/SRR1810912.out
samtools view -u ~/midterm/SRR1810912.sam | samtools sort -o ~/midterm/SRR1810912.bam
samtools index ~/midterm/SRR1810912.bam

bowtie2 -x ~/midterm/refgenome/GRCh38_noalt_as -U ~/midterm/SRR1810918.fastq.gz -S ~/midterm/SRR1810918.sam --threads 6 -k 1 2> ~/midterm/SRR1810918.out
samtools view -u ~/midterm/SRR1810918.sam | samtools sort -o ~/midterm/SRR1810918.bam
samtools index ~/midterm/SRR1810918.bam

bowtie2 -x new_refgenome/GRCh38_noalt_as -U ~/midterm/h3k27ac_trimmed.fastq.gz -S ~/midterm/h3k27ac_trimmed.sam --threads 6 -k 1 -q 2> ~/midterm/h3k27ac_trimmed.out
samtools view -u ~/midterm/h3k27ac_trimmed.sam | samtools sort -o ~/midterm/h3k27ac_trimmed.bam
samtools index ~/midterm/h3k27ac_trimmed.bam

bowtie2 -x new_refgenome/GRCh38_noalt_as -U ~/midterm/h3k4me3_cleaned_trimmed.fastq.gz -S ~/midterm/h3k4me3_cleaned_trimmed.sam --threads 6 -k 1 -q 2> ~/midterm/h3k4me3_cleaned_trimmed.out
samtools view -u ~/midterm/h3k4me3_cleaned_trimmed.sam | samtools sort -o ~/midterm/h3k4me3_cleaned_trimmed.bam
samtools index ~/midterm/h3k4me3_cleaned_trimmed.bam

bowtie2 -x new_refgenome/GRCh38_noalt_as -U ~/midterm/h3k4me1.fastq.gz -S ~/midterm/h3k4me1.sam --threads 6 -k 1 -q 2> ~/midterm/h3k4me1.out
samtools view -u ~/midterm/h3k4me1.sam | samtools sort -o ~/midterm/h3k4me1.bam
samtools index ~/midterm/h3k4me1.bam

