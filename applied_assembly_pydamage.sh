#this is a shell script for assembly towards the pydamage project

# http://bio331.devbioinformatics.org/assembly.html
conda deactivate 
mkdir genome_assembly 
cd genome_assembly
cd /projectnb2/ct-shbioinf/ykovalski/pyDamage_project/data


#used selected raw data in the files used for project 28 folder

#ASSEMBLE##

#spades
module load spades
#gunzip files_project_28/ERR{somenumber}.fastq.gz

gunzip files_project_28/*.fastq.gz
gunzip files_project_28/ERR3678595_2.fastq.gz

# https://github.com/ablab/spades
metaspades.py -t 8 -k 21,33,55,77 -1 files_project_28/ERR3678595_1.fastq -2 files_project_28/ERR3678595_2.fastq -o spades_output

metaspades.py -t 8 -k 21,33,45 --careful -1 files_project_28/all_trim_reads.R1.fastq -2 files_project_28/all_trim_reads.R2.fastq -o spades_output






#!/usr/bin/env Rscript
Spades_assem_V <- function(){
  library(ggplot2)
  library(Biostrings)
  contigs = readDNAStringSet('spades_output/contigs.fasta', format='fasta')
  #print contigs output presetns the DNAStringsSet object and the length.
print(contigs)

maxcontig=max(contigs@ranges@width)
maxcontigSTR = contigs[which(contigs@ranges@width == maxcontig)]
#this would show the node length and the sequence. 
as.character(maxcontigSTR)

#the following is to visualize the contig_length vs counts in a histogram plot
widths = data.frame(contigs@ranges@width)
names(widths) = 'contig_length'
ggplot(widths) + 
  geom_histogram(aes(x=contig_length)) + 
  scale_x_log10()
}

Spades_assem_V()


#VARIANT CALLING proceedure.


cd data

mkdir -p results_pyDamage/sam results_pyDamage/bam results_pyDamage/vcf

#index the reference genome
bwa index /ecoli_rel606.fasta

# align reads (reading )
bwa mem data_YK/ecoli_rel606.fasta \
data_YK/SRR2584863_1.trim.fastq data_YK/SRR2584863_2.trim.fastq \
> results_YK/sam/SRR2584863.sam

head results_YK/sam/SRR2584863.sam






# --------

export OMP_NUM_THREADS=32  
spades.py -t 16 --continue -o spades_output # now that we fixed the OMP from 1 to 8, attempt this again with -t 16 