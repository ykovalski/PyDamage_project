#this is for the project 
#https://cran.r-project.org/web/packages/rentrez/rentrez.pdf


library(readr)
library(ShortRead)
library(Biostrings)
library(tidyverse)
library(tidyr)
library(rBLAST)
library(rentrez)
library(stringr)
library(ggplot2)

# fastq-dump the raw data
#setup####

fastq_d_function <- function(ERR_vector_input){
  
  for(i in 1:length(ERR_vector_input))
  {
    #next line for debuging
    cat("current SRR input is -",ERR_vector_input[i])
    
    fastq_d_system <-paste("parallel-fastq-dump --sra-id", ERR_vector_input[i], "--threads 8 --outdir ./ --split-files --gzip",sep=" ")
    system(fastq_d_system)
  }
  
}


fun<- function(){
  print("start fun()")
  
  xtrct_idNum_from_PRJNAinputRdLn <- paste(str_extract_all(readline("type in your PRJEB id number in here (i.e. PRJEB33577) "), "[0-9]", simplify = TRUE), collapse = "")
  
  PRJEB_number <- paste("PRJEB",xtrct_idNum_from_PRJNAinputRdLn,sep="")
  
  
  PRJEB_storage <- readline("type in directory name or desired path:\nExample 1: /projectnb2/ct-shbioinf/ykovalski/pyDamage_project/new_name_of_dir\nExample 2: new_name_of_dir")
  
  system(paste("mkdir",PRJEB_storage,sep = " "))
  
  setwd(PRJEB_storage)
  
  input_system_PRJEB_id <- paste("esearch -db bioproject -query",PRJEB_number,"| efetch -format runinfo -mode xml | xtract -pattern ProjectID -element CenterID", sep = " ") 
  id_insert <- system(input_system_PRJEB_id, intern = TRUE)
  get_ERR_files_names <- paste("esearch -db sra -query",id_insert,"| efetch -format runinfo -mode xml | xtract -pattern SraRunInfo -element Run", sep= " ")
  ERR_vector_preSplit <- system(get_ERR_files_names, intern = TRUE)
  ERR_vector_postSplit<- strsplit(ERR_vector_preSplit, '\t')
  
  fastq_d_function(ERR_vector_postSplit)
}

fun()

# converted into fastq
gunzip files_project_28/*.fastq.gz


# than trimmed and cleaned data with automate wiht a loop
for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
trimmomatic PE -threads 4 ${base}_1.fastq ${base}_2.fastq \
${base}_1.trim.fastq ${base}_1.untrim.fastq \
${base}_2.trim.fastq ${base}_2.untrim.fastq \
SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done

#than I cat them all together. 
cat files_project_28/*._1.trim.fastq >files_project_28/all_trim_reads.R1.fastq
cat files_project_28/*._2.trim.fastq >files_project_28/all_trim_reads.R2.fastq

# De novo assembly with metaspades (see our assembly 
# unit on using spades for help with this, metaspades is accessible in much the same way)

#the content is in a shell form to mimic the class room. 





#ASSEMBLE spades####

#spades
module load spades
#gunzip files_project_28/ERR{somenumber}.fastq.gz






# https://github.com/ablab/spades
module load spades

metaspades.py -t 10 -k 21,33,45 -1 files_project_28/ERR3678595_1.fastq -2 files_project_28/ERR3678595_2.fastq -o spades_output



library(ggplot2)
library(Biostrings)

Spades_assem_V <- function(){

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


#VARIANT CALLING proceedure.####
#Align to one of the bacterial genomes they write about in the paper with bwa as in the variant calling workflow.#


# https://www.ncbi.nlm.nih.gov/nuccore/AUBL00000000

#https://www.ncbi.nlm.nih.gov/nuccore/KE386864.1?report=fasta
#download fasta to a data folder
esearch -db nucleotide -query "KE386864.1" | efetch -format fasta > KE386864.1.fasta


mkdir -p results_595_28/sam results_595_28/bam results_595_28/vcf

#index the reference genome #if outside the data folder use data/KE386864.1.fasta
bwa index KE386864.1.fasta

# align reads (reading )
bwa mem KE386864.1.fasta spades_output/K45/final_contigs.fasta> results_595_28/sam/K45_contigs.sam

#head results_595_28/sam/ERR3678595.sam


#working with sam (bam, BCF and VCF) files

#working with SAM(BAM , BCF and VCF) files


#ERR=ERR3678595

module load samtools

#convert SAM to BAM

samtools view -S -b results_595_28/sam/K45_contigs.sam > results_595_28/bam/K45_contigs.bam

#sort bam file (by genomic coordinates)
samtools sort -o results_595_28/bam/K45_contigs.aligned.sorted.bam results_595_28/bam/K45_contigs.bam

samtools index results_595_28/bam/K45_contigs.aligned.sorted.bam
#check out our alignment
samtools flagstat results_595_28/bam/K45_contigs.aligned.sorted.bam 
#this just shows that we are on the right tract 


#Then:#####
# Convert to bam with samtools as in the variant calling workflow

# Send the bam file to pydamage###



#testing contigs out using Meghit
contigs = readDNAStringSet('dirMEGAHIT/final.contigs.fa', format='fasta')
print(contigs)

maxcontig=max(contigs@ranges@width)
maxcontigSTR = contigs[which(contigs@ranges@width == maxcontig)]
as.character(maxcontigSTR)

library(ggplot2)
widths = data.frame(contigs@ranges@width)
names(widths) = 'contig_length'
ggplot(widths) + 
  geom_histogram(aes(x=contig_length)) + 
  scale_x_log10()


pydamage analyze results_595_28/bam/K45_contigs.aligned.sorted.bam
