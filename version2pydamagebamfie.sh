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

######fastq-dump the raw data we downloaded from ERRs via rentrez library######

#it should be noted that we used parallel fastq dump for speeding up the process. 
fastq_d_function <- function(ERR_vector_input){
  
  for(i in 1:length(ERR_vector_input))
  {
    #next line for debuging
    cat("current ERR input is -",ERR_vector_input[i])
    
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


# than trimmed and cleaned data with automate with a loop hoping to ease the process of spades. 
for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
trimmomatic PE -threads 4 ${base}_1.fastq ${base}_2.fastq \
${base}_1.trim.fastq ${base}_1.untrim.fastq \
${base}_2.trim.fastq ${base}_2.untrim.fastq \
SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done

#than I cat the 4 ERRs out of the 6 that the researches did becuase of imbalance in number alignments
# (ENA run accession codes ERR3678595, ERR3678598, ERR3678602, ERR3678603, 122 and ERR3678613)

# which through erros in the spades protocol them all together. 
cat files_project_28/*._1.trim.fastq >files_project_28/all_trim_reads.R1.fastq
cat files_project_28/*._2.trim.fastq >files_project_28/all_trim_reads.R2.fastq

# De novo assembly with metaspades (see our assembly 
# unit on using spades for help with this, metaspades is accessible in much the same way)
#I faced complication in the process of correcting mismatches via spades. but I did my best with what
# I had 


#ASSEMBLE spades####
#spades
module load spades
#gunzip files_project_28/ERR{somenumber}.fastq.gz

metaspades.py -t 10 -k 21,33,45 -1 files_project_28/all_trim_reads.R1.fastq-2 files_project_28/all_trim_reads.R2.fastq -o spades_output

#ran the metaspades as the above, but later did the following since we didnt recieve enough contigs for a valid bam file towards the pydamage applicatoin 
metaspades.py -t 8 -k 21,33,45 --careful -1 files_project_28/all_trim_reads.R1.fastq -2 files_project_28/all_trim_reads.R2.fastq -o spades_output

#I used the ability to return to the protocol based on the last checkpoint of spades thanks to the manuel
#because one of the debugging issues had been the process was extreme long and the sessions not as long. 
# do to a drop in threads, but we suspect that though we correct OMP number of threads, that the increase in the 
# threads in the pipline of spads going from 1 to 8 created a buttleneck effect or that the I/O disk interactions with 
# the data isnt sufficient in the cluster settings. hence we were unable to complete the mismatch corrections
# in the pipline these checkpoints can be used as a reference to where to return or continue from if there is an error
# or other issues.

export OMP_NUM_THREADS=32  
spades.py -t 16 --continue -o spades_output # now that we fixed the OMP from 1 to 8, attempt this again with -t 16

# we were able to visialize our contigs.fasta results to see if the contigs lengths match what we were striving for. 
# which is longer length of 1000 based on the paper. 
# given the failure to reach the end level of the spades and we unsure about the process the researches took, we are unable to 
# tell how and what they discovered and its methods/protocol.

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

# 
# so first we applied genome assembly protocol and concepts, and now we are moving on to the variant calling aspect of 
# the workflow to be followed by pydamage

# we created a data and results folders 
mkdir data
mkdir -p results_pyDamage/sam results_pyDamage/bam results_pyDamage/vcf

#another point of error solving one thought was that the fasta reference used in the documentation of the paper
# of Actinomyces dentalis, 72% GC was taken from a smaller sample of the reference sequence. and therefore had 
# no alignments for referenece to offer. it is unclear if that was an adiional cause of error, 
# however based on the plots recieved from the pydamage software later it does show different plot results possibly. 

# so we used these code lines to fetch teh fasta file and later bwa index it.  

esearch -db nucleotide -query "KE386864.1" | efetch -format fasta > KE386864.1.fasta




#index referene genome
bwa index KE386864.1.fasta

# but later repeated the process with the following reference 

GCF_000429225.1_ASM42922v1_genomic.fna

#index second time referene genome

bwa index GCF_000429225.1_ASM42922v1_genomic.fna

# align reads (reading )
bwa mem GCF_000429225.1_ASM42922v1_genomic.fna spades_output/K45/final_contigs.fasta> results_595_28/sam/K45_contigs.sam


#than we were working with SAM(BAM , BCF and VCF) files

#at first we used all_reads.sam and later repeated with the biggest K45 contigues we had from spades, K45_contigs.sam
#we had more data to use, but eventully that too was deemed to have faulties. 

module load samtools




#convert SAM to BAM

samtools view -S -b results_595_28/sam/K45_contigs.sam > results_595_28/bam/K45_contigs.bam



#sort bam file (by genomic coordinates)
samtools sort -o results_595_28/bam/K45_contigs.aligned.sorted.bam results_595_28/bam/K45_contigs.bam


samtools index results_595_28/bam/K45_contigs.aligned.sorted.bam


#check out our alignment
samtools flagstat results_595_28/bam/K45_contigs.aligned.sorted.bam 


./bam trimBam results_595_28/bam/K45_contigs.sam results_595_28/bam/after_trimBam.sam --clip 135
#this just shows that we are on the right tract 

#need to test this out. seen this on github possibly can assist and help 
#samtools calmd -b results_595_28/bam/K45_contigs.aligned.sorted.bam index GCF_000429225.1_ASM42922v1_genomic.fna > aln.bam

# we set up pydamge python based software in teh following manner :

cd ~/src


#copy the url for the git clone in the future. 
git clone https://github.com/maxibor/pydamage
cd pydamage
git checkout dev
conda env create -f environment.yml
conda activate pydamage
pip install -e .
#then test it
pydamage --help

#run 
# and at the end ran once we had the right formated bam file the next code -> pydamage analyze aligned.bam

#first time with results_595_28/bam/all_reads.aligned.sorted.bam second time with the following code
# with an aditional --plot flag for a constructed graphs and reualts. 
conda activate pydamage
pydamage analyze results_595_28/bam/K45_contigs.aligned.sorted.bam --plot 


# --plots
# The visual output are PNG files, one per reference contig. They show the frequency of observed C to T,
# and G to A transition at the 5' end of the sequencing data and overlay it with the fitted models for 
# both the null and the damage model, including 95% confidence intervals. Furthermore, it provides a 
# "residuals versus fitted" plot to help evaluate the fit of the pydamage damage model. Finally, the plot 
# contains informtion on the average coverage along the reference and the p-value calculated from the 
# likelihood-ratio test-statistic using a chi-squared distribution.

# observations:
# Coverage (or depth) in DNA sequencing is the number of unique reads that include a given nucleotide
# in the reconstructed sequence. Deep sequencing refers to the general concept of aiming for high
# number of unique reads of each region of a sequence. or in other words
# average number of reads that align to, or "cover," known reference bases. The sequencing 
# coverage level often determines whether variant discovery can be made with a certain degree 
# of confidence at particular base positions.ideal lower than 25% or higher than 75%
# 
# was at 0.2% in the example done in the paper


#differences could possibly be by the version or variation of the data via the id
# even though it was in march 

# taking apart the following 2 pictures. the NZ_AUBL01000094.1 and NZ_AUBL01000048.1 
# where in the first one we a p value of 0.004 which is great than 0.001 but is still less than 0.005 
# could argue to still have some sort of significance along side the RMSE value of residuals vs fitted at 0.009
# in comparison to the 0.01 in the example the paper conducted. 
When conducting a residual analysis, a "residuals versus fits plot" is the most frequently created plot. It is 
a scatter plot of residuals on the y axis and fitted values (estimated responses) on the x axis.
The plot is used to detect non-linearity, unequal error variances, and outliers.
the smaller that value the better. 


