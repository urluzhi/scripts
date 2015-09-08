#### The Yeast RNA-seq data were downloaded from GSE42983, random sample 1M reads per sample,
###  wild-type samples, two replicates (wt1.fq and wt2.fq)
### 0.5 hours H2O2 treatment (oxidative stress) samples, two replicates (wt1Ox.fq and wt2Ox.fq)

#!bin/bash
## A simple protocol for mapping reads to genome, assemble transcriptome and differentially expression analysis
## Created by Chao Di, 2013.11.1
##
####################   pre-process of raw data, sample 100k reads #############################
## download raw data from sra
#wget -bc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR636/SRR636633/SRR636633.sra # wt-1
#wget -bc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR636/SRR636634/SRR636634.sra # wt-2
#wget -bc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR636/SRR636636/SRR636636.sra # wt-1-OX
#wget -bc ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR636/SRR636637/SRR636637.sra # wt-2-OX

## convert to fq
#for i in `ls|grep "sra$"`;do echo $i; fastq-dump $i; done
## random select N reads

#for i in `ls ../../raw | grep fastq$`;do
#	j=${i%.*}
#	echo $j
#	# Starting FASTQ files
#	export FQ1=../../raw/$j.fastq;
#	# Output random reads file
#	export FQ1SUBSET=$j.rand.fq;
#	# How many random reads do we want?
#	export N=100000
#	# "linearize" the two mates into a single record.  Add a random number to the front of each line
#    cat $FQ1 | awk 'BEGIN{srand()}; {OFS="\t"; getline seqs; getline sep; getline quals; print rand(),$0,seqs,sep,quals}' | \
#	# sort by the random number
# 						 sort -k1,1 | \
#	# grab the first N records
#  						head -n $N | \
#	# Convert the stream back to FASTQ files.
#	 awk 'BEGIN{FS="\t";OFS="\n"}; {print $2,$3,$4,$5 >> ENVIRON["FQ1SUBSET"]}'
#done

## rename
#mv SRR636633.rand.fq wt1.fq
#mv SRR636634.rand.fq wt2.fq
#mv SRR636636.rand.fq wt1Ox.fq
#mv SRR636637.rand.fq wt2Ox.fq

