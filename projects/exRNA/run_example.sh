#!/bin/bash
if [ $# = 0 ]; then echo "Run as: ./run.sh #";exit;fi

proj="~/projects/exRNA/example"
hg38=~/genomes/human_hg38/sequence/GRCh38.p12.genome.fa
gtf=~/genomes/human_hg38/gtf
RNAs="miRNA
piRNA
Y_RNA_exon
snRNA
srpRNA_exon
tRNA
lncRNA_exon
mRNA_exon
"

#############################################################
# Setp 1.1: Index
#############################################################
# Step 1.1.1: gtf2fasta 
#############################################################
if [ "$1" = "1.1.1" ]; then
	for i in $RNAs ;do
 		echo "start $i :"
		bedtools getfasta -fi $hg38 -bed $gtf/"$i".gff -fo RNAs_fasta/"$i".fa
		echo "$i finished."
	done

#############################################################
# Topic: 1.1.2: fasta2bt3 
#############################################################
elif [ "$1" = "1.1.2" ]; then
	for i in $RNAs ; do
 		echo "start $i :"
		bowtie2-build RNAs_fasta/$i.fa RNAs_fasta/$i
		echo "$i finished."
	done

#############################################################
# Topic: 1.2: fastaqc 
#############################################################
elif [ "$1" = "1.2" ]; then
	for i in `ls example/fastq/*.fastq`  ; do
 		echo "start $i :"
		fastqc -q -o example/fastqc $i 
		echo "$i finished."
	done

#############################################################
# Topic: 1.3: cut adaptor 
#############################################################
elif [ "$1" = "1.3" ]; then
	for i in `ls example/fastq/*.fastq`  ; do
 		echo "start $i :"
		cutadapt -u -100 -q 30,30 --trim-n -m 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $i.cutadapt  $i 
		echo "$i finished."
	done



#end of the jobs	
###################################
else
	echo "$hg38"
	echo $hg38.fa 
	echo "Type in options"
#	awk -F '\t' '{if(NR>1 &&$5>0&&$6>0) print $0}' $i.txt > $i.filter.txt
#	cut -f 4-6 $i.filter.txt |sort -k1,1 -k2,2n -k3,3n > $i.filter.bed

fi



