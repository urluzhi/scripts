#!/bin/bash
if [ $# = 0 ]; then echo "Run as: ./run.sh #";exit;fi

proj="$HOME/proj_exRNA/example"
hg38=$HOME/genomes/human_hg38/sequence/GRCh38.p10.genome.fa
gtf=$HOME/genomes/human_hg38/anno/gtf
index=$HOME/genomes/human_hg38/index/transcriptome_rsem_bowtie2
RNAs=( miRNA.gencode27 )


#############################################################
# Setp 0: Index Genome/Transcriptome
#############################################################
# 0 - RSEM index for transcritome 
# RSEM will take exons only
# only if you prepare gtf files in the right format
#############################################################
if [ "$1" = "0" ]; then
	for i in ${RNAs[@]} ; do
 		echo "start $i.gtf:"
		rsem-prepare-reference --gtf $gtf/$i.gtf --bowtie2 $hg38 RNA_index/$i
		echo "$i finished."
	done


#############################################################
# Topic: 1.1: fastaqc 
#############################################################
elif [ "$1" = "1.1" ]; then
	for i in `ls example/fastq/*.fastq`  ; do
 		echo "start $i :"
		fastqc -q -o example/fastqc $i 
		echo "$i finished."
	done

#############################################################
# Topic: 1.2: trim and cut adaptor 
#############################################################
elif [ "$1" = "1.2" ]; then
	for i in `ls example/fastq/*.fastq`  ; do
 		echo "start $i :"
		cutadapt -u -100 -q 30,30 --trim-n -m 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $i.trim  $i 
		echo "$i finished."
	done

#############################################################
# Topic: 1.3: fastaqc trimmed reads
#############################################################
elif [ "$1" = "1.3" ]; then
	for i in `ls example/fastq/*.fastq.trim`  ; do
 		echo "start $i :"
		fastqc -q -o example/fastqc $i 
		echo "$i finished."
	done

#############################################################
# Topic: 1.4: remove rRNA reads 
#############################################################
elif [ "$1" = "1.4" ]; then

	cd $proj/fastq 
	for i in `ls *.fastq`  ; do
		 		
		sample="${i/.fastq/}"
		sam=$sample.rRNA.sam
		bam=$sample.rRNA.bam
		fqin=$sample.fastq.trim
		fqout=$sample.no_rRNA.fq
		RNA=$index/rRNA
		echo "start $sample :"
		
		bowtie2 -p 4 --sensitive-local --norc --no-unal --un ../mapped/$fqout \
    		 -x $RNA $fqin -S ../mapped/$sam 
		#map to rRNAs (exons only, sense strand only)
		# (default) look for multiple alignments, report best, with MAPQ
		# --sensitive-local(default): allow no mismatch, etc
		# --norc:  do not align reverse-complement version of read
		# --no-unal: suppress SAM records for unaligned reads
		# --un: store unmapped reads
		# -x: indexed genome/transcriptome 
		# -S: output file foramt as sam

		samtools view -b ../mapped/$sam > ../mapped/$bam
		#convert sam to bam for mapped reads on rRNAs (sense strand)
	
		echo "$sample finished."
	done

#############################################################
# Topic: 1.5: map miRNAs first; 
# then do the same for other RNAs sequentially 
#############################################################
elif [ "$1" = "1.5" ]; then

	cd $proj/fastq 
	for i in `ls *.fastq`  ; do
		 		
		sample="${i/.fastq/}"
		fqin=$sample.no_rRNA.fq
		fqout=$sample.no_miRNA.fq
		sam=$sample.miRNA.sam
		bam=$sample.miRNA.bam
		RNA=$index/miRNA
		echo "start $sample :"
		
		bowtie2 -p 4 --sensitive-local --norc --no-unal --un ../mapped/$fqout \
    		 -x $RNA ../mapped/$fqin -S ../mapped/$sam 
		#map to miRNAs (pre-miRNA or pri-miRNA, sense strand only)
		# (default) look for multiple alignments, report best, with MAPQ
		# --sensitive-local(default): allow no mismatch, etc
		# --norc:  do not align reverse-complement version of read
		# --no-unal: suppress SAM records for unaligned reads
		# --un: store unmapped reads
		# -x: indexed genome/transcriptome 
		# -S: output file foramt as sam

		samtools view -b ../mapped/$sam > ../mapped/$bam
		#convert sam to bam for mapped reads on miRNAs (sense strand)
	
		echo "$sample finished."
	done



#############################################################
## Deprecated: 
# Topic: 1.4-old: remove rRNA reads 
#############################################################
elif [ "$1" = "1.4old" ]; then

	cd $proj/fastq 
	for i in `ls *.fastq`  ; do
		 		
		sample="${i/.fastq/}"
		sam=$sample.rRNA.sam
		bam=$sample.rRNA.bam
		fq=$sample.no_rRNA.fq
		echo "start $sample :"
		
		bowtie2 -p 4 --sensitive-local --no-unal --un ../mapped/$fq.raw \
    		 -x $index/rRNA $sample.fastq.trim -S ../mapped/$sam.raw 
	

		#map to rRNAs (exons only, but including antisense as well)
		# (default) look for multiple alignments, report best, with MAPQ
		# --sensitive-local(default): allow no mismatch, etc
		# --no-unal: suppress SAM records for unaligned reads
		# --un: store unmapped reads
		# -x: indexed genome/transcriptome 
		# -S: output file foramt as sam
	
		cd ../mapped

		samtools view -b -F 16 $sam.raw > $bam
		#clean reads mapped to the antisense of rRNAs
		#this is the final bam of mapped reads on rRNAs (sense strand)

		samtools view -b -f 16 $sam.raw > $bam.rev
		#collect reads mapped to the antisense of rRNAs
		samtools fastq $bam.rev > $bam.rev.fq
		#convert this bam to fastq
		cat $fq.raw $bam.rev.fq > $fq
		#merge the antisense reads to the unmapped raw reads
		
		cd ../fastq
		echo "$sample finished."
	done




#end of the jobs	
###################################
else
	echo "Type in options"
#	awk -F '\t' '{if(NR>1 &&$5>0&&$6>0) print $0}' $i.txt > $i.filter.txt
#	cut -f 4-6 $i.filter.txt |sort -k1,1 -k2,2n -k3,3n > $i.filter.bed

fi



