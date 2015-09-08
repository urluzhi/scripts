#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe make 4

### The following steps are adapted from (Trapnell et al., Nature Protocols, 2012)
cd /data/Protocol/diff_exp
source bin_path 

#Part I: Align the RNA-seq reads to the genome
#1. Assemble transcripts for each sample:
tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o wt1_thout bowtie_index/YeastGenome Raw_reads_100k/wt1.fq 
tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o wt2_thout bowtie_index/YeastGenome Raw_reads_100k/wt2.fq 
tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o wt1X_thout bowtie_index/YeastGenome Raw_reads_100k/wt1X.fq 
tophat -p 4 -G yeast_annotation.gff --no-coverage-search -o wt2X_thout  bowtie_index/YeastGenome Raw_reads_100k/wt2X.fq 

echo "=============================================="
echo "Tophat mapping and cufflinks assembling done ! "
echo "=============================================="

#2. Extract mapped reads on chr I  only (for quick running)
#index bam file first
bamtools index -in wt1_thout/accepted_hits.bam 
bamtools index -in wt2_thout/accepted_hits.bam 
bamtools index -in wt1X_thout/accepted_hits.bam
bamtools index -in wt2X_thout/accepted_hits.bam 
#extract
bamtools filter -region chrI -in wt1_thout/accepted_hits.bam -out wt1_thout/chrI.bam
bamtools filter -region chrI -in wt2_thout/accepted_hits.bam -out wt2_thout/chrI.bam
bamtools filter -region chrI -in wt1X_thout/accepted_hits.bam -out wt1X_thout/chrI.bam
bamtools filter -region chrI -in wt2X_thout/accepted_hits.bam -out wt2X_thout/chrI.bam

echo "=============================================="
echo "bamtools filter is done ! "
echo "=============================================="

#Part II: Assemble transcripts on Chrom I by Cufflinks 
#1. Assemble transcripts for each sample:
cufflinks -p 4 -o wt1_clout wt1_thout/chrI.bam 
cufflinks -p 4 -o wt2_clout wt2_thout/chrI.bam 
cufflinks -p 4 -o wt1X_clout wt1X_thout/chrI.bam 
cufflinks -p 4 -o wt2X_clout wt2X_thout/chrI.bam 

echo "=============================================="
echo "cufflinks assemble is done ! "
echo "=============================================="



#2. Create a file called assemblies.txt that lists the assembly files for each sample. The file should contain the following lines:
#################
#wt1X_clout/transcripts.gtf
#wt1_clout/transcripts.gtf
#wt2X_clout/transcripts.gtf
#wt2_clout/transcripts.gtf
#################
# You can use a command to create this file:
rm -f assemblies.txt;
for i in `ls | grep clout`;do echo ./$i/transcripts.gtf >> assemblies.txt;done


#3. Merge all assemblies to one file containing merged transcripts: 

cuffmerge -g yeast_chrI_annotation.gff -s bowtie_index/YeastGenome.fa  -p 4 assemblies.txt  

echo "=============================================="
echo "cuffmerge is done ! "
echo "=============================================="


#Part III: Identify differentially expressed genes and transcripts 
#(we'll use the output from 100K_reads, not 10K_reads)

cuffdiff -o diff_out -b bowtie_index/YeastGenome.fa -p 4 -u merged_asm/merged.gtf ./wt1_thout/chrI.bam,./wt2_thout/chrI.bam   ./wt1X_thout/chrI.bam,./wt2X_thout/chrI.bam  

echo "=============================================="
echo "cuffdiff is done ! "
echo "=============================================="


#Part IV: Explore differential analysis results with CummeRbund TIMING variable
#(we'll use the output from 1M_reads, neither 100K nor 10K)

#Rscript plot_DE_chart.R 


echo "=============================================="
echo "R plot is done ! "
echo "=============================================="










