#!/bin/bash

####################################
#Header: ncrna_matrix_worm.sh
#Collecting data for ncRNA in Worm
####################################

####################################
#Files: File requirement
#genome_ncrna.pl
####################################

####################################
#Variables:
ARRAY_DIR=/home1/zl222/Projects/modENCODE/worm/tilingarray
sRNA_DIR=/home1/zl222/Projects/sRNApgene/sRNAworm/eland
RNA_DIR=/home1/zl222/Projects/TFworm/rnaseq
RNA_DIR2=/home1/zl222/Projects/modENCODE/worm/rnaseq/Waterson
BIN_DIR=/home1/zl222/Projects/ncRNA/Worm/matrix
GFF_DIR=/home1/zl222/Projects/modENCODE/worm/gff/processed
####################################

####################################
#Section: I. Collecting Matrix Data
###################################

#Topic: 1.Average singnals from tilting array (6 stages)
#option=-array
###################################
if [ "$1" = "-array" ]; then
	for i in emb-reference L2 L2poA L3 L4 YA 
	do	
		oldwig=$ARRAY_DIR/$i"_pm_mm_norm_min1.0_smooth110_log2.wig"
		mywig=$ARRAY_DIR/$i"_pm_mm_norm_min1.0_smooth110_log2.w25.wig"
		#awk '{if(NR>1) print $1"\t"$2"\t"$2+25-1"\t"$4}' $oldwig > $mywig
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.plus $mywig |  ncrna_matrix.pl -arraywig $mywig | awk '{print $1"\t"$2"\t"$3"\tarray\t"$4"\t+"}' > array.$i.bed
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.minus $mywig |  ncrna_matrix.pl -arraywig $mywig | awk '{print $1"\t"$2"\t"$3"\tarray\t"$4"\t-"}' >> array.$i.bed
 		#BinOverlapper --noHeader /home1/zl222/Projects/gff/worm/silverRNA.coor  $mywig | ncrna_matrix.pl -arraywig $mywig | awk '{print $1"\t"$2"\t"$3"\tRNA\t"$4"\t."}' > silverRNA.$i.bed
	done

#Topic: 2.Average singnals from small RNAseq in 11 stages
#option=-srna
###################################
elif [ "$1" = "-srna" ]; then
	dir=( Egg_s_3 L1_s_5 L2_s_6 L3_s_7 L4_s_8 yMale_dpyhim_s_6 Adult_s_1 DAY0_spe-9_s_1 DAY5_spe-9_s_2 DAY8_spe-9_s_3 DAY12_spe-9_s_5 )
	for i in ${dir[@]} ; do 
		SGR=$sRNA_DIR/$i.plus.sgr 	
		cat $BIN_DIR/bins.plus |mapsrna.pl -seqsgr $SGR > srna.$i.bed
		SGR=$sRNA_DIR/$i.minus.sgr 	
		cat $BIN_DIR/bins.minus |mapsrna.pl -seqsgr $SGR >> srna.$i.bed
	done

#Topic: 3.Average singnals from poly RNAseq in emb and starved L1
#option=-rnaseq
###################################
elif [ "$1" = "-rnaseq" ]; then
	dir=( PHA4EMB  PHA4L1 )
	for i in ${dir[@]} ; do 
		SGR=$RNA_DIR/$i.unique.wig 	
		cat $BIN_DIR/bins |mapsrna.pl -seqwigf $SGR > rnaseq.$i.bed
	done
elif [ "$1" = "-rnaseq2" ]; then
	dir=( L2 L3 L4 YA )
	for i in ${dir[@]} ; do 
		SGR=$RNA_DIR2/sgr/$i.sgr	
		cat $BIN_DIR/bins |mapsrna.pl -seqsgr $SGR > rnaseq.$i.bed
	done


#Topic: 4.RPKM for poly and short RNAseq
#option=-rpkm
###################################
elif [ "$1" = "-rpkm" ]; then
	dir=( Egg_s_3 L1_s_5 L2_s_6 L3_s_7 L4_s_8 yMale_dpyhim_s_6 Adult_s_1 DAY0_spe-9_s_1 DAY5_spe-9_s_2 DAY8_spe-9_s_3 DAY12_spe-9_s_5 )
	for i in ${dir[@]} ; do 
		mapsrna.pl -rpkm  srna.$i.bed $sRNA_DIR/mapped_reads > srna.$i.rpkm.bed
	done
	dir=( PHA4EMB  PHA4L1 )
	for i in ${dir[@]} ; do 
		mapsrna.pl -rpkm rnaseq.$i.bed $RNA_DIR/mapped_reads > rnaseq.$i.rpkm.bed
	done
elif [ "$1" = "-rpkm2" ]; then
	dir=( L2 L3 L4 YA )
	for i in ${dir[@]} ; do 
		mapsrna.pl -rpkm  rnaseq.$i.bed $RNA_DIR2/mapped_reads_no_rDNA > rnaseq.$i.rpkm.bed
	done


#Topic: 5 collect the matrix
#option=-matrix
###################################
elif [ "$1" = "-matrix" ]; then
	#	cat ../pairs/all.id.sense ../pairs/all.id.rev | awk '{print $1"\t"$2"\t"$3"\tid\t"$5"\t"$4}'  |sort -k6r,6 -k1,1 -k2,2n -k3,3n >  ../matrix/identities.bed
	#cat bins  | awk '{print $1"\t"$2"\t"$3"\tl\t"$3-$2+1"\t"$4}'   >  length.bed
	dir0=( length.bed identities.bed dyn.both.bed rnaz.both.bed ) 
	dir1=( array.emb-reference.bed array.L2.bed array.L2poA.bed array.L3.bed array.L4.bed array.YA.bed )
	dir2=( rnaseq.PHA4EMB.rpkm.bed  rnaseq.PHA4L1.rpkm.bed )
	dir3=( srna.Egg_s_3.rpkm.bed srna.L1_s_5.rpkm.bed srna.L2_s_6.rpkm.bed srna.L3_s_7.rpkm.bed srna.L4_s_8.rpkm.bed srna.yMale_dpyhim_s_6.rpkm.bed srna.Adult_s_1.rpkm.bed srna.DAY0_spe-9_s_1.rpkm.bed srna.DAY5_spe-9_s_2.rpkm.bed srna.DAY8_spe-9_s_3.rpkm.bed srna.DAY12_spe-9_s_5.rpkm.bed )

	paste ${dir0[@]}  ${dir1[@]} ${dir2[@]}   |awk '{for(i=5;i<=72;i+=6) {printf ",%s",$i} print "" }' > foo1 
	paste ${dir3[@]}  |awk '{for(i=5;i<=66;i+=6) {printf ",%s",$i} print "" }' > foo2

	paste bins foo1 foo2 |sed 's/\t/,/g' |sed 's/,,/,/g' > features.matrix	

#Topic: 6 collect the matrixheader
#option=-matrixheader
###################################
elif [ "$1" = "-matrixheader" ]; then
	dir0=( length.bed identities.bed dyn.both.bed rnaz.both.bed ) 
	dir1=( array.emb-reference.bed array.L2.bed array.L2poA.bed array.L3.bed array.L4.bed array.YA.bed )
	dir2=( rnaseq.PHA4EMB.rpkm.bed  rnaseq.PHA4L1.rpkm.bed )
	dir3=( srna.Egg_s_3.rpkm.bed srna.L1_s_5.rpkm.bed srna.L2_s_6.rpkm.bed srna.L3_s_7.rpkm.bed srna.L4_s_8.rpkm.bed srna.yMale_dpyhim_s_6.rpkm.bed srna.Adult_s_1.rpkm.bed srna.DAY0_spe-9_s_1.rpkm.bed srna.DAY5_spe-9_s_2.rpkm.bed srna.DAY8_spe-9_s_3.rpkm.bed srna.DAY12_spe-9_s_5.rpkm.bed )

for i in 'chromosome' 'start' 'end' 'strand'   ${dir0[@]}  ${dir1[@]} ${dir2[@]}  ${dir3[@]} ; do 
  		echo -n "$i,"	
	done
#Topic: 7 add matrix
#option=-addm
###################################
elif [ "$1" = "-addm" ]; then
	#dir0=( ../output/igb/zscore.bed ../output/igb/sci.bed ../output/igb/tblastx.bed ) 
#	paste ${dir0[@]}   |awk '{for(i=5;i<=18;i+=6) {printf ",%s",$i} print "" }' > foo1 
	cut -f 5 ../output/igb/gc.bed > foo1 
	mv -f  features.matrix features.matrix.bak
	paste features.matrix.bak foo1 |sed 's/\t/,/g' |sed 's/,,/,/g' > features.matrix	
	
#Topic: 7 Fix extreme values
#option=-fix
###################################
elif [ "$1" = "-fix" ]; then
#	awk '{if($1> -0.9294948+14.59844*2 ) print "28.26739"; else if($1< -0.9294948-14.59844*2) print "-30.12637"; else print $1 }' zscore >  zscore_fix
#	paste tblastx identity | awk '{if($2==0) print 0 ; else print $1/$2}' > tblastx_fix
#	awk '{if($1>1) print "1"; else print $1}' sci > sci_fix
	paste zscore_fix sci_fix tblastx_fix > foo1 
	mv -f  features.matrix features.matrix.bak
	paste features.matrix.bak foo1 |sed 's/\t/,/g' |sed 's/,,/,/g' > features.matrix	

#Topic: 8 Get Maximum values of expressions
#option=-max
###################################
elif [ "$1" = "-max" ]; then

awk -F ',' '{max=$9; for(i=10;i<=14;i++){if(max<$i) max=$i} print max}' features.matrix > array.max
awk -F ',' '{max=$15; for(i=16;i<=16;i++){if(max<$i) max=$i} print max}' features.matrix > rnaseq.Snyder.max
awk -F ',' '{max=$17; for(i=18;i<=27;i++){if(max<$i) max=$i} print max}' features.matrix > srna.max


###################################
#Section: II. Annotate the Matrix 
###################################

#Topic: 0.Sample 63/631 of tRNA
#option=-annot
###################################
elif [ "$1" = "-annot0" ]; then
	grep  tRNA ncRNA.gff > foo1 
	grep -v tRNA ncRNA.gff > foo2 
 	shuffle_file.pl foo1 | head -63 |sort -k1,1 -k4,4n -k5,5 > foo3
 	cat foo2 foo3 > ncRNAsample.gff
	\rm foo?
#Topic: 1.Annotate the matrix
#option=-annot
###################################
elif [ "$1" = "-annot1" ]; then
#for i in exonic  intronic  ; do
#	for i in predicted_ncRNA; do
#	for i in pseudogene; do
#	for i in otherRNA; do
#	for i in ncRNAsample; do
	for j in novel_Masa; do
		i=$GFF_DIR/$j
		coverage=0.5
#cut -f 1,4,5,7  $i.gff |grep + > $i.plus
#		cut -f 1,4,5,7  $i.gff |grep - > $i.minus
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.plus $i.plus |  ncrna_matrix.pl -annot  $i.gff + $coverage 2ways > $j.bed
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.minus $i.minus |  ncrna_matrix.pl -annot  $i.gff - $coverage 2ways >> $j.bed
		echo "done"
	done	
elif [ "$1" = "-annot2" ]; then
	for j in intergenic ; do
		i=$GFF_DIR/$j
		sed s/+/-/g  intergenic_original.gff > foo.gff
		cat intergenic_original.gff foo.gff > intergenic.gff
		cut -f 1,4,5,7  $i.gff |grep + > $i.plus
		cut -f 1,4,5,7  $i.gff |grep - > $i.minus
		coverage=0.9
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.plus $i.plus |  ncrna_matrix.pl -annot  $i.gff + $coverage > $i.bed
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.minus $i.minus |  ncrna_matrix.pl -annot  $i.gff - $coverage >> $i.bed
	done	
elif [ "$1" = "-annot3" ]; then
#for j in ccds ; do
	for j in genefinder_CDS genefinder_exon genefinder_intron twinscan_coding_exon  wormbase_exon ; do
		i=$GFF_DIR/$j
		coverage=0.9
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.plus $i.plus |  ncrna_matrix.pl -annot  $i.gff + $coverage > $j.bed
		BinOverlapper --noHeader --minOverlap=10  $BIN_DIR/bins.minus $i.minus |  ncrna_matrix.pl -annot  $i.gff - $coverage >> $j.bed
	done	


#Topic: 2 collect the annotation
#option=-ca
###################################
elif [ "$1" = "-ca" ]; then
	dir=( genefinder_CDS.bed genefinder_exon.bed genefinder_intron.bed twinscan_coding_exon.bed  wormbase_exon.bed )
#	dir=( ncRNA.bed predicted_ncRNA.bed  intergenic.bed  exonic.bed  intronic.bed ccds.bed )
	paste ${dir[@]}  |awk '{for(i=4;i<=30;i+=6) {printf ",%s",$i} print "" }' > foo
#	paste ../matrix/bins foo |sed 's/\t/,/g' |sed 's/,,/,/g' > annot.matrix	
#add=pseudogene.bed
#	add=novel_Masa.bed
#	cut -f 4 $add > foo
	paste annot.matrix foo |sed 's/\t/,/' | sed 's/,,/,/g' > foo1
	mv -f annot.matrix annot.matrix.bak
	mv -f foo1 annot.matrix
	

elif [ "$1" = "-test" ]; then
	mydir=/home1/zl222/Projects/ncRNA/Worm/training
	for m in intergenic exonic; do
	cut -f 1,4,5,7 $mydir/$m.gff | mapsrna.pl -seqwigf $SGR > test.$m.bed
	for i in Egg_s_3 ; do 
		SGR=$sRNA_DIR/$i.plus.sgr 	
		grep + $mydir/$m.gff|cut -f 1,4,5,7 |mapsrna.pl -seqsgr $SGR > testsrna.$m.bed
		SGR=$sRNA_DIR/$i.minus.sgr 	
		grep - $mydir/$m.gff|cut -f 1,4,5,7 |mapsrna.pl -seqsgr $SGR >> testsrna.$m.bed
	done

done
	
	

#Topic: 3.Calculating sensitivity
#option=-stat
###################################
elif [ "$1" = "-stat" ]; then
#paste L2.array.bed  L2poA.array.bed  L3.array.bed  L4.array.bed  YA.array.bed |cut -f 1-3,5,11,17,23,29  > expressed 

 	BinOverlapper --noHeader /home1/zl222/Projects/gff/worm/silverRNA.coor  | ncrna_matrix.pl -arraywig $mywig | awk '{print $1"\t"$2"\t"$3"\tRNA\t"$4"\t."}' > silverRNA.$i.bed
	
for j in dyn rnaz; do 
	paste silverRNAset.bed $j.both.bed expressed | cut -f 1-4,11,16,17,18,19,20 > $j.matrix 
	cat $j.matrix | ncrna_matrix.pl -stat > $j"_exp.stat"	
	for i in `seq 0 7`;do awk -v t=$i '{if($2==t) print $0}' $j"_exp.stat" > $j"_exp"$i.stat;done
	for i in  0 `seq 10 20 90` 95  99 ;do awk -v t=$i '{if($1==t) print $0}' $j"_exp.stat" > $j$i"_exp.stat" ;done
done
	
###################################
#Section: III. Parse kevin's supervised training
###################################

#Topic: 1.ROC curves
#option=-roc
###################################
elif [ "$1" = "-roc" ]; then
	for i in 7_7_7; do
		dir=with_gc
		dir1=kevin/$dir/$i
		dir2=results/$dir/$i
		for j in SVM_with_2nd_degree_polynomial_kernel Bayes_Net Logistic_Regression Random_Forest Decision_Tree Naive_Bayes ; do
#	tail -439825 $dir1/$j.txt |  sort -t $'\t' -k8,8 -k5,5 -k6,6n -k7,7n > $dir2/$j.out
#		awk '{if($4~/ncRNA_sampled/ || $4~/exon_CDS/ || $4~/intergenic_gold/) print $0}' $dir2/$j.out > $dir2/$j.training.out
#cut -f 1,4 $dir2/$j.training.out | ncrna_matrix.pl -roc ncRNA_sampled > $dir2/$j.ncrna.roc
		cut -f 2,4 $dir2/$j.training.out | ncrna_matrix.pl -roc exon_CDS > $dir2/$j.cds.roc
		done
	done
#Topic: 1a.other ROC curves
#option=-roc
###################################

elif [ "$1" = "-roc1" ]; then
	feature_dir=/home1/zl222/Projects/ncRNA/Worm/training/features_raw
	for dir in 7_7_7; do
		cd results/with_gc/$dir
		#cut -f 4 Bayes_Net.out > annot4training
		for i in array.max rnaseq.Snyder.max srna.max;do  
#paste $feature_dir/$i  annot4training |awk -F '\t' '{if($2~/ncRNA_sampled/||$2~/exon_CDS/||$2~/intergenic_gold/) print $1"\t"$2}' | ncrna_matrix.pl -lognormroc ncRNA_sampled  > $i.ncrna.roc
			paste $feature_dir/$i  annot4training |awk -F '\t' '{if($2~/ncRNA_sampled/||$2~/exon_CDS/||$2~/intergenic_gold/) print $1"\t"$2}' | ncrna_matrix.pl -lognormroc exon_CDS  > $i.cds.roc
		done
	
#		for i in identity;do  
#			paste $feature_dir/$i  annot4training |awk -F '\t' '{if($2~/ncRNA_sampled/||$2~/exon_CDS/||$2~/intergenic_gold/) print $1"\t"$2}' | ncrna_matrix.pl -roc ncRNA_sampled  > $i.ncrna.roc
#		done
#		for i in dyn.both rnaz.both;do  
#			paste $feature_dir/$i  annot4training |awk -F '\t' '{if($2~/ncRNA_sampled/||$2~/exon_CDS/||$2~/intergenic_gold/) print $1/1000"\t"$2}' | ncrna_matrix.pl -roc ncRNA_sampled  > $i.ncrna.roc
#		done
#		for i in rnaseq.PHA4EMB srna.egg;do  
#			paste $dir2/$i annot4training |awk -F '\t' '{if($2~/ncrna/||$2~/neg/||$2~/cds/) print $1/200"\t"$2}' | ncrna_matrix.pl -roc ncrna  > $i.ncrna.roc
#		done
#		for i in array.emb ;do  
#			paste $dir2/$i annot4training |awk -F '\t' '{if($2~/ncrna/||$2~/neg/||$2~/cds/) print $1/20"\t"$2}' | ncrna_matrix.pl -roc ncrna  > $i.ncrna.roc
#		done
		
		cd ../../..
	done
#Topic: 2.Merge the annotations for each group
#option=-mannot
###################################
elif [ "$1" = "-mannot" ]; then
#	awk -F ',' '{if($5~/non/ && $14~/non/ &&  $16~/non/ && $17~/non/ && $18~/non/ && $19~/non/ && $21~/non/ && $22~/non/ ) print $1"\t"$15 }'  training/annot.matrix > annot.final 
	mv -f training/annot.final training/annot.final.bak
	perl -n -e '@t=split ("[,\n]",$_);
		if($t[13] ne 'non') {print "ncRNA_sampled";}
		elsif ($t[4] ne 'non') {print "ncRNA";} 
		elsif($t[9] ne 'non' && $t[12] ne 'non') {print "exon_CDS";}
		elsif($t[9] ne 'non') {print "exon_coding_exons";}
		elsif($t[12] ne 'non') {print "exon_ccds2";}
		elsif($t[7] ne 'non') {print "exonic";}
		elsif($t[5] ne 'non') {print "ncRNA_predicted_21Uetc";}
		elsif($t[14] ne 'non') {print "ncRNA_Masa";}
		elsif($t[11] ne 'non') {print "ncRNA_blast";}
		elsif($t[8] ne 'non') {print "intronic";}
		elsif($t[10] ne 'non') {print "pseudogene";}
		elsif($t[6] ne 'non') {print "intergenic_gold";}
		
		elsif($t[19] ne 'non') {print "wormbase_unconfirmed_exon";}
		elsif($t[15] ne 'non' && $t[18] ne 'non' ) {print "genefinder_twinscan_CDS";}
		elsif($t[15] ne 'non') {print "genefinder_CDS";}
		elsif($t[18] ne 'non') {print "twinscan_coding_exon";}
		elsif($t[16] ne 'non') {print "genefinder_exon";}
		elsif($t[17] ne 'non') {print "genefinder_intron";}
		
		else {print "intergenic";}
		print "\n"; ' training/annot.matrix > training/annot.final

#Topic: 2.Find positive ration for each group
#option=-pr
###################################
elif [ "$1" = "-pr" ]; then
	dir=rnaseq_preprocessed_with_cds_ncrna_sampled/7_7_7
	for j in Logistic_Regression Random_Forest ;do
	for i in  ncRNA_sampled ncRNA exon_CDS exon_coding_exons exon_ccds2 exonic ncRNA_predicted_21Uetc ncRNA_blast intronic pseudogene intergenic_gold intergenic; do
	for x in 0 0.2 0.5 0.9  ; do
			echo -n "$j	$i	"
			paste training/annot.final results/$dir/$j.out |awk -v var=$i '{if($1==var) print $0}'| sed 's/,/\t/g' | awk -F '\t' -v var=$x '{all++;if($2>=var)sum++}END{print "cutoff_"var":\t"sum"\t"sum/all*100"%"}' 

#paste results/$dir/$j.out training/annot.matrix | sed 's/,/\t/g' | awk -F '\t' -v var=$x '{if($1>=var && $13~/non/ && $14~/non/ &&  $16~/non/ && $17~/non/ && $18~/non/ && $19~/non/ && $21~/non/ && $22~/non/ ) print $1"\t"$15 }' |wc -l  
#paste results/$dir/$j.out training/annot.matrix | sed 's/,/\t/g' | awk -F '\t' -v var=$x '{if($1>=var && $14!~/non/ ) print $1"\t"$15 }' |wc -l 
#paste results/$dir/$j.out training/annot.matrix | sed 's/,/\t/g' | awk -F '\t' -v var=$x '{if($1>=var && $18!~/non/ ) print $1"\t"$15 }' |wc -l 
   done
	done   
   done
#Topic: 3.Distribution
#option=-distr
###################################
elif [ "$1" = "-distr" ]; then
	dir=results/rnaseq_preprocessed_with_cds_ncrna_sampled/7_7_7
	for j in Logistic_Regression Random_Forest ;do
#for i in  ncRNA_sampled exon_CDS intergenic_gold intergenic; do
	for i in  intergenic_gold ; do
			paste $dir/$j.out training/annot.final |awk -F '\t' -v var=$i '{if($9~var) print $0}' > R/$j.$i.distr
			cut -f 1 R/$j.$i.distr > R/$j.$i.distr.prob
	done
	done

###################################
#Section: III. Select the Candidates 
###################################
#Topic: 1. making gff files for bins
#option=-mgff
elif [ "$1" = "-mgff" ]; then
	awk '{print $1}' foo > foo.gff

#Topic: 2. generate 
#option=-pick
elif [ "$1" = "-pick" ]; then
	out=results/rnaseq_preprocessed_with_cds_ncrna_sampled/7_7_7/Random_Forest.out
	paste training/annot.final $out training/features.matrix |awk '{if ($1~/intergenic/ && $2 >=0.9 ) print $1"\t"$2"\t"$9}' | sed 's/,/\t/g' > ncRNA.list 


#end of the jobs	
###################################
else 
	echo "Type in options"
fi


