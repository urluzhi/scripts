#!/bin/bash

####################################
#Variables:
Bin_DIR=/home1/zl222/Projects/ncRNA/worm/bins/loc
TAR_DIR=$Bin_DIR
my_tars=( bins )
GFF_DIR=/home1/zl222/Projects/modENCODE/worm/gff/processed
ARRAY_wig=/home1/zl222/Projects/modENCODE/worm/expression/tilingarray/wig
sRNA_sgr=/home1/zl222/Projects/modENCODE/worm/expression/srnaseq/sgr
RNA_wig=/home1/zl222/Projects/modENCODE/worm/expression/rnaseq/Snyder/wig
RNA_sgr=/home1/zl222/Projects/modENCODE/worm/expression/rnaseq/Waterson/sgr
sRNA_DIR=/home1/zl222/Projects/modENCODE/worm/expression/srnaseq
RNA_DIR=/home1/zl222/Projects/modENCODE/worm/expression/rnaseq/Snyder
RNA_DIR2=/home1/zl222/Projects/modENCODE/worm/expression/rnaseq/Waterson
Prediction_DIR=/home1/zl222/Projects/ncRNA/worm/results/new_rnaseq/7_7_7
Prediction_TAR_DIR=/home1/zl222/Projects/ncRNA/worm/results/tar_without_length
Rfam_DIR=/home1/zl222/Projects/ncRNA/worm/rfam
Feature_DIR=/home1/zl222/Projects/ncRNA/worm/training
Pha4_DIR=/home1/zl222/Projects/ncRNA/worm/Chipseq/pha4
PolII_DIR=/home1/zl222/Projects/ncRNA/worm/Chipseq/polII
all_arrays="N2EE
N2LE
L1-NDT
L2
L3
L4
YA
JKL4-NDT
MALE
Db-24
Db-48
Ef-24
OP50-24
OP50-48
Pl-24
L2-reference
L3-L4-reference
emb-AVA
emb-GABA
emb-bwm-v2
emb-coelomocytes
emb-dop
emb-intestine
emb-panneural
emb-A-class
emb-GLP
emb-hypodermis
emb-reference
YA-ref
L2poA
GON
L2-glr
L2-A-class
L2-GABA_neurons
L2-bwm
L2-excretory_cell
L2-intestine
L2-panneural
L3-L4-PVD_OLL
L3-L4-dop
YA-CEPsh
"
###################################
#Section: I. Annotate the TARs 
###################################
#Topic: 1a.Annotate the TAR
#option=-annot1 (50% 2ways)
###################################
if [ "$1" = "-annot1" ]; then
for j in ncRNAsample ncRNA  predicted_ncRNA novel_Masa otherRNA; do
	for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0.5
		ways=0
		BinOverlapper $i.plus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage $ways > annot/$k-$j.annot
		BinOverlapper $i.minus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways >> annot/$k-$j.annot
	done	
	done

#Topic: 1b.Annotate the TAR
#option=-annot2 (50% 1way)
###################################
elif [ "$1" = "-annot2" ]; then
for j in exonic intronic pseudogene genefinder_CDS genefinder_exon genefinder_intron twinscan_coding_exon  wormbase_exon; do
	for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0.5
		ways=1
		BinOverlapper $i.plus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage $ways > annot/$k-$j.annot
		BinOverlapper $i.minus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways >> annot/$k-$j.annot
	done	
	done
#Topic: 1bb.Annotate the TAR
#option=-annot2b (50% 1way)
###################################
elif [ "$1" = "-annot2b" ]; then
for j in five_prime_UTR three_prime_UTR; do
	for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0.5
		ways=1
		BinOverlapper $i.plus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage $ways > annot/$k-$j.annot
		BinOverlapper $i.minus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways >> annot/$k-$j.annot
	done	
	done

#Topic: 1c.Annotate the TAR
#option=-annot3 (90% 1way)
###################################
elif [ "$1" = "-annot3" ]; then
	for j in pgene_exon  ; do
#	for j in coding_exons ccds intergenic ; do
		for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0.9
		ways=1
		BinOverlapper $i.plus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage  $ways > annot/$k-$j.annot
		BinOverlapper $i.minus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways >> annot/$k-$j.annot
	done	
	done

#Topic: 1d.Annotate the TAR
#option=-annot4 (1 base overlap for exclusive(exons) and for inclusive (introns,pseudogene..))
###################################
elif [ "$1" = "-annot4" ]; then
	for j in exclusive inclusive; do
		for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0
		ways=0
		BinOverlapper $i.plus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage  $ways > annot/$k-$j.annot
		BinOverlapper $i.minus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways >> annot/$k-$j.annot
	done	
	done


#Topic: 2a collect the annotation into matrix
#option=-ca
###################################
elif [ "$1" = "-ca" ]; then
	for j in ${my_tars[@]}; do
		a=""
		for i in  ncRNAsample ncRNA  coding_exons ccds exonic predicted_ncRNA novel_Masa otherRNA  intronic pseudogene intergenic wormbase_exon twinscan_coding_exon genefinder_CDS genefinder_exon genefinder_intron inclusive exclusive three_prime_UTR five_prime_UTR pgene_exon;do 
			a=$a" annot/$j-$i.annot"
   		done
		paste $a |sed 's/\t/,/g' | sed 's/,,/,/g' > annot/$j.annot.matrix
	done	
		

#Topic: 2b collapse the annotation into one column
#option=-mannot
###################################
elif [ "$1" = "-mannot" ]; then
	for j in ${my_tars[@]}; do
	perl -n -e '@t=split ("[,\n]",$_);
		if($t[0] ne 'non') {print "ncRNA_sampled";}
		elsif ($t[1] ne 'non') {print "ncRNA";} 
		elsif($t[2] ne 'non' && $t[3] ne 'non') {print "exon_CDS";}
		elsif($t[2] ne 'non') {print "exon_coding_exons";}
		elsif($t[18] ne 'non') {print "three_prime_UTR";}
		elsif($t[19] ne 'non') {print "five_prime_UTR";}
		elsif($t[4] ne 'non') {print "exonic";}
		elsif($t[5] ne 'non') {print "ncRNA_predicted_21Uetc";}
		elsif($t[6] ne 'non') {print "ncRNA_Masa";}
		elsif($t[7] ne 'non') {print "ncRNA_blat_RNAi_etc";}
		elsif($t[20] ne 'non') {print "pgene_exon";}
		elsif($t[8] ne 'non') {print "intronic";}
		elsif($t[9] ne 'non') {print "pseudogene";}
		elsif($t[10] ne 'non') {print "intergenic_gold";}
		
		elsif($t[11] ne 'non') {print "wormbase_unconfirmed_exon";}
		elsif($t[12] ne 'non') {print "twinscan_coding_exon";}
		elsif($t[14] ne 'non') {print "genefinder_exon";}
		elsif($t[15] ne 'non') {print "genefinder_intron";}
		elsif($t[3] ne 'non') {print "CCDS";}
		elsif($t[13] ne 'non') {print "genefinder_CDS";}
		else {print "intergenic";}
		print "\n"; ' annot/$j.annot.matrix > annot/$j.annot.final
	
	
	perl -n -e '@t=split ("[,\n]",$_);
		if($t[16] ne 'non') {print $t[16];}
		else {print "non";}
		print "\t";
		if($t[17] eq 'non') {print "non";}
		else {print "$t[17]";}
		print "\n"; ' annot/$j.annot.matrix > annot/$j.annot.novel
	
	done

#Topic: 2c summary the annotations 
#option=-sumannot
elif [ "$1" = "-sumannot" ]; then
	for j in ${my_tars[@]}; do
		sort  annot/$j.annot.final |uniq -c |awk '{
			all+=$1
		
			if($2~/^ncRNA_sampled/) {ncrna+=$1; ncrna_sampled+=$1 }
			else if($2~/^exon_CDS/) {exon+=$1; exon_CDS+=$1 }
			else if($2~/^three_prime_UTR/) {exon+=$1; three_prime_UTR+=$1 }
			else if($2~/^five_prime_UTR/) {exon+=$1; five_prime_UTR+=$1 }
			else if($2~/^pgene_exon/) {pgene_exon+=$1;pgene+=$1 }
			else if($2~/^exon/) {exon+=$1 }
			else if($2~/intron/) {intron+=$1}
			else if($2~/intergenic_gold/) {inter_gold+=$1}
			else if($2~/genefinder/ || $2~/twinscan/ || $2~/wormbase_unconfirmed/) {exon2+=$1 }
			else if($2~/CCDS/) {intron+=$1}
			else if($2~/pseudogene/) {pgene+=$1}
			else if($2~/ncRNA_blat/) {ncrna_blat+=$1}
			else if($2~/ncRNA/) {ncrna+=$1}
			else {inter+=$1}
			}
			END {
				print "Type\t# of TARs\tSource"
				print "ALL\t"all
				print "Exon/Confirmed Exon/Confirmed Exon_CDS\t"exon+exon2"/"exon"/"exon_CDS"\tWormbase(confirmed and unconfirmed); Twinscan; GeneFinder"
				print "Five_prime_end/Three_prime_end\t"five_prime_UTR"/"three_prime_UTR"\tWormbase"
				print "ncRNA/Sampled ncRNA\t"ncrna+ncrna_blat"/"ncrna_sampled"\tWormbase(known types,unknown type,21U-RNA); novel miRNA/21U-RNA from Masa; blat,RNAi"
				print "Pseudogene/Pseudogene_Exon\t"pgene"/"pgene_exon"\tPseudopipe"
				print "Intron\t"intron"\tWormbase(confirmed); GeneFinder"
				print "Intergenic/Selected intergenic region\t"inter"/"inter_gold"\tUnannotated Region/Ladeana selected intergenic set"
			}'  > annot/$j.annot.sum
	done
#Topic: 2e.if assigned by array TAR TAR
#option=-assign (1 base overlap )
###################################
elif [ "$1" = "-assign" ]; then
	for j in $other_prediction; do
		for k in ${my_tars[@]} ; do 
		coverage=0
		ways=0
		BinOverlapper $OTHER_PREDICTION_DIR/$j.loc  $TAR_DIR/$k.loc |  ncrna_matrix.pl -assign  $j.csv  $coverage  $ways > annot/$k-$j.assign
	done	
	done




###################################
#Section: II. Annotate the antisense strand 
###################################

#Topic: 1.Antisense of the TAR
#option=-anti (50% 2ways)
###################################
elif [ "$1" = "-anti1" ]; then
for j in exclusive; do
	for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0.5
		ways=0
		BinOverlapper $i.minus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways > annot/$k-$j.anti
		BinOverlapper $i.plus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage $ways >> annot/$k-$j.anti
	done	
	done


#Topic: 2.Antisense of the TAR
#option=-anti (1 base)
###################################
elif [ "$1" = "-anti2" ]; then
for j in exclusive inclusive; do
	for k in ${my_tars[@]} ; do 
		i=$GFF_DIR/$j
		coverage=0
		ways=0
		BinOverlapper $i.minus  $TAR_DIR/$k.plus |  ncrna_matrix.pl -annot_simple  $i.gff - $coverage $ways > annot/$k-$j.anti2
		BinOverlapper $i.plus $TAR_DIR/$k.minus |  ncrna_matrix.pl -annot_simple  $i.gff + $coverage $ways >> annot/$k-$j.anti2
	done	
	done




	
##################################t#
#Section: III. Collecting Feature Matrix 
###################################

#Topic: 1.Average singnals from tilting array (6 stages)
#option=-array
###################################
elif [ "$1" = "-array" ]; then
	for i in $all_arrays; do 
		mywig=$ARRAY_wig/$i.norm.wig
		for j in ${my_tars[@]}; do
			echo $i
			BinOverlapper  $mywig $TAR_DIR/$j.plus |  ncrna_matrix.pl -arraywig $mywig | awk '{print $4}' > features/array.$j.$i
			echo $i.minus
			BinOverlapper  $mywig $TAR_DIR/$j.minus |  ncrna_matrix.pl -arraywig $mywig | awk '{print $4}' >> features/array.$j.$i
		done
	done

#Topic: 2.Average singnals from small RNAseq in 11 stages
#option=-srna
###################################
elif [ "$1" = "-srna" ]; then
	dir=( Egg_s_3 L1_s_5 L2_s_6 L3_s_7 L4_s_8 yMale_dpyhim_s_6 Adult_s_1 DAY0_spe-9_s_1 DAY5_spe-9_s_2 DAY8_spe-9_s_3 DAY12_spe-9_s_5 )
	for i in ${dir[@]} ; do 
		for j in ${my_tars[@]}; do
			SGR=$sRNA_sgr/$i.plus.sgr 	
			cat $TAR_DIR/$j.plus |mapsrna.pl -seqsgr $SGR | cut -f 5 > features/srna.$j.$i
			SGR=$sRNA_sgr/$i.minus.sgr 	
			cat $TAR_DIR/$j.minus |mapsrna.pl -seqsgr $SGR| cut -f 5  >> features/srna.$j.$i
		done
	done

#Topic: 3.Average singnals from poly RNAseq in emb and starved L1
#option=-rnaseq
###################################
elif [ "$1" = "-rnaseq" ]; then
	dir=( PHA4EMB  PHA4L1 )
	for i in ${dir[@]} ; do 
		SGR=$RNA_wig/$i.unique.wig 	
		for j in ${my_tars[@]}; do
			cat $TAR_DIR/$j.plus $TAR_DIR/$j.minus |mapsrna.pl -seqwigf $SGR |cut -f 5 > features/rnaseq.$j.$i
		done
	done

#Topic: 3a .Average singnals from poly RNAseq (Waterson Group) 
#option=-rnaseq
###################################
elif [ "$1" = "-rnaseq2" ]; then
    dir=( L2 L3 L4 YA )
    for i in ${dir[@]} ; do
	    SGR=$RNA_sgr/$i.sgr    
		for j in ${my_tars[@]}; do
			cat $TAR_DIR/$j.plus $TAR_DIR/$j.minus |mapsrna.pl -seqsgr $SGR |cut -f 5 > features/rnaseq.$j.$i
		done
	done
							

	
#Topic: 4.DCPM for poly and short RNAseq
#option=-dcpm
###################################
elif [ "$1" = "-dcpm" ]; then
	for j in ${my_tars[@]}; do
	dir=( Egg_s_3 L1_s_5 L2_s_6 L3_s_7 L4_s_8 yMale_dpyhim_s_6 Adult_s_1 DAY0_spe-9_s_1 DAY5_spe-9_s_2 DAY8_spe-9_s_3 DAY12_spe-9_s_5 )
	for i in ${dir[@]} ; do 
		mapsrna.pl -dcpm  features/srna.$j.$i $sRNA_DIR/mapped_reads $i > features/srna.$j.$i.dcpm
	done
	dir=( PHA4EMB  PHA4L1 )
	for i in ${dir[@]} ; do 
		mapsrna.pl -dcpm features/rnaseq.$j.$i $RNA_DIR/mapped_reads_no_rDNA $i > features/rnaseq.$j.$i.dcpm
	done
	dir=( L2 L3 L4 YA )
	for i in ${dir[@]} ; do 
		mapsrna.pl -dcpm  features/rnaseq.$j.$i $RNA_DIR2/mapped_reads_no_rDNA $i > features/rnaseq.$j.$i.dcpm
	done
	done				   
					
#Topic: 5 collect the matrix
#option=-expression
###################################
elif [ "$1" = "-expression" ]; then
	cd features
	for j in ${my_tars[@]}; do
		features=''
		header=''
		n=4
		for i in  PHA4EMB  PHA4L1 L2 L3 L4 YA ;do
			let n=$n+1
			features=$features" rnaseq."$j"."$i".dcpm"
			header=$header" "$n".rnaseq."$i
		done
		for i in  Egg_s_3 L1_s_5 L2_s_6 L3_s_7 L4_s_8 yMale_dpyhim_s_6 Adult_s_1 DAY0_spe-9_s_1 DAY5_spe-9_s_2 DAY8_spe-9_s_3 DAY12_spe-9_s_5;do
			let n=$n+1
			features=$features" srna."$j"."$i".dcpm"
			header=$header" "$n".srna."$i
		done
		for i in $all_arrays;do
			let n=$n+1
			features=$features" array."$j"."$i
			header=$header" "$n".array."$i
		done		
		echo "1.chromosome 2.start 3.end 4.strand" $header |sed 's/ /,/g' > ../output/expression.header	
		paste $TAR_DIR/$j.loc  $features |sed 's/\t/,/g' > $j.expression
	done	
	cd ..
#Topic 6: expression maximum
#option=-expressioni_max
###################################
elif [ "$1" = "-expression_max" ]; then
 echo "1.polyA_RNAseq_max_all,2.small_RNAseq_max_all,3.Array_max_all,4.Array_max_totalRNA,5.Array_max_polyA,6.polyA_RNAseq_max_Snyder,7.polyA_RNAseq_max_Waterson" > features/expression.max.header
	cd features
	for j in ${my_tars[@]}; do
		cat $j.expression |awk -F ','  '{ 
			max=$5; for(i=6;i<=10;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$11; for(i=12;i<=21;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$22; for(i=23;i<=62;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$22; for(i=23;i<=50;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$51; for(i=52;i<=62;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$5; for(i=6;i<=6;i++) {if(max<$i)max=$i} printf "%f,",max; 
			max=$7; for(i=8;i<=10;i++) {if(max<$i)max=$i} printf "%f\n",max; 
			}' > $j.expression.max
	done	
	cd ..



#Topic: 7 Fix and normalize features
#option=-fix
###################################
elif [ "$1" = "-fix" ]; then
	cd features
	awk '{if($1> -0.9294948+14.59844*2 ) print "28.26739"; else if($1< -0.9294948-14.59844*2) print "-30.12637"; else print $1 }' zscore >  zscore_fix
	paste tblastx identity | awk '{if($2==0) print 0 ; else print $1/$2}' > tblastx_fix
	awk '{if($1>1) print "1"; else print $1}' sci > sci_fix
	cd ..
#Topic: 8 remove mtDNA
#option=-rm_mt
###################################
elif [ "$1" = "-rm_mt" ]; then
	for i in dyn.both  length     sci      tblastx      zscore identity  rnaz.both  sci_fix  tblastx_fix  zscore_fix GC; do
		paste loc/bins_with_MtDNA.loc ../bins_old/training/features/$i |grep -v DNA |cut -f 5 > features/$i	
	done

#Topic: 9 collect the matrix
#option=-matrix
###################################
elif [ "$1" = "-matrix" ]; then
	cd features 
	echo "00.Annotation,0.ID,1.chromosome,2.start,3.end,4.strand,5.length,6.GC%,7.identities,8.zscore_fix,9.sci_fix,10.tblastx_fix,11.RNAz,12.DynaFind" > features.header
	paste features.header expression.max.header | \
   	sed -r 's/[0-9]+\.//g'|sed 's/,/\t/g' | awk -F '\t' 'END{for(i=1;i<NF;i++) {printf i"."$i","} printf NF"."$NF"\n"}'> matrix.header
	cd ..
	cd features
	for j in ${my_tars[@]}; do
		paste ../annot/$j.annot.final ../output/$j.id length GC identity zscore_fix sci_fix tblastx_fix rnaz.both dyn.both $j.expression.max | \
		sed 's/\t/,/g' |sed 's/,,/,/g' > $j.matrix
		cat matrix.header $j.matrix > ../output/$j.training.csv
	done
	cd ..
	
#end of the jobs	
###################################
else 
	echo "Type in options"
fi


