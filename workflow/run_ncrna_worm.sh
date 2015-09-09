#!/bin/bash

####################################
#Header: run_ncrna_worm.sh
#Running ncRNA for Worm
#Running Dynafind and RNAz on Briggsae198 and Elegans170
#Fine Documents.
####################################


####################################
#Files: File requirement
#genome_ncrna.pl
####################################

####################################
#Variables: 
HOME=/home1/zl222/Projects/ncRNA/Worm
ALIGN1=briggsae
ALIGN2=elegans

SVMZ_MODEL=../SVMZ_MODEL/model
SVMZ_SCALER=../SVMZ_MODEL/scaler
AVE_MODEL=../AVE_SD_SVM_MODEL/scaledAverageOutput-all.model
AVE_SCALER=../AVE_SD_SVM_MODEL/averageOutput-all.scaler
SD_MODEL=../AVE_SD_SVM_MODEL/scaledSDOutput-all.model
SD_SCALER=../AVE_SD_SVM_MODEL/sdOutput-all.scaler
####################################


#if  [ $1 ]; then echo "type something"; exit; fi 

####################################
#Section: I. Prepare alignments
#prepare aligments (MUMmer) briggsae198-elegan170
###################################

#Topic: 1.Align with MUMmer
#option=0
###################################
if [ $1 = 0 ]; then
	cd /home1/zl222/Projects/genomes/worm/align2
	nucmer -p out -b 1600 -c 10 ../cbriggsae/c_briggsae.WS198.dna.fa ../celegans/elegans.WS170.dna.fa > brig-elegan.out 2> brig-elegan.log
	for i in chrI chrII chrIII chrIV chrV chrX; do for j in I II III IV V X MtDNA; do
		show-aligns out.delta $i $j > $i-$j.align
	done; done
#Topic: 2. make dir
#option=-mkdir	
###################################
elif [ "$1" = "-mkdir" ] ;then
	for i in sequences dynalign maf jobs dynalign/ct  splits pairs svm4ave_sd svmz single single/ct rnaz  rnaz/out  output
	do if [ -d $i ] ;then echo -n ""; else echo "making $i";mkdir $i; fi; done 
fi
####################################
#Section: II.Prepare the files for folding
####################################
for i in chrI chrII chrIII chrIV chrV chrX; do for j in I II III IV V X MtDNA; do
	chr=$i-$j
	ALIGN=alignments/$chr.align
	ALIGN_BOTH=alignments/"$chr"_both.align
	REV=alignments/"$chr"_rev.align
	SPLIT=splits/$chr.split
	REV_SPLIT=splits/"$chr"_rev.split
	PAIR=pairs/$chr.pair
	LOC=pairs/$chr.loc
	SEQ=sequences/$chr.seq
	MAF=maf/$chr.maf
#Topic: 1.revese strand of the alignments 
####################################
if [ "$1" = "-rev" ]; then
	echo  "make rev:"
    genome_worm.pl -revmum $<  $REV 
	cat $ALIGN $REV > $ALIGN_BOTH
#Topic: 2.split alignment into windows
####################################
elif [ "$1" = "-split" ]; then
    echo "make split:"
	genome_worm.pl -spmum $ALIGN $SPLIT > $SPLIT.log
    genome_worm.pl -spmum $REV $REV_SPLIT > $REV_SPLIT.log
#Topic: 3. making paired sequences and maf
####################################
elif [ "$1" = "-seq" ]; then
	echo "$chr"
	genome_worm.pl -seq $SPLIT $PAIR.sense $SEQ.sense  $MAF.sense> $PAIR.log
	genome_worm.pl -seq $REV_SPLIT $PAIR.rev $SEQ.rev $MAF.rev > $PAIR.rev.log
#Topic: 4.calc. coordinates of each window
####################################
elif [ "$1" = "-loc" ]; then
	date
	echo "$PAIR.sense"
	genome_worm.pl -loc $PAIR.sense $ALIGN > $LOC.sense
	date
	genome_worm.pl -loc $PAIR.rev $REV > $LOC.rev
	date
	
fi
done
done
###################################a
#Section: III. distribute the jobs
####################################
if [ "$1" = "-jobseq" ]; then
	for j in all.sense all.rev; do	
		SL=`wc -l $j |cut -d " " -f 1`
		let size=(SL/80/3+1)*3
		echo $size
		split  -d -l $size  $j  $j
	done
elif [ "$1" = "-jobmaf" ]; then
	for j in all.maf.sense all.maf.rev; do	
		SL=`wc -l $j |cut -d " " -f 1`
		let size=(SL/80/4+1)*4
		echo $size
		split  -d -l $size  $j  $j
	done
elif [ "$1" = "-pbs" ]; then
		i=$2
		j=rev
		echo "#!/bin/bash" 
		echo "#PBS -q gerstein -W group_list=gerstein" 
		echo "cd /home1/zl222/Projects/ncRNA/Worm" 
		echo "export DATAPATH=/home1/zl222/svn/dynalign/" 
		echo "date" 
		echo "mmfold sequences/all.$j$i > single/ct/all.$j$i.out" 
		echo "date" 
		echo "RNAz maf/all.maf.$j$i > rnaz/out/all.$j$i.out" 
		echo "date" 
		echo "mdynalign sequences/all.$j$i > dynalign/ct/all.$j$i.out"
		echo "date" 
fi
####################################
#Section: IV. Run SVMz
####################################
for i in 00 01 02 03 04 05 06 07 08 09  `seq 10 79`;do
for j in rev  ; do
	chr=all.$j$i
#Topic: 1.predict ave and sd for randomized seq.
####################################
if [ "$1" = "-svm-ave-sd" ]; then
		echo $chr
		date
		genome_worm2.pl -svminput sequences/$chr > svm4ave_sd/$chr.ave.svminput
		cp svm4ave_sd/$chr.ave.svminput svm4ave_sd/$chr.sd.svminput
		#calculate Ave. and sd. 
		svm-scale -r $AVE_SCALER svm4ave_sd/$chr.ave.svminput > svm4ave_sd/$chr.ave.svminput.scale 
		svm-predict svm4ave_sd/$chr.ave.svminput.scale  $AVE_MODEL  svm4ave_sd/$chr.ave.predict > svm4ave_sd/$chr.ave.log   
    	svm-scale -r $SD_SCALER  svm4ave_sd/$chr.sd.svminput > svm4ave_sd/$chr.sd.svminput.scale
		svm-predict svm4ave_sd/$chr.sd.svminput.scale $SD_MODEL svm4ave_sd/$chr.sd.predict > svm4ave_sd/$chr.sd.log 
		paste  svm4ave_sd/$chr.ave.predict svm4ave_sd/$chr.sd.predict > svm4ave_sd/$chr.ave_sd 
		date
#Topic: 2.inputs for SVMz
####################################
elif [ "$1" = "-svmzinput" ]; then
	SCI=svmz/$chr.sci
	SVM=svmz/$chr.svmz
	ZSCORE=svmz/$chr.zscore
	paste dynalign/ct/$chr.out single/ct/$chr.out | awk '{if(($4+$5)==0) {print $1"\t0"} else{ print $1"\t"$2/($4+$5)} }' > $SCI
	paste dynalign/ct/$chr.out svm4ave_sd/$chr.ave_sd | awk '{print $1"\t"($2-$3)/$4}' > $ZSCORE
	cut -d " " -f 2-9 svm4ave_sd/$chr.ave.svminput | awk '{print "1 "$0}' > $SVM.tmp1
	paste $ZSCORE $SCI |awk '{print " 9:"$2" 10:"$4}' > $SVM.tmp2
	paste $SVM.tmp1 $SVM.tmp2 > $SVM 
	rm $SVM.tmp1; rm $SVM.tmp2
#Topic: 3.run svmz 
####################################
elif [ "$1" = "-svmz" ]; then
	SVM=svmz/$chr.svmz
	svm-scale -r $SVMZ_SCALER $SVM > $SVM.scale
	svm-predict -b 1 $SVM.scale $SVMZ_MODEL $SVM.predict > $SVM.log
	genome_worm.pl -readsvm $SVM.predict $SVM.prob
####################################
#Section: V. Collect the result 
####################################

#Topic: 1. parse output of RNAz
####################################
elif [ "$1" = "-rnaz" ]; then
			grep "^##maf" maf/all.maf.$j$i |cut -d " " -f 4 > maf/all.id.$j$i
		    grep "SVM RNA-class probability:"  rnaz/out/all.$j$i.out |cut -d " " -f 5 > rnaz/out/all.$j$i.prob
#Topic: 2. collect output of mfold, dynalign and rnaz
####################################
elif [ "$1" = "-collect" ]; then
		paste maf/all.id.$j$i rnaz/out/all.$j$i.prob >> output/rnaz.out.$j    
		paste dynalign/ct/all.$j$i.out single/ct/all.$j$i.out |sed 's/[>]//g' >> output/dynalign.energy.$j   
#		paste svmz/all.$j$i.zscore svmz/all.$j$i.sci svmz/all.$j$i.svmz.prob |sed 's/[>]//g' >> output/svmz.out.$j   
fi
done
done	
####################################
#Section: VI. Analyze the result 
####################################

#Topic: 1. annalyze the result
#option=-a
###################################
if [ "$1" = "-a" ]; then
	for j in both;do
	#generate prob file for sensitivity calculation
#paste rnaz.out.$j  all.loc.$j|cut -f 2,8,9,10,11 > rnaz.prob.$j
#	paste svmz.out.$j  all.loc.$j  | cut -f 5,11,12,13,14 > svmz.prob.$j
	#generage igb files
	awk '{if($1>=0)print "chr"$2"\t"$3"\t"$4"\tDynRNA\t"$1*1000"\t"$5}' svmz.prob.$j > igb/svmz.bed
	awk '{if($1>=0)print "chr"$2"\t"$3"\t"$4"\tRNAz\t"$1*1000"\t"$5}' rnaz.prob.$j > igb/rnaz.bed
	done
	#cut -f 1,3,4,5,7 ~/Projects/array-seq/gff/RNA.gff.merged | genome_worm.pl -stat svmz.prob.sense >&log

#Topic: 2. annotate the bins (overlap large than 20nt)
#option=-annot
###################################
elif [ "$1" = -annot ]; then
	WIG_DIR=/home1/zl222/Projects/modENCODE/worm/tilingarray
	SEG_DIR=/home1/zl222/Projects/ncRNA/Worm/matrix
	GFF_DIR=/home1/zl222/Projects/gff/worm
	for i in silverRNA 
	do	
		gff=$GFF_DIR/$i.gff
	    coor=$GFF_DIR/$i.coor  
		awk '{print "chr"$1"\t"$4"\t"$5}' $gff > $coor
   		BinOverlapper --noHeader --minOverlap=20 $SEG_DIR/dyn.plus.bed $coor | genome_worm.pl -annot $gff | awk '{print $0"\t+"}' > $i"set.bed"
   		BinOverlapper --noHeader --minOverlap=20 $SEG_DIR/dyn.minus.bed $coor  | genome_worm.pl -annot $gff | awk '{print $0"\t-"}' >> $i"set.bed"
	done
fi


