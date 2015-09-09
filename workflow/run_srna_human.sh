#!/bin/bash
if [ $# = 0 ]; then echo "Run as: ./run.sh #";exit;fi

#############################################################
# Topic: 1. clean the maps
#############################################################
if [ "$1" = "-c" ]; then
	#sort -t ','  -k5,5 -k2,2 -k3,3n -k4,4n records.txt > records.sorted 
	#awk -F ','  '{sum=0;for(i=6;i<=42;i++) {sum+=$i} print $2"\t"$3"\t"$4"\ts"NR"\t"sum"\t"$5}' records.sorted  > records.bed
#grep '+' records.bed |awk '{if($5>0) print $0}'> records.plus.bed
	grep '-' records.bed |awk '{if($5>0) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6}'> records.minus.bed
	
#############################################################
# Topic: 2.  merge the mapped locations
#############################################################
elif [ "$1" = "-m" ]; then
	id=records
	for i in plus minus;do 
		sort -k1,1  -k2,2n -k3,3n $id.$i.bed |awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"_x"$5}' | mapsrna.pl -merge 10 > $id.$i.mer10 
	done
#############################################################
# Topic: 3. prepare pgene
#############################################################
elif [ "$1" = "-pgene" ]; then	
	for i in pipeline  zhang; do
		#awk -F '\t' '{if(NR>1 && $7~/+/)print $1"\t"$4"\t"$5"\t"$7}'  $i.gtf  |sort -k1,1 -k2,2n -k3,3n > $i.plus.bed
		#awk -F '\t' '{if(NR>1 && $7~/-/)print $1"\t"$4"\t"$5"\t"$7}'  $i.gtf  |sort -k1,1 -k2,2n -k3,3n > $i.minus.bed
		awk -F '\t' '{if(NR>1 && $7~/+/)print $1"\t"$4"\t"$5"\t"$7"\t"$9}'  $i.gtf  |sort -k1,1 -k2,2n -k3,3n  > $i.plus.gtf
			
		awk -F '\t' '{if(NR>1 && $7~/-/)print $1"\t"$4"\t"$5"\t"$7"\t"$9}'  $i.gtf  |sort -k1,1 -k2,2n -k3,3n  > $i.minus.gtf

	
		done
#############################################################
# Topic: 4. map to pgene and score the density,coverage
#############################################################
elif [ "$1" = "-map" ]; then
	id=records
	for i in pipeline zhang;do
		for l in plus minus ; do for m in plus minus;do
			#overlaper first
			BinOverlapper --noHeader pgene/$i.$l.bed reads/$id.$m.mer10  |awk '{if($4) {print $0} else {print $0"NA,"} }' > pgene/$i-$id.mer10.$l-$m	
			#fetch the merged sequences
			mapsrna.pl -fetch pgene/$i-$id.mer10.$l-$m  reads/$id.$m.mer10 > pgene/$i-$id.mer10.$l-$m.f
			#calculate the score
			cat  pgene/$i-$id.mer10.$l-$m.f | mapsrna.pl -score > pgene/$i-$id.mer10.$l-$m.score
 		done;done
done
#############################################################
# Topic: 5.analyze result
#############################################################
elif [ "$1" = "-a1" ]; then
	id=records
	for i in pipeline zhang;do
		for l in plus minus; do  
		for m in plus minus;do
		awk -F '\t' '{if($0~/NA/){ print "0\t0\t0\t0\t0\t0\t0\t"$1"\t"$2"\t"$3} else{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11} }' $i-$id.mer10.$l-$m.score > $i-$id.mer10.$l-$m.score.s 
######################
#1.allsequences 2.coverage/100nt 3.maxread/100nt 4.numofwindows (win_name) 
#5.max_density_merge 6.max_Length_mer 7.max_reads_mer 8.chr 9.start 10.end
######################
		done
#paste $id.mer10.u$j.$l-plus.score.s $id.mer10.u$j.$l-minus.score.s pgene.$l | awk '{if ( ($2>0&&$3>=4)||($12>0&&$13>=4) ) print $3"\t"$13"\t"$0}'   
		#paste $id.mer10.u$j.$l-plus.score $id.mer10.u$j.$l-minus.score |grep -v "NA.*NA"  
		done
	done
elif [ "$1" = "-a2" ]; then
	id=records
	for i in zhang  ;do
	paste $i-$id.mer10.plus-plus.score.s $i-$id.mer10.plus-minus.score.s $i.plus.gtf| awk '{if ( $3>=4 || $13 >=4 ) print $0}'   
	paste $i-$id.mer10.minus-minus.score.s $i-$id.mer10.minus-plus.score.s $i.minus.gtf| awk '{if ( $3>=4 || $13 >=4 ) print $0}'   

#awk '{print "chr"$1"\tPseudopipe\t"$6"\t"$2"\t"$3"\t0\t"$4"\t.\t"$7}' Pgenes_without_seqs.txt
	done 
#############################################################
# 9. send list for Nitin
#############################################################
elif [ $1 = 11 ]; then
	cut -f 25 pgene13.out |cut -d ";" -f 4 |cut -d " " -f 3 |cut -d "." -f 1 > genename 
	id=records
	for i in pipeline ;do
	paste $i-$id.mer10.plus-plus.score.s $i-$id.mer10.plus-minus.score.s $i.plus.bed | awk '{if ( $3>=4 || $13 >=4 ) print $0}'   
	paste $i-$id.mer10.minus-minus.score.s $i-$id.mer10.minus-plus.score.s $i.minus.bed| awk '{if ( $3>=4 || $13 >=4 ) print $0}'   

done
elif [ $1 = 12 ] ;then
	for l in plus minus;do
		sed 's/Mus_musculus.NCBIM37.49.dna.chromosome.//g' pgene.$l.coor | mapsrna.pl -subseq ../genomes/bowtie/Mus_musculus.NCBIM37.49.dna.chromosome.all.fa > pgene.$l.fa
	done
fi

