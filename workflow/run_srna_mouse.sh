#!/bin/bash
if [ $# = 0 ]; then echo "Run as: ./run.sh #";exit;fi

dir=( uplow all ) 
#############################################################
#I. MAPPING(ELAND & Bowtie)
#############################################################
# 1. clean the reads
#############################################################
if [ $1 = 1 ]; then
	awk 'BEGIN {i=0}{ print $0"\ts"i;i++ }' GSE10364.reads > GSE10364.name
	awk -F '\t' '{if($2)  print ">"$5"_x"$2"\n"$1}' GSE10364.name > low.fa
	awk -F '\t' '{if($3)  print ">"$5"_x"$3"\n"$1}' GSE10364.name > up.fa
	awk -F '\t' '{if($4)  print ">"$5"_x"$4"\n"$1}' GSE10364.name > third.fa
	awk -F '\t' '{if($2||$3)  print ">"$5"_x"$2+$3"\n"$1}' GSE10364.name > uplow.fa
	awk -F '\t' '{if($2||$3||$4)  print ">"$5"_x"$2+$3+$4"\n"$1}' GSE10364.name > all.fa
#############################################################
# 2. split by length
#############################################################
elif [ $1 = 2 ]; then
		splitByLength.py uplow.fa > uplow.log
		splitByLength.py all.fa >all.log
#############################################################
# 3. Eland: pbs jobs, 1 or 25 hits maximum
#############################################################
elif [ $1 = 3 ]; then
		g=/home1/zl222/Projects/genomes/mouse/ensemble_mm9/squash/
   		cd /home1/zl222/Projects/sRNApgene/sRNAmouse/eland	
		i=uplow1 
			echo "#!/bin/bash"  > $i.pbs
			echo "cd $(pwd)"  >> $i.pbs
			for f in $(find . -type f -name 'uplow.fa_*');	do
				bn=$(basename $f)
				ln=${bn/uplow.fa_}
				lln=$(( $ln > 32 ? 32 : $ln ))
	#echo "eland_$lln $bn $g $i"_"$ln.out --multi=25,25,25 >& $i"_"$ln.log" >> $i.pbs
	echo "eland_$lln $bn $g $i"_"$ln.out  >& $i"_"$ln.log" >>$i.pbs
		done
#############################################################
# 4. bowtie: pbs jobs, 1 or 5 maximum 
#############################################################
elif [ $1 = 4 ]; then
		cd /home1/zl222/Projects/sRNApgene/sRNAmouse/bowtie	
		g=/home1/zl222/Projects/genomes/mouse/ensemble_mm9/bowtie/Mus_musculus.NCBIM37.49.dna.chromosome
		for i in uplow.dup  all.dup ;do 
			for n in 1 5;do
				echo "#!/bin/bash"  > $i.$n.pbs
				echo "#PBS -l nodes=1:ncpus=4" >>$i.$n.pbs
				echo "cd $(pwd)"  >> $i.$n.pbs
				echo "bowtie -p 4 -m $n -k $n -f $g ../reads/$i.fa $i.$n.map" >>$i.$n.pbs 
				echo "awk  '{if(\$2~/+/) print \$0}' $i.$n.map > $i.$n.plus.map" >>$i.$n.pbs
				echo "awk  '{if(\$2~/-/) print \$0}' $i.$n.map > $i.$n.minus.map" >>$i.$n.pbs
			done
		done

fi
#############################################################
#############################################################
# II. MERGE, CLUSTER and FEATURE
#############################################################
ELAND=17-32
#############################################################
# 5. collect the result
#############################################################
if [ $1 = 5 ]; then
	wc -l $ELAND.out 
#############################################################
# 6. extract semi-unique mapping (up to 2snp) locations/coordinates 
#############################################################
elif [ $1 = 6 ]; then
	for i in ${dir[@]} ; do 
		id=$ELAND.$i
		grep  -P "\t1:\d+:\d+\t" $id.out  > $id.u0
		grep  -P "\t0:1:\d+\t" $id.out > $id.u1
		grep  -P "\t0:0:1\t" $id.out > $id.u2

		mapsrna.pl -list_locations $id.u0 0 > $id.u0.loc 
		mapsrna.pl -list_locations $id.u1 1 > $id.u1.loc 
		mapsrna.pl -list_locations $id.u2 2 > $id.u2.loc 
		cat $id.u0.loc  $id.u1.loc > $id.u01.loc
		cat $id.u01.loc $id.u2.loc > $id.u012.loc
	done
#############################################################
# 7.  merge the mapped locations
#############################################################
elif [ $1 = 7 ]; then
for i in ${dir[@]} ; do 
	id=$ELAND.$i
	for j in 0 01 012;do
		#merge the reads
		grep "+" $id.u$j.loc > $id.u$j.loc.plus
		grep "-" $id.u$j.loc > $id.u$j.loc.minus
		sort -k1,1  -k2,2n -k3,3n $id.u$j.loc.plus | mapsrna.pl -inclu > $id.u$j.loc.plus.inclu 
		sort -k1,1  -k2,2n -k3,3n $id.u$j.loc.minus | mapsrna.pl -inclu > $id.u$j.loc.minus.inclu 
		sort -k1,1  -k2,2n -k3,3n $id.u$j.loc.plus | mapsrna.pl -merge 10 > $id.u$j.loc.plus.mer10 
		sort -k1,1  -k2,2n -k3,3n $id.u$j.loc.minus | mapsrna.pl -merge 10 > $id.u$j.loc.minus.mer10 
	done
done
#############################################################
# 8.  prepare pgene
#############################################################
elif [ $1 = 8 ]; then
	awk '{if(NR>1 && $4~/+/)print $0}'  Pgenes_without_seqs.txt |sort -k1,1 -k2,2n -k3,3n > pgene.plus
	awk '{if(NR>1 && $4~/-/)print $0}'  Pgenes_without_seqs.txt |sort -k1,1 -k2,2n -k3,3n > pgene.minus
	awk '{print "Mus_musculus.NCBIM37.49.dna.chromosome."$1"\t"$2"\t"$3}' pgene.plus > pgene.plus.coor	
	awk '{print "Mus_musculus.NCBIM37.49.dna.chromosome."$1"\t"$2"\t"$3}' pgene.minus > pgene.minus.coor	
#############################################################
# 9. map to psedogene and score the density,coverage
#############################################################
elif [ $1 = 9 ]; then
	for i in ${dir[@]} ; do for j in  0 01 012;do
		id=$ELAND.$i
		for l in plus minus ; do for m in plus minus;do
			#overlaper first
			#BinOverlapper pgene.$l.coor  ../eland/$id.u$j.loc.$m.inclu -noheader |awk '{if($4) {print $0} else {print $0"NA,"} }' > $id.inclu.u$j.$l-$m	
			BinOverlapper pgene.$l.coor  ../eland/$id.u$j.loc.$m.mer10 -noheader |awk '{if($4) {print $0} else {print $0"NA,"} }' > $id.mer10.u$j.$l-$m	
			#fetch the merged sequences
			mapsrna.pl -fetch $id.mer10.u$j.$l-$m  ../eland/$id.u$j.loc.$m.mer10 > $id.mer10.u$j.$l-$m.f
			#calculate the score
			cat  $id.mer10.u$j.$l-$m.f | mapsrna.pl -score > $id.mer10.u$j.$l-$m.score
 		done;done
done; done
#############################################################
# 10.analyze result
#############################################################
elif [ $1 = 10 ]; then
for i in all; do
	id=$ELAND.$i
	for j in  012;
	do for l in plus minus; do  
		for m in plus minus;do
#awk -F '\t' '{if($0~/NA/){ print "0\t0\t0\t0\t0\t0\t0\t"$1"\t"$2"\t"$3} else{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11} }' $id.mer10.u$j.$l-$m.score > $id.mer10.u$j.$l-$m.score.s 
######################
#1.allreads 2.coverage/100nt 3.maxread/100nt 4.numofwindows (win_name) 
#5.max_density_merge 6.max_Length_mer 7.max_reads_mer 8.chr 9.start 10.end
######################
			goo=good
		done
#paste $id.mer10.u$j.$l-plus.score.s $id.mer10.u$j.$l-minus.score.s pgene.$l | awk '{if ( ($2>0&&$3>=4)||($12>0&&$13>=4) ) print $3"\t"$13"\t"$0}'   
		#paste $id.mer10.u$j.$l-plus.score $id.mer10.u$j.$l-minus.score |grep -v "NA.*NA"  
	done
	paste $id.mer10.u$j.plus-plus.score.s $id.mer10.u$j.plus-minus.score.s pgene.plus | awk '{if ( ($2>=50&&$3>=4)||($12>=50&&$13>=4) ) print $3"\t"$13"\t"$0}'   
	paste $id.mer10.u$j.minus-minus.score.s $id.mer10.u$j.minus-plus.score.s pgene.minus | awk '{if ( ($2>=50&&$3>=4)||($12>=50&&$13>=4) ) print $3"\t"$13"\t"$0}'   


#awk '{print "chr"$1"\tPseudopipe\t"$6"\t"$2"\t"$3"\t0\t"$4"\t.\t"$7}' Pgenes_without_seqs.txt
done; done
#############################################################
# 9. send list for Nitin
#############################################################
elif [ $1 = 11 ]; then
for i in all; do
	id=$ELAND.$i
	for j in  012;do for l in plus minus; do  
		for m in plus minus;do
#paste $id.mer10.u$j.$l-$m pgene.$l |grep -v NA |cut -f 4- >  pgenelist.$l-$m 
			paste $id.mer10.u$j.$l-$m.score.s  pgene.$l |awk '{if($2>0&&$3>4) print $0}' 
		done
	done; done
done
elif [ $1 = 12 ] ;then
	for l in plus minus;do
		sed 's/Mus_musculus.NCBIM37.49.dna.chromosome.//g' pgene.$l.coor | mapsrna.pl -subseq ../genomes/bowtie/Mus_musculus.NCBIM37.49.dna.chromosome.all.fa > pgene.$l.fa
	done
fi

