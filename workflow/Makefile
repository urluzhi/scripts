# make file for scanning genome for identifying ncRNA

#specific definition here
#######
HOME=/home1/zl222/Projects/ncRNA/Worm
ALIGN1=briggsae
ALIGN2=elegans
# two sets of mouse genome in the alignments , using  mouse not mouse_leur
TAG=${CHR}

#general definition
#######
SVMZ_MODEL=../SVMZ_MODEL/model
SVMZ_SCALER=../SVMZ_MODEL/scaler
AVE_MODEL=../AVE_SD_SVM_MODEL/scaledAverageOutput-all.model
AVE_SCALER=../AVE_SD_SVM_MODEL/averageOutput-all.scaler
SD_MODEL=../AVE_SD_SVM_MODEL/scaledSDOutput-all.model
SD_SCALER=../AVE_SD_SVM_MODEL/sdOutput-all.scaler

DIR=  sequences dynalign dynalign/ct dynalign/jobs splits pairs zscore svm4ave_sd SVM single single/jobs single/ct rnaz rnaz/jobs rnaz/pairs rnaz/out rnaz/maf output
ALIGN= alignments/${TAG}.align
ALIGN_BOTH= alignments/${TAG}_both.align
REV= alignments/${TAG}_rev.align
SPLIT= splits/${TAG}.split
REV_SPLIT= splits/${TAG}_rev.split
PAIR= pairs/${TAG}.pair
SEQ= sequences/${TAG}.seq
ENERGY= dynalign/${TAG}.energy
ZSCORE= zscore/${TAG}.zscore
SCI = single/${TAG}.sci
SVM = SVM/${TAG}.svm
AVE_SD= svm4ave_sd/${TAG}.ave_sd
DYN_JOBS_NUM= 1000
SINGLE_JOBS_NUM= 1000
RNAZ_JOBS_NUM= 1000

#####################
all: 
	@echo 'Type make [options]'
tag:
	@echo "The chromosome is ${TAG}"
	
#peparing the directory structure
#####################
dir:
	for i in ${DIR};do if [ -d $$i ] ;then echo -n ""; else echo "making $$i";mkdir $$i; fi; done 

#prepare windows for dynalyn
######################
#reverse alignment
rev: ${ALIGN}
	@echo "make rev:"
	genome_encode.pl -revmum $<  ${REV}	
	cat ${ALIGN} ${REV} > ${ALIGN_BOTH}
#split alignment into windows
split: ${ALIGN} ${REV}
	@echo "make split:"
	genome_encode.pl -spmum ${ALIGN} ${SPLIT} > ${SPLIT}.log
	genome_encode.pl -spmum ${REV} ${REV_SPLIT} > ${REV_SPLIT}.log
#making paired sequences
seq: ${SPLIT} ${REV_SPLIT}
	@echo "make seq:"
	genome_encode.pl -seq ${SPLIT} ${PAIR}.sense ${SEQ}.sense  > ${PAIR}.log
	genome_encode.pl -seq ${REV_SPLIT} ${PAIR}.rev ${SEQ}.rev  > ${PAIR}_rev.log
	cat ${PAIR}.sense ${PAIR}.rev > ${PAIR}
#make dynalign jobs
dyn: ${PAIR}
	@echo "make dyn:"
	genome_encode.pl -dyn ${PAIR}  ${DYN_JOBS_NUM} dynalign/jobs/${TAG} ${HOME}	
#make sinlge folding jobs
fold: ${PAIR}
	@echo "make single fold:"
	genome_encode.pl -single ${PAIR}  ${SINGLE_JOBS_NUM} single/jobs/${TAG} ${HOME}	


#zscore calculation
######################
#svminput for ave and sd
svm-ave-sd: ${PAIR} 
	@echo "make svm-ave-sd:"
	genome_encode.pl -svminput ${PAIR} > svm4ave_sd/${TAG}.ave.svminput
	cp svm4ave_sd/${TAG}.ave.svminput svm4ave_sd/${TAG}.sd.svminput
#calculate Ave. and sd. 
ave_sd:  svm4ave_sd/${TAG}.ave.svminput  svm4ave_sd/${TAG}.sd.svminput
	@echo "make ave_sd:"
	svm-scale -r ${AVE_SCALER} svm4ave_sd/${TAG}.ave.svminput > svm4ave_sd/${TAG}.ave.svminput.scale 
	svm-predict svm4ave_sd/${TAG}.ave.svminput.scale  ${AVE_MODEL}	svm4ave_sd/${TAG}.ave.predict > svm4ave_sd/${TAG}.ave.log	
	svm-scale -r ${SD_SCALER}  svm4ave_sd/${TAG}.sd.svminput > svm4ave_sd/${TAG}.sd.svminput.scale
	svm-predict svm4ave_sd/${TAG}.sd.svminput.scale ${SD_MODEL} svm4ave_sd/${TAG}.sd.predict > svm4ave_sd/${TAG}.sd.log	
	paste  svm4ave_sd/${TAG}.ave.predict svm4ave_sd/${TAG}.sd.predict > ${AVE_SD} 
#fetch free energy from ct files of dynalign
energy: ${PAIR} 
	@echo "make energy:"
	genome_encode.pl -e ${PAIR} > ${ENERGY}
#calculate zscore
Zscore: ${ENERGY} ${AVE_SD} 
	@echo "make Zscore:"
	genome_encode.pl -z ${AVE_SD} ${ENERGY} > ${ZSCORE}


#SVM prediction
###########################
#read energy sum from ct files of single folding
energy1: ${PAIR}
	@echo "make energy1:"
	genome_encode.pl -e1 ${PAIR} > single/${TAG}.singleenergy
#calculate sci
sci: ${ENERGY} single/${TAG}.singleenergy
	@echo "make SCI:"
	paste ${ENERGY} single/${TAG}.singleenergy |awk '{if($$2==0){print "0";}else {print $$1/$$2 } }' > ${SCI}
#make SVM inputs
svm: ${ZSCORE} ${SCI} 
	@echo "make svm:"
	cut -d " " -f 2-9 svm4ave_sd/${TAG}.ave.svminput | awk '{print "1 "$$0}' > ${SVM}.tmp1
	paste ${ZSCORE} ${SCI} |awk '{print " 9:"$$1" 10:"$$2}' > ${SVM}.tmp2
	paste ${SVM}.tmp1 ${SVM}.tmp2 > ${SVM} 
	rm ${SVM}.tmp1; rm ${SVM}.tmp2
# predict probility of being ncRNA from SVM
predict: ${SVM}
	@echo "make predict:"
	svm-scale -r ${SVMZ_SCALER} ${SVM} > ${SVM}.scale
	svm-predict -b 1 ${SVM}.scale ${SVMZ_MODEL} ${SVM}.predict > ${SVM}.log
	genome_encode.pl -readsvm ${SVM}.predict SVM/${TAG}.prob
#localize the position in genome
loc: ${ZSCORE} ${ALIGN} SVM/${TAG}.prob
	@echo "make loc:"
	genome_encode.pl -loc ${PAIR} ${ZSCORE} SVM/${TAG}.prob ${ALIGN_BOTH} > output/${TAG}.loc


# rnaz
###################
maf4rnaz: ${SPLIT}
	@echo "make maf4rnaz:"
	genome_encode.pl -maf ${SPLIT} ${TAG} rnaz/maf rnaz/pairs/${TAG}.pair.sense > rnaz/pairs/${TAG}.log
	genome_encode.pl -maf ${REV_SPLIT} ${TAG} rnaz/maf rnaz/pairs/${TAG}.pair.rev >> rnaz/pairs/${TAG}.log
	cat rnaz/pairs/${TAG}.pair.sense  rnaz/pairs/${TAG}.pair.rev > rnaz/pairs/${TAG}.pair 
RNAZ: ${PAIR} 
	@echo "make rnaz:"
	genome_encode.pl -rnaz ${PAIR}  ${RNAZ_JOBS_NUM} rnaz/jobs/${TAG} ${HOME}
RNAZ_OUT: ${PAIR}
	genome_encode.pl -rnaz_out ${PAIR}  > rnaz/${TAG}.prob

