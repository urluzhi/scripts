#!/bin/bash

####################################
#Header: run_tf_worm_mab5.sh
#Running Transcription factor mab-5
####################################


####################################
#Files: File requirement
#NA
####################################

####################################
#Variables: 
DIR=mab5_L3_41608_42308-mab5_L3_GFP-0206091753-results
####################################


####################################
#Section: I. Running ashish's program
#finding targets from Guoneng's peakseq hits
###################################

#Topic: 1.Finding targets of mab-5
#option=-t
###################################
if [ "$1" = "-t" ]; then
	for i in 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1;do
		~aa544/bin/findTargets -i $DIR/all/hits/normalized_filtered_hits.bed -o $DIR/targets/pval$i --pval=$i  
	done
#Topic: 2. make dir
#option=-mkdir	
###################################
elif [ "$1" = "-mkdir" ] ;then
	echo "ok"
else
	echo "type something"
fi

