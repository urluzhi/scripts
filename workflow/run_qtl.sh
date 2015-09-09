#!/bin/bash
if [ $# = 0 ]; then echo "Run as: ./run.sh #";exit;fi

#############################################################
# Topic: 1. filter the unqualified terms in qtl
#############################################################
if [ "$1" = "-f" ]; then
	i=qtl_human
	awk -F '\t' '{if(NR>1 &&$5>0&&$6>0) print $0}' $i.txt > $i.filter.txt
	cut -f 4-6 $i.filter.txt |sort -k1,1 -k2,2n -k3,3n > $i.filter.bed
#############################################################
# Topic: 2. map pseudogene with 
#############################################################
	

fi

