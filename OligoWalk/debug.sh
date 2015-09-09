#!/bin/bash
export DATAPATH="NN_data/"
#OligoWalk -seq example_long.seq -type dd -o nofile -m 2 -st 1 -en 0 -M 1 -N 0 -l 25 -s 0 -co 1  -unit -7  -fi 1 -fold 76 -test 4 2 5 160 4438
#OligoWalk -seq /home1/zl222/FileServerBackup/array-seq/sequences/transcripts/seq/I+11499_16828.seq -type dd -o nofile -m 2 -st 1 -en 0 -M 1 -N 0 -l 25 -s 0 -co 1  -unit -7  -fi 1 -fold 76 -test 4 2 5 160 4438
   

exe/OligoWalk -type  r -seq  example.seq -o nofile -m 2 -st 1 -en 0 -M  1 -N  0 -l  19 -s  2 -co  1 -unit  -7 -fi  0 -fold 60 -score 
   

