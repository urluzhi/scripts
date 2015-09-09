#!/bin/bash
sed 's/,/\t/g' |sort -k4,4 -k1,1 -k2,2n -k3,3n |mergeloc.pl | awk '{sum+=$3-$2+1}END{print sum}' 


