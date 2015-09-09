#!/bin/bash
for i in I II III IV V X ;do 
echo $i
compare_fasta.pl /home2/yhl3/data/ensembl/caenorhabditis_elegans_50_190/dna/Caenorhabditis_elegans.WS190.50.dna.chromosome.$i".fa" /home1/zl222/Projects/sRNApgene/sRNAworm/genome/WS190/$i".fa"
done
# $Id: $


