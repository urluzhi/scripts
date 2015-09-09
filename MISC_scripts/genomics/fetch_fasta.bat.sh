#!/bin/bash
sed 's/chr//g' $1 | awk '{cmd="~/svn/John_scripts/seqmanipulator.pl -fetch /home1/zl222/Projects/genomes/worm/elegans/fasta/elegans.WS170.dna.fa "$0;system(cmd) }'
# $Id: $


