#! /usr/bin/env python

#shuffleSequenceUsingAltschulErikson.py
#P. Clote, Oct 2003

#------------------------------------------------------------------
#Input RNAs in FASTA file, and compute NUM many shufflings of RNA sequence
#using Altman-Erikson randomly shuffled dinucleotide method.
#------------------------------------------------------------------

PRINT   = 0
LINELEN = 70

import sys,os,stats,string
from altschulEriksonDinuclShuffle import dinuclShuffle
import computeRNAfoldEnergyForRNAsInFile




def file2string(fileName):
  file = open(fileName,"r")
  L = []
  line = file.readline()
  while line:
    while line[0]==">":  # treat lines beginning with '>' as comment and skip 
      line = file.readline()
      continue
    else: 
      line = line[:-1]
      L.append(line)
      line = file.readline()
  text = string.join(L,"")
  return text


def main(fileName,NUM):
  seq = file2string(fileName)
  for i in range(NUM):
    shuffledSeq = dinuclShuffle(seq) 
    sys.stdout.write(">%d\n" % (i+1))
    sys.stdout.write("%s\n" % shuffledSeq)



  
if __name__ == '__main__':  
  if len(sys.argv) < 3 :
     print "Usage: %s RNAs.faa NUM" %  sys.argv[0]
     text = """
            1) RNA.faa is FASTA file of ONE RNA sequence
            2) NUM is number of random sequences to generate by
               shuffling the dinucleotides of RNAs input
     Script to compute Altman-Erikson randomly shuffled dinucleotides.
            """
     print text
     sys.exit(1)
  fileName = sys.argv[1]
  NUM      = int(sys.argv[2])
  main(fileName,NUM)
  
