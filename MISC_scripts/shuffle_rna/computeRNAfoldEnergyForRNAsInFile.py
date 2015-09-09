#! /usr/bin/env python

# computeRNAfoldEnergyForRNAsInFile.py
# P. Clote, Sep 2003


RNAfold = "/usr/src/ViennaRNA-1.4/Progs/RNAfold "
import sys,os
PRINT = 0



def main(fileName):
  file = open(fileName)
  line = file.readline()
  L = []
  while line:
    tmpfile = open("tmp0","w")
    FASTAcomment = line
    tmpfile.write(FASTAcomment)
    tRNA = file.readline()
    tmpfile.write(tRNA)
    tmpfile.close()
    cmd  = RNAfold+" < tmp0"
    foo  = os.popen(cmd,"r")
    line = foo.readline()
    while line:
      oldLine = line
      line = foo.readline()
    foo.close()
    words = oldLine.split()
    x = words[-1]
    if x[-1]==")": x = x[:-1]
    if x[0] =="(": x = x[1:]
    x = float(x)
    if PRINT: print "%f" %x 
    L.append(x)
    line = file.readline() 
  return L


if __name__ == '__main__':
  if len(sys.argv) != 2:
    print "Usage: %s file" % sys.argv[0]
    sys.exit(1)
  PRINT = 1
  main(sys.argv[1]) 
