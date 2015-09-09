#! /usr/bin/python
# stats.py
# P. Clote

import sys,os,tempfile,string,math


def getStats(L):
  #L is list
  if L==[]:
    return (0,0,0,0)
  max = -sys.maxint
  min = sys.maxint
  
  count = 0
  sum   = 0.0
  sumSquares = 0.0
  for value in L:
    x = float(value)
    if x<min: min=x
    if x>max: max=x
    count = count+1
    sum = sum+x
    sumSquares = sumSquares+x**2
  mean = sum/count
  variance = sumSquares/count - mean**2
  stdev    = math.sqrt(variance)
  #print "Mean:%f\tStDev:%f\tMax:%f\tMin:%f" % (mean,stdev,max,min)
  return (mean,stdev,max,min)


if __name__ == '__main__':
  L = sys.argv[1:] 
  (mean,stdev,max,min) = getStats(L)
  print "Mean:%f\tStDev:%f\tMax:%f\tMin:%f" % (mean,stdev,max,min)

