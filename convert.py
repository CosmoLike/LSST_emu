#!/usr/bin/python
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i","--ifile", help="Input File Datavector - assumed to be masked")
parser.add_argument("-o","--ofile", help="Ouput File Datavector - it will have 0.2 length")
args = parser.parse_args()

idata = np.genfromtxt(args.ifile)
xx = np.shape(idata)
start= int(0.8*xx[0])
odata = np.zeros((xx[0]-start,xx[1]+1))
odata[:,1] = idata[start:,-1]
odata[:,0] = 1
for i in range(xx[1]-1):
   odata[:,i+2] = idata[start:,i]
np.savetxt(args.ofile, odata)
