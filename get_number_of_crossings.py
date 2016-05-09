#!/usr/bin/env python

Usage="""get number of rsite-rsite pairs within a given set of log distance bins

args: rsites.dat distance_bins_per_decade outfile.dat

rsites.dat should have its second column be a sorted list of rsite locations
distance bins will be log-spaced according to distance_bins_per_decade
outfile.dat is decade start and number of crossings
"""

import sys,os
import numpy as np

try:
    rsites_file=sys.argv[1]
    nbins_per_decade=int(sys.argv[2])
    assert nbins_per_decade>0
    outfile=sys.argv[3]
except:
    sys.exit(Usage)


#read rsites list
rsites=[i.split()[1] for i in open(rsites_file)]
try:
    int(rsites[0])
except:
    rsites=rsites[1:]
rsites=map(int,rsites)
rsites=np.array(rsites)

#build bins
stepsz=1./nbins_per_decade
dbins=[]
upper=np.log10(max(rsites)-min(rsites))+stepsz
cut=0.
while cut<=upper:
    dbins.append(10**cut)
    cut+=stepsz
dbins=np.array(dbins)

#loop over rsites
counts=np.zeros(len(dbins)-1, dtype=int)
for i in xrange(len(rsites)-1):
    print "{}\r".format(i),
    counts += np.histogram(rsites[(i+1):]-rsites[i], bins=dbins)[0]

fl=open(outfile,'w')
for i in xrange(len(counts)):
    fl.write("{}\t{}\n".format(dbins[i],counts[i]))

