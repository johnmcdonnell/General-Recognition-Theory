
import os
import pylab as pl
import numpy as np

import grt

tests = [ grt.nullmodel, grt.interceptmodel, grt.onedfit, grt.twodfit ]
datafolder = 'johndata/nolabeled'
datafolder = 'data'

def load_data_raw( fn ):
    pipe = os.popen( "awk '$4==1 { print $5, \" \", $6, \" \", $10 }' %s | sed s/ch1/0/ | sed s/ch2/1/" % fn )
    data = pl.loadtxt( pipe )
    return data

datafiles = [os.path.join( datafolder, fn ) for fn in os.listdir(datafolder)
	     if fn[-4:]==".dat" ]


fits = {}
bics = {}
retdict = dict()
for df in datafiles:
    if "alldata.dat" in df:
        continue
    print "Fitting ", df
    data = load_data_raw( df )
    fits[df] = [ test(data) for test in tests ]


for df in datafiles:
    print df
    nlls = [ fit.negloglike for fit in fits[df] ]
    bics = [ fit.bic for fit in fits[df] ]
    print "-log(likelihood):", nlls
    print "bics:", bics
    winner = np.argmin( bics )
    print winner
    if winner==2:
        print "Dimension: ", fits[df][2].dim
    
    #retdict[df] = (bics,np.argmin(bics), fits[2].dim, fits)


