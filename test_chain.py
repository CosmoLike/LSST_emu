import numpy as np
from numpy.linalg import inv
from numpy.linalg import det
import time
import matplotlib.pyplot as plt
from scipy.stats import chi2
### setting parameters here

# where you save your data
datadir="./"
# the filename of your list file, you can omit some files
# by adding '#' at the beginning of the line
fnames_list = "chain_list"
# number of bins: the total chain is breaked into num_subchain
# bins, then calculate FoM of each bin
num_subchain = 15
n_x,n_y = 5,5
dangerous = set([1,5,13])
# if you have already computed the FoM, just read existing
# results, do not need compute again.
READ_FILE = False

### Main program
# read filenames of the chains
fnames = []
with open(fnames_list) as fp:
    for line in fp.readlines():
        if line[0]!='#':
            line = line.strip()
            fnames.append(line)
Nchains = len(fnames)
print("number of chains: ",Nchains)
print fnames

for i,fname in enumerate(fnames):
    if not READ_FILE:
        print("Reading %s"%fname)
        t1 = time.time()
        data = np.genfromtxt(datadir+fname, skip_header=0,skip_footer=0,names=True, dtype=None)
        t2 = time.time()
        # merge sort
        data_ = np.sort(data, order='log_like',kind='mergesort')
        t3 = time.time()
        print("Length of the file = %d\n "
              "Time used in reading data = %d s\n "
              "Time used in sorting data = %d s"%(len(data), 
                                                  t2-t1, t3-t2))
        # clip outlier
        for j in range(len(data)):
            if np.isfinite(data_['log_like'][j]):
                startid = j
                break
        len_68 = int((len(data_)-startid-1)*0.6827)

        chi2_thresh = data_['log_like'][startid+len_68]
        print("log_likelihood threshold is %.3f (p-value = %.3f)"%(-2*chi2_thresh, 1-chi2.cdf(-chi2_thresh*2,len(data_.dtype)-1)))
        # clip out terms with too large chi2
        data_DE = np.zeros([len_68,2])
        n, j = 0, 0
        while n<len_68:
            if np.isfinite(data['log_like'][j]) and \
            data['log_like'][j]<=chi2_thresh:
                data_DE[n,0] = data['w0'][j]
                data_DE[n,1] = data['wa'][j]
                n += 1
            j += 1
        # calculate FoM
        steps, fom = np.zeros(15), np.zeros(15)
        subchain_len = int(len_68/15)
        FoM = np.zeros(15)
        for j in range(15):
            cov=np.cov(data_DE[j*subchain_len:(j+1)*subchain_len,:],rowvar=False)
            FoM[j]=(np.power(det(inv(cov)),1./2))
            print j*subchain_len,(j+1)*subchain_len,FoM[j]
            
        steps = [(j+1)*subchain_len/1000000 for j in range(15)]

        # write results
        data_write = np.vstack([steps, FoM]).T
        np.savetxt(fname+'_cov', data_write, fmt='%.18e %.18e', delimiter=' ',
                   header='# step FoM', comments='# chi2 threshold = %f\n'%(-2*chi2_thresh))
    else:
        print("Reading existing results from %s"%(fname+"_cov"))
        data_write = np.genfromtxt(fname+"_cov", names=True, skip_header=1, dtype=None)
        steps = data_write['step']
        FoM = data_write['FoM']
