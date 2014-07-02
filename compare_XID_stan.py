import numpy as np
import pickle
import pylab as plt
from astropy.io import fits
from scipy.io import readsav

# define a function to get percentile for a particular parameter
def quantileGet(q,param):
    #q is quantile
    #param is array (nsamples,nparameters)
    # make a list to store the quantiles
    quants = []
 
    # for every predicted value
    for i in range(param.shape[1]):
        # make a vector to store the predictions from each chain
        val = []
 
        # next go down the rows and store the values
        for j in range(param.shape[0]):
            val.append(param[j,i])
 
        # return the quantile for the predictions.
        quants.append(np.percentile(val, q))
 
    return quants




file='eg5'

fit=pickle.load(open('XIDfit_noise'+file+'.p','r'))
hdulist=fits.open(file+'_src.fits')
src_f=hdulist[1].data['F']
median_fits=np.empty((src_f.size))
low_lim=quantileGet(2.5,fit['src_f'])
up_lim=quantileGet(97.5,fit['src_f'])

print src_f.size,fit['src_f'][0,:].shape
for i in range(0,src_f.size):
    median_fits[i]=np.median(fit['src_f'][:,i])

for i in range(0,src_f.size):
    print src_f[i],low_lim[i],median_fits[i],up_lim[i]
plt.errorbar(src_f,median_fits,yerr=[median_fits-low_lim,up_lim-median_fits],fmt='o')
plt.yscale('log')
plt.xscale('log')
#plt.xlim(0,20)
#plt.ylim(0,20)
plt.show()
s=readsav(file+'.sav')
Nsig=s.sig_map.flatten()
plt.hist(Nsig)
plt.show()

