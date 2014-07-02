import numpy as np
import scipy
import astropy.io.ascii as ascii
from scipy.sparse import coo_matrix, dia_matrix
from scipy.io import readsav
from scipy.sparse import linalg

#---File name-----
file='eg5'
s=readsav(file+'.sav')
#-- create data array
db=s.map.flatten()
npix=db.shape[0]
#-- create diag sig array
Nsig=s.sig_map.flatten()
Nsig[:]=0.00001
#--Read image matrix
Asp=np.loadtxt(file+'.matrix')
A=coo_matrix((Asp[:,2],(Asp[:,0],Asp[:,1])),shape=(db.shape[0],100))
Amat=A.todense()

import pystan
print np.max(Asp[:,1]),np.max(Asp[:,0])
XID_data={'npix':npix,
          'nsrc':100,
          'nnz':Asp[:,2].size,
          'db':db,
          'sigma':Nsig,
          'Val':Asp[:,2],
          'Row': Asp[:,0].astype(long),
          'Col': Asp[:,1].astype(long)}

fit=pystan.stan(file='XIDfit.stan',data=XID_data,iter=1000,chains=4,verbose=True)
import pickle
pickle.dump(fit.extract(),open('XIDfit_nnoise'+file+'.p','wb'))




