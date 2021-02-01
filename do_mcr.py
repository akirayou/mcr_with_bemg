import matplotlib.pyplot as plt
import numpy as np
from constraint import ConstraintBEMG
from pymcr.constraints import ConstraintNonneg, ConstraintNorm

import logging
import sys
logger = logging.getLogger('pymcr')
logger.setLevel(logging.DEBUG)
stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_format = logging.Formatter('%(message)s')  # Just a basic message akin to print statements
stdout_handler.setFormatter(stdout_format)
logger.addHandler(stdout_handler)

d=np.load("dAB.npz")
dA=d["dA"]
dB=d["dB"]
dB*=0.5
D=dA+dB +  np.random.randn(*dA.shape)*0.001
gA=np.sum(dA,axis=1)#Grand Truth
gB=np.sum(dB,axis=1)#Grand Truth
gsA=np.sum(dA,axis=0)#Grand Truth
gsB=np.sum(dB,axis=0)#Grand Truth


from sklearn.decomposition import FastICA
ica=FastICA(n_components=2)
C_est=ica.fit_transform(D)
#flip negative peak
C_est[:,np.max(C_est,axis=0)  < -1*np.min(C_est,axis=0)]*=-1
#remove base line
C_est -= np.min(C_est,axis=0)

plt.figure()
plt.title("Initial estimate by FastICA")
plt.plot(C_est)
plt.pause(0.1)

from pymcr.mcr import McrAR
mcrar = McrAR()
mcrar.fit(D, C=C_est,verbose=True)

plt.figure()
plt.title("result by MCR")
plt.plot(mcrar.C_opt_  * np.sum(mcrar.ST_opt_,axis=1)  )
plt.plot(gA,label="trueA")
plt.plot(gB,label="trueB")
plt.legend()
plt.pause(0.1)


mcrar_c = McrAR(c_regr='NNLS', st_regr='NNLS',st_constraints=[],c_constraints=[ConstraintBEMG(range(2)) ])
C_est_c=ConstraintBEMG(range(2)).transform(C_est)

mcrar_c.fit(D, C=C_est_c,verbose=True)
#mcrar_c.fit(D, ST=mcrar.ST_opt_,verbose=True)

plt.figure()
plt.title("result by MCR with BEMG")
print(np.sum(mcrar.ST_opt_,axis=1).shape)

plt.plot(mcrar_c.C_opt_* np.sum(mcrar_c.ST_opt_,axis=1).reshape(1,-1))
plt.plot(gA,label="trueA")
plt.plot(gB,label="trueB")
plt.legend()
plt.pause(0.1)
#plt.show()

plt.figure()
plt.title("result by MCR with BEMG / Spectrum")
plt.plot(mcrar_c.ST_opt_.transpose()* np.sum(mcrar_c.C_opt_,axis=0))
plt.plot(gsA,label="trueA")
plt.plot(gsB,label="trueB")
plt.legend()
plt.show()
