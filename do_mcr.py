import matplotlib.pyplot as plt
import numpy as np
from constraint import ConstraintBEMG
from pymcr.constraints import ConstraintNonneg, ConstraintNorm
d=np.load("dAB.npz")
dA=d["dA"]
dB=d["dB"]
dB*=1
D=dA+dB + np.random.randn(*dA.shape)*0.001
gA=np.sum(dA,axis=1)#Grand Truth
gB=np.sum(dB,axis=1)#Grand Truth


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


mcrar_c = McrAR(c_constraints=[ConstraintNonneg(),ConstraintBEMG(range(2)) ])
mcrar_c.fit(D, C=C_est,verbose=True)
#mcrar_c.fit(D, ST=mcrar.ST_opt_,verbose=True)

plt.figure()
plt.title("result by MCR with BEMG")
plt.plot(mcrar_c.C_opt_* np.sum(mcrar.ST_opt_,axis=1))
plt.plot(gA,label="trueA")
plt.plot(gB,label="trueB")
plt.legend()
plt.show()
