import matplotlib.pyplot as plt
import numpy as np
from constraint import ConstraintBEMG,RegArray,FunBEMG,FunLin
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
dB*=0.4
D=dA+dB +  np.random.randn(*dA.shape)*0.01
gA=np.sum(dA,axis=1)#Grand Truth
gB=np.sum(dB,axis=1)#Grand Truth
gsA=np.sum(dA,axis=0)#Grand Truth
gsB=np.sum(dB,axis=0)#Grand Truth


from sklearn.decomposition import FastICA,NMF
ica=FastICA(n_components=2)
#ica=NMF(n_components=2)
DN=np.copy(D)
DN[DN<1e-100]=1e-100
C_est=ica.fit_transform(DN)
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




if True:
    #BEMG関数制約にregressorを使う
    #処理がとっても重いが、収束はする　　しかし、最適解にいくとは限らない
    c_regr=RegArray([FunBEMG(),FunBEMG()])
    mcrar_c = McrAR(c_regr=c_regr, st_regr='NNLS',st_constraints=[],c_constraints=[])
    #C_est_c=ConstraintBEMG(range(2)).transform(C_est)
    mcrar_c.fit(D, C=C_est,verbose=True)

    plt.figure()
    plt.title("result by MCR regresssion with BEMG regressoijn")
    print(np.sum(mcrar.ST_opt_,axis=1).shape)

    plt.plot(mcrar_c.C_opt_* np.sum(mcrar_c.ST_opt_,axis=1).reshape(1,-1))
    plt.plot(gA,label="trueA")
    plt.plot(gB,label="trueB")
    plt.legend()
    plt.pause(0.1)
    #plt.show()

    plt.figure()
    plt.title("result by MCR with BEMG regression / Spectrum")
    plt.plot(mcrar_c.ST_opt_.transpose()* np.sum(mcrar_c.C_opt_,axis=0))
    plt.plot(gsA,label="trueA")
    plt.plot(gsB,label="trueB")
    plt.legend()
    plt.pause(0.1)

if True:
    #Constrintを使うほうが処理が早いが、収束が保証されないので一旦誤差量が上がる
    mcrar_c = McrAR(c_regr='NNLS', st_regr='NNLS',st_constraints=[],c_constraints=[ConstraintBEMG(range(2)) ],tol_increase=100)
    #C_est_c=ConstraintBEMG(range(2)).transform(C_est)
    mcrar_c.fit(D, C=C_est,verbose=True)

    plt.figure()
    plt.title("result by MCR with BEMG constraint")
    print(np.sum(mcrar.ST_opt_,axis=1).shape)

    plt.plot(mcrar_c.C_opt_* np.sum(mcrar_c.ST_opt_,axis=1).reshape(1,-1))
    plt.plot(gA,label="trueA")
    plt.plot(gB,label="trueB")
    plt.legend()
    plt.pause(0.1)
    #plt.show()

    plt.figure()
    plt.title("result by MCR with BEMG constraint/ Spectrum")
    plt.plot(mcrar_c.ST_opt_.transpose()* np.sum(mcrar_c.C_opt_,axis=0))
    plt.plot(gsA,label="trueA")
    plt.plot(gsB,label="trueB")
    plt.legend()
    plt.pause(0.1)

plt.show()
