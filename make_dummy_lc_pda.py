from sim_chrom import SimChrom,IsoRP,IsoGLang,IsoBET,IsoLang
import matplotlib.pyplot as plt
import numpy as np 

#Make simulalated chromatgarm

inject=1

C0=100.1  #satulation 
W1=2
rt=1.0 
a=rt/W1*C0
iso=IsoBET(a,W1,C0,eps=1e-3,alpha=1e-2)

Ws=2.0
a=rt/Ws
iso=IsoLang(a,Ws)


chrom=SimChrom(100,iso,dead_volume=0,pass_rate=0,inject=inject)
chrom=np.array([ chrom.step() for i in range(300)])
chromB=chrom
chromA=np.roll(chrom,-60)

c=np.linspace(0,1,100)


fig=plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
plt.title("Isotherm ")
plt.plot(c,iso.iso(c))

ax2 = fig.add_subplot(2, 1, 2)
plt.plot(chromB,label="chrom")

#plt.yscale("log")
plt.legend()
plt.pause(0.1)
plt.show()
exit(0)
#Make dummpy spectrum
x=np.linspace(0,1,512)
spA=np.exp(-2*x)*(2+np.sin(x*15))
spB=spA+ np.exp(-2*x)*(np.cos(x*20)*0.3)
plt.figure()
plt.plot(spA,label="spA")
plt.plot(spB,label="spB")
plt.legend()
plt.show()

#save simulated chromagram for pyMCR's D
dA=chromA.reshape(-1,1)*spA.reshape(1,-1)
dB=chromB.reshape(-1,1)*spB.reshape(1,-1)

np.savez("dAB",dA=dA,dB=dB)