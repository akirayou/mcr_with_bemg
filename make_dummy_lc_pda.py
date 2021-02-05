from sim_chrom import SimChrom,IsoRP,IsoGLang
import matplotlib.pyplot as plt
import numpy as np 

#Make simulalated chromatgarm
#1/a+1/b=1
#1/b=1/(1-1/a)
b=2
Ws=1
n=1.0
m=1.2

chrom=SimChrom(100,IsoGLang(b,Ws,n,m),dead_volume=30,pass_rate=0)
chrom=np.array([ chrom.step() for i in range(300)])
chromB=chrom
chromA=np.roll(chrom,-60)

plt.figure()
#plt.plot(chromA,label="chromA")
plt.plot(chromB,label="chromB")
plt.yscale("log")
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