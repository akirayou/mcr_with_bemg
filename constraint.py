
import pymcr.constraints
from scipy.special import erfc
from scipy.optimize import curve_fit
import numpy as np
import math
import matplotlib.pyplot as plt
def _bemg(x_in,la,lb, s, u,A):
    a=math.exp(la) #la,lb default 2
    b=math.exp(lb)
    x=(x_in-u)/s
    Eb=b*x+b**2/4
    Ea=a**2/4-a*x
    Eb[Eb>700]=700 #to avoid NaN
    Ea[Ea>700]=700
    ret=np.exp(Eb)*erfc(x+b/2) +np.exp(Ea)*erfc(-x+a/2)*3.4
    return A*ret

def _fit_bemg(y):
    N=len(y)
    x=np.linspace(-1,1,N)

    la=2
    lb=2
    ui=np.argmax(y)
    u=x[ui]
    A=y[ui]
    
    si=ui
    while(0<si and A/2<y[si]):si-=1
    ei=ui
    while(ei<N-2 and A/2 <y[ei] ): ei+=1
    s=(x[ei]-x[si])/2.86 

    p0=[la,lb,s,u,A]
    #print("BEMG intial p0",p0)
    p,_dummy = curve_fit(_bemg,x,y,p0=p0)
    #print("BEMG estimated p",p)
    yy=_bemg(x,*p)
    if False:
        plt.figure()
        plt.plot(y,label="in")
        plt.plot(yy,label="fited")
        plt.legend()    
        plt.pause(0.1)
    return yy
class ConstraintBEMG(pymcr.constraints.Constraint):
    """
    Parameters
    ----------
    copy : bool
        Make copy of input data, A; otherwise, overwrite (if mutable)
    """
    def __init__(self, peak_indexes , copy=False ):
        """ A must be non-negative"""
        super().__init__(copy)
        self.indexes=peak_indexes

    def transform(self, A):
        """ Apply nonnegative constraint"""
        C=A
        if self.copy:
            C=np.copy(A)
        
        for i in self.indexes:
            C[:,i]=_fit_bemg(C[:,i])
        return C


class FunBEMG():
    def __init__(self):
        pass
    def get_p0(self,x,y):
        N=len(y)
        la=2
        lb=2
        ui=np.argmax(y)
        u=x[ui]
        A=y[ui]
        
        si=ui
        while(0<si and A/2<y[si]):si-=1
        ei=ui
        while(ei<N-2 and A/2 <y[ei] ): ei+=1
        s=(x[ei]-x[si])/2.86     
        return la, lb, s, u, A
    def nof_p(self):
        return 5
    def fun(self,x,p):
        return _bemg(x,*p)

class FunLin():
    def __init__(self,deg):
        self.deg=deg

    def get_p0(self,x,y):
        p=np.polyfit(x,y,self.deg)
        return p.tolist()
    def nof_p(self):
        return self.deg+1
    def fun(self,x,p):
        return np.poy1d(x,p)

from pymcr.regressors import OLS, NNLS
class RegArray:
    coeff_ = None 
    def __init__(self,funs):
        self.funs=funs
        self.initial_reg=OLS()
    def fit(self,STT,DT):
        #初期値を
        self.initial_reg.fit(STT, DT)
        C_temp = self.initial_reg.coef_
        p0=[]
        N=C_temp.shape[0]
        x=np.linspace(-1,1,N)
        for n in range(C_temp.shape[1]):
            p0.extend(funs[i],get_p0(x,C_temp[:,n]))
        
        #p0に入った初期値を使ってcurve_fit
        ##########WIP########

        pass


if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    x=np.arange(100)*0.01
    plt.figure()
    plt.plot(_bemg(x,2,2,0.1,0.5,1))
    plt.show()