
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
    u=np.argmax(y)
    A=y[u]
    s=u
    while(0<s and A/2<y[s]):s-=1
    e=u
    while(e<N-2 and A/2 <y[e] ): e+=1
    s=(e-s)/2.86 * (2/N)
    p0=[2,2,s,x[u],A]
    #print("BEMG intial p0",p0)
    p,_dummy = curve_fit(_bemg,x,y,p0=p0)
    print("BEMG estimated p",p)
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
if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    x=np.arange(100)*0.01
    plt.figure()
    plt.plot(_bemg(x,2,2,0.1,0.5,1))
    plt.show()