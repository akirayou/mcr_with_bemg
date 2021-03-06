
import pymcr.constraints
from scipy.special import erfcx
from scipy.optimize import curve_fit,minimize
import numpy as np
import math
import matplotlib.pyplot as plt
def _gauss(xi,u,s,A):
    x=(xi-u)/s
    return A*np.exp(-x**2)
def _fit_gauss(y):
    N=len(y)
    x=np.linspace(-1,1,N)
    ui=np.argmax(y)
    u=x[ui]
    A=y[ui]
    
    si=ui
    while(0<si and A/2<y[si]):si-=1
    ei=ui
    while(ei<N-2 and A/2 <y[ei] ): ei+=1
    s=(x[ei]-x[si])/2.355 

    p0=[u,s,A]
    p,_dummy = curve_fit(_gauss,x,y,p0=p0,bounds=( (-np.inf,1e-2,0 ), (np.inf,np.inf,np.inf)  ))
    yy=_gauss(x,*p)
    if False:
        plt.figure()
        plt.plot(y,label="in")
        plt.plot(yy,label="fited")
        plt.legend()    
        plt.pause(0.1)
    return yy
class ConstraintGauss(pymcr.constraints.Constraint):
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
            C[:,i]=_fit_gauss(C[:,i])
        return C

def _bemg(x_in,la,lb, u, s ,A):
    a=math.exp(la) #la,lb default 2
    b=math.exp(lb)
    x=(x_in-u)/s
    Ea=(a/2+x)
    Eb=(b/2-x)
    Ea[Ea<-25]=-25
    Eb[Eb<-25]=-25
    return A*a*b/(a+b)* np.exp(-x**2)* (erfcx(Eb)+erfcx(Ea))/1.1273434771994788

def _fit_bemg(y):
    N=len(y)
    x=np.linspace(-1,1,N)

    la=6
    lb=6
    ui=np.argmax(y)
    u=x[ui]
    A=y[ui]
    
    si=ui
    while(0<si and A/2<y[si]):si-=1
    ei=ui
    while(ei<N-2 and A/2 <y[ei] ): ei+=1
    s=(x[ei]-x[si])/2.355  

    p0=[la,lb,u,s,A]
    #print("BEMG intial p0",p0)
    p,_dummy = curve_fit(_bemg,x,y,p0=p0,bounds=( (-5,-5,-np.inf,1e-2,0 ), (10,10,np.inf,np.inf,np.inf)  ))
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

class FunGauss():
    def __init__(self):
        pass
    def get_p0(self,x,y):
        N=len(y)
        ui=np.argmax(y)
        u=x[ui]
        A=y[ui]
        
        si=ui
        while(0<si and A/2<y[si]):si-=1
        ei=ui
        while(ei<N-2 and A/2 <y[ei] ): ei+=1
        s=(x[ei]-x[si])/2.355    
        return  u, s, A
    def upper_bounds(self):
        return [np.inf,np.inf,np.inf]
    def lower_bounds(self):
        return [-np.inf,1e-2,0]
    def nof_p(self):
        return 3
    def fun(self,x,p):
        return _gauss(x,*p)

class FunBEMG():
    def __init__(self):
        pass
    def get_p0(self,x,y):
        N=len(y)
        la=6
        lb=6
        ui=np.argmax(y)
        u=x[ui]
        A=y[ui]
        
        si=ui
        while(0<si and A/2<y[si]):si-=1
        ei=ui
        while(ei<N-2 and A/2 <y[ei] ): ei+=1
        s=(x[ei]-x[si])/2.355      
        return la, lb, u, s, A
    def upper_bounds(self):
        return [10,10,np.inf,np.inf,np.inf]
    def lower_bounds(self):
        return [-5,-5,-np.inf,1e-2,0]
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
    def upper_bounds(self):
        return [np.inf]*self.nof_p()
    def lower_bounds(self):
        return [-np.inf]*self.nof_p()
        
    def fun(self,x,p):
        return np.poly1d(p)(x)
from pymcr.regressors import OLS, NNLS
class RegArray:
    coef_ = None 
    def __init__(self,funs):
        self.funs=funs
        self.initial_reg=NNLS()
    def fit(self,STT,DT):
        #初期値を
        self.initial_reg.fit(STT, DT)
        C_temp = self.initial_reg.coef_
        p0=[]
        N=C_temp.shape[0]
        x=np.linspace(-1,1,N)
        bounds=[[],[]]
        for n in range(C_temp.shape[1]):
            p0.extend(self.funs[n].get_p0(x,C_temp[:,n]))
            #plt.figure("R{}".format(n))
            #plt.clf()
            #plt.plot(C_temp[:,n])
            #plt.plot(self.funs[n].fun(x,p0[-self.funs[n].nof_p():  ]))
            #plt.pause(0.1)
            bounds[0].extend(self.funs[n].lower_bounds())
            bounds[1].extend(self.funs[n].upper_bounds())
        def fun(x,*p):
            pos=0
            C=np.zeros((DT.shape[1],STT.shape[1]))
            for i,fun in enumerate(self.funs):
                pp=p[pos:pos+fun.nof_p()]
                C[:,i]=fun.fun(x,pp)
                pos+=fun.nof_p()
            return C
        def funD(x,*p):
            return np.dot(fun(x,*p),STT.T).reshape(-1)
        
        #p=p0        
        p=curve_fit(funD,x,DT.T.reshape(-1),p0=p0,bounds=bounds)[0]
        #print("p" ,p)
        C=fun(x,*p)
        #plt.figure("C")
        #plt.clf()
        #plt.plot(C)
        #plt.pause(10)
        self.coef_=C


if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    x=np.arange(100)*0.01
    plt.figure()
    plt.plot(_bemg(x,2,2,0.1,0.5,1))
    plt.show()