"""
isotherm functions
See: https://www.jstage.jst.go.jp/article/oleoscience/2/5/2_275/_pdf

Important variables
C:濃度
W:吸着量 
S=W+C
"""

from abc import (ABC as _ABC, abstractmethod as _abstractmethod)
import numpy as np
class Isotherm(_ABC):
    """ 
    Isotherm settings
    """
    def __init__(self,eps=1e-10,alpha=0.1):
        self.eps=eps
        self.alpha=alpha
        
    @_abstractmethod
    def iso(self, W):
        """ Calc W from C
        
        Parameters
        ----------
        W : real
            Weight of absorbed amount of static phase(normalized)
        
        Returns
        -------
        C : real
            Concentration of moving phase (normalized)
        """

    def _CfromS(self, S):
        """ Calc C from S
        Override this if you have analytic solution.
        This function calcs simple (but slow) numerical method.
        To tune this, self.eps for precision, self.alpha for calculation stabirity and speed

        Parameters
        ----------
        S : Sum of W[moving phase] and C[static phase]

        Returns
        -------
        C : real
            Concentration of moving phase (normalized)
        """
        oldC=S
        C=S
        while(True):
            W=self.iso(oldC)
            C=S-W
            C[C<0]=0 # workaround for calculation error
            res=np.max(np.abs(oldC-C))
            if(np.isnan(res)): raise(ValueError("NaN at S="+str(S)))
            #print(res)
            if( res <self.eps):break
            
            C=oldC*(1-self.alpha)+C*self.alpha
            oldC=C
        return C


    
    def redist(self,S):
        """ Redistribute sample to W and C
        
        Parameters
        ----------
        S : Sum of W[moving phase] and C[static phase]

        Returns
        -------
        C : real
            Concentration of moving phase (normalized)
        W : real
            Concentration of moving phase (normalized)
        """
        C = self._CfromS(S)
        W = S - C
        return C, W


class IsoRP(Isotherm):
    """Radke-Prausnitz Isotherm"""
    def __init__(self,a,b,m,eps=1e-7,alpha=0.1):
        super().__init__(eps,alpha)
        self.a=a
        self.b=b
        self.m=m
    
    def iso(self,C):
        return  (C**(self.m+1)*self.a*self.b)/(C**self.m*self.b+C*self.a+1e-200)

class IsoLang(Isotherm):
    """Langmuir"""
    def __init__(self,a,Ws,eps=1e-7,alpha=0.1):
        super().__init__(eps,alpha)
        self.a = a
        self.Ws = Ws
    
    def iso(self,C):
        return self.a*self.Ws*C*(1+self.a*C)

class IsoFre(Isotherm):
    """Freundlich"""
    def __init__(self,K,n,eps=1e-7,alpha=0.1):
        super().__init__(eps,alpha)
        self.K = K
        self.n = n
    
    def iso(self,C):
        return self.K*C**(1/self.n)


class IsoGLang(Isotherm):
    """Generaized Langmuir"""
    def __init__(self,b,Ws,n,m,eps=1e-7,alpha=0.1):
        super().__init__(eps,alpha)
        self.b = b
        self.Ws = Ws
        self.n=n  # 0<n<=1
        self.m=m  # 0<m<=1
    def iso(self,C):
        return self.Ws*    ((self.b*C)**self.n /(1+ (self.b*C)**self.n ))**(self.m/self.n) 
if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    import numpy as np
    c=np.linspace(0,1,10000)
    
    isos=[IsoRP(1,1,0.9),IsoLang(1,1),IsoFre(1,1.2)]
    labels=["RP","Lang","Fre"]
    plt.figure()
    for iso,l in zip(isos,labels):
        plt.plot(c,iso.iso(c),label=l)
    plt.legend()
    plt.show()

    plt.figure()
    for iso,l in zip(isos,labels):
        plt.plot(c,iso.redist(c)[1]  ,label=l)
    plt.legend()
    plt.show()



