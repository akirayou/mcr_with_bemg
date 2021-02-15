"""
Simulate chromagram peak shape wiht plate theory
"""
import sim_chrom.isotherm as isotherm
import numpy as np
class SimChrom:
    def _dv_input(self):
        if(self.inject_width>0):self.inject_width-=1

        if(self.dv<1):
            ret=self._input
            if(self.inject_width<=0):self._input=0
        else:
            ret=self._input / self.dv
            if(self.inject_width<=0):self._input-=ret
        return ret

      
    def __init__(self,n_plates : int, iso : isotherm.Isotherm,dead_volume : float =0,pass_rate : float =0,inject=1,inject_width=1 ):
        """
        Parameters
        -----------
        n_plates : int 
            number of plates
        iso : isotherm.Isotherm
            isotherm model
        dead_voume : float
            dead volume tailing  simulation (under 1 means no dead volume. )
        """
        self.dv=dead_volume
        self._input=inject
        self.inject_width=inject_width
        self._output=0
        self.iso=iso
        self.moving=np.zeros(n_plates)
        self.static=np.zeros(n_plates)
        self.moving[0]=1
        self.pass_rate=pass_rate
        for _i in range(n_plates-1):
            self.step()
    
    def step(self):
        """
        get one sample intensity.
        """
        self.moving=np.roll(self.moving,1)
        o=self.moving[0]
        self.moving[0]=self._dv_input()
        passing=self.moving*self.pass_rate
        S=self.moving*(1-self.pass_rate)+self.static
        self.moving,self.static=self.iso.redist(S)
        self.moving+=passing
        return o

if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    chrom=SimChrom(1000,isotherm.IsoRP(1,1,1.2),dead_volume=100)
    ret=[ chrom.step() for i in range(1000)]
    print(ret)
    plt.figure()
    #plt.yscale("log")
    plt.plot(ret)
    plt.show()


        