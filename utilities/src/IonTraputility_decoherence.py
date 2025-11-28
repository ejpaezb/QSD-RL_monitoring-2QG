from __future__ import division
import numpy as np
from scipy.constants import value, c, h, hbar, pi, k, epsilon_0

from fractions import Fraction
from scipy.interpolate import interp1d

am=value('atomic mass constant')
a0=value('Bohr radius')
ec=value('elementary charge')

"This notebook pins down the most noticeable decoherence processes, scattering rates and voltage fluctuations and its impact on the trap frequencies "


class decoherence(object):
    "All units are in SI compliance"
    
    def __init__(self,params):
        self.V = params['V']
        self.U = params['U']
        self.FS = params['FS']
        self.alpha = ec**2/(hbar*c*4*pi*epsilon_0)
        self.N = params['N']
        self.mass = params['mass']
        self.r = params['r0']
        self.rp = params['rp'] #radious of each patch potentials
        self.Z0 = params['Z0'] #characteristic impedance
        self.Rs = params['Rs'] #source impedance
        self.rms = params['rms']
        self.Cov = params['Cov']
        self.filteratt = params['filteratt']
        self.BW = params['BW']
        self.Q = params['Q']
        self.phi = params['phi']
        
        
        #self.ftrap[0]+1/self.ftrap[1])
        
    def raman_scatt(self,rabi,eta,gamma=2*pi*20e6,det=2*pi*33e12): #Rabi is given in Mrad/s, eta is Lamb-Dicke parameter, gamma is the linewidth of the intermediate state. Here taking linewidth of P1/2=P3/2
        rabi=1e6*rabi
        f = 1/290   #branching ratio to D-Levels from P levels in Yb 171 [1]
        eps = 1     #to avoid division by zero warnings
        FS = self.FS+eps #Fine spliting between P levels in Yb
        N = 2 #number of ions
        
        rate = 4/3*gamma*rabi*FS/np.abs(det*(det-FS))*(2*N) # the last factor includes a pair of Raman lasers and N ions  
        
        #for a time of pi/2, the probability is
        p = rate*pi/(2*rabi*eta) #it is equivalent to the error gate
        
        p_inf = 3*pi*gamma*f/FS*(2*N)/eta #probability of scattering to D level at infinite detuning
        
        ratio = 3*f/2*(2*det**2+(det-FS)**2)/FS**2 #ratio of probability of scattering D/S level
        
        
        return rate,p, p_inf, ratio
    
    def raman_scatt_power(self,rabi,waist,Power): #probability of error for a given fixed rabi, waist and laser power
        rabi=1e6*rabi

        eps_s = 2*pi*rabi*hbar*w**3*waist**2/(3*c**2*Power)
        
        return eps_s
    
    def rayleigh_scatt(self,rabi,eta,gamma=2*pi*20e6,det=2*pi*33e12): #rabi in rad/s, eta is lamb-dickey parameter
        rabi=1e6*rabi

        N = 2
        eps = 1
        FS = self.FS + eps
        rate = 4*N*rabi*gamma/FS*(3*det**2-2*det*FS+FS**2/3)/np.abs(det*(det-FS))
        
        p = rate*pi/(2*rabi*eta)
        
        error = p*5/12*eta**2
        
        ratio = 5/16*eta**2*(3*det**2-2*det*FS+FS**2/3)/FS**2
         
        p_inf = 4/eta*2*pi*gamma/FS
        
        error_inf = 5/2*pi*eta*gamma/FS
        
    #    print('rate = %0.4f' % rate, '1/s')
     #   print('prob = %0.7f' % p)
      #  print('prob_inf = %0.7f' % p_inf)
       # print('10^4 error = %0.7f' % (10**4*error))
        #print('10^4 error_inf = %0.7f' % (10**4*error_inf))
       # print('ratio Rayleigh/S level= %0.5f' % ratio)


        return rate,p,p_inf,error,error_inf,ratio
    
    
    def raman_scatt_error(self,power,time,freq,waist): #power in W, waist in m, freq in Hz, time in s.
        error = 8*pi**2*waist**2*(2*pi*c/329e-9)*171*am*(2*pi*freq)/(3*power*time)
        
        return error
    
    def cross_kerr(self,ws,wr,wz, nr=1): # ws: stretch mode, wr: radial, wz: COM of axial modes, nr: phonon number ocupation in the radial motional mode.
        chi = -ws*(2*hbar*wz/(self.alpha**2*171*am*c**2))**(1.0/3.0)*(wz/wr)*(1/2+ws**2/2/(4*wr**2-ws**2))*nr
 
        return chi
 

    def volt_fluct(self, freq, exp_rate = True): #coverage of patch potential regarding a characteristic sphere of radius rp, supression factor form filter, bandwidth, quality factor of trap circuit elements, average patch potential (less than 1) experimental rate if given
        
        #assuming white noise for voltage fluctuation and an efective patch potential \phi, there is a heating rate due to voltage fluctuation, the PSD is calculated rms/BW. [2]
        
        S_U = self.rms**2/self.BW*self.filteratt**2
        S_V = (self.Z0/(2*pi*self.self.Rs*self.Q))*(1*self.rms)**2/self.BW
        
        Gamma_U = self.N*ec**2*S_U/(4*self.mass*am*hbar*2*np.pi*freq[2]*self.Z0**2)*(1+(self.phi/self.U)**2)

        Gamma_V = 3*self.Cov*self.rp**2*S_V/(16*np.pi*self.r**4)*self.N*ec**2/(4*self.mass*am*hbar)(1/freq[0]+1/freq[1])/(2*pi)
        
        if exp_rate: rfvolt_fluct = None        
        else: rfvolt_fluct = np.sqrt(self.BW*exp_rate*(16*np.pi*self.r**4*4*self.mass*am*hbar*2*np.pi)/(1/freq[0]+1/freq[1]) /(3*self.N*ec**2*self.Cov*self.rp**2)) #if the experimental rate is given then calculate the corresponding voltage fluctuation
        
        return Gamma_U, Gamma_V, rfvolt_fluct #not including BBR heating rate due to resistive losses, the heating rate due to patch potentials and RF voltage fluctuations
        
    def whitenoise(self):
        V0 = self.V + np.random.normal(0,self.rms**2)
        return V0

# [1]: Errors in trapped-ion quantum gate due to spontaneous photon scattering. PRA 75, 042329 (2007)
# [2]: Ion trap measurements of electric filed noise near surfaces. Brownnutt et. al. Rev. Mod. Phys., Vol 87, No 4, p. 1452 eqn. 71, 2015