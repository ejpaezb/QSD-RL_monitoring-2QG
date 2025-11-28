from __future__ import division
import numpy as np
from scipy.constants import value, c, h, hbar, pi, k, epsilon_0
from fractions import Fraction
from scipy.interpolate import interp1d

from openpyxl import load_workbook
from sympy.functions.special.tensor_functions import KroneckerDelta

am=value('atomic mass constant')
a0=value('Bohr radius')
ec=value('elementary charge')

"Calculates the stability region of a Linear Paul trap useful for trap design and evaluation of trapping stability under common experimental fluctations. This design is meant for Quantum gates experiment"

class IonCrystal:
    "All units are in SI compliance"
    def __init__(self, params):
        self.N = params['N'] #number of trapped ions
        self.f = params['f'] #Define the RF frequency [Hz]
        self.mass = params['mass'] #Define mass [a.u]
        self.wx = None
        self.wy = None
        self.wz = None
        self.Laser_conf = params['Laser_conf']
        self.HFS = params['HFS']
        self.lambd = params['lambd']
        self.angle = params['angle']


        #self.alpha=(self.wx/self.wz)**2 #anisotropy of the trap
 
    def parameters(self): # Delivers few constants necessary to calculate the ion equilibrium position in self.equi_position
        A = 0.436*self.N**(0.596)
        B = 0.0375*self.N**(-0.178)
        alpha = np.sqrt(4*A/(3*B))
        beta = np.sqrt(27*B/(4*A**3))
        
        return alpha,beta
     
    def equi_position(self,f): #see T.P Meyrath, D.F.V James phys rev lett A 240 (1998): "theoretical and numerical studies of the positions of cold trapped ions"
        length = (ec**2/(4*pi*epsilon_0*self.mass*am*(2*pi*f[2])**2))**(1/3) #Lenght scale of the trap. Only depends on the confinement on 'Z'
        alpha = (f[2]/f[0])**2
        u_n = np.empty(shape=(self.N))
        
        if self.N==2:
            u_n[0]=-(1/2)**(2/3)
            u_n[1]=(1/2)**(2/3)
        elif self.N==3:
            u_n[0]=-(5/4)**(1/3)
            u_n[1]=0
            u_n[2]=(5/4)**(1/3)
        else:
            k1,k2=self.parameters()
            for i in range(self.N):
                u_n[i]=k1*np.sin(1/3*np.arcsin(k2*(i+1-(self.N+1)/2)))
        
        u_min=2.29*self.N**(-0.596)
        self.u_n = u_n
        self.alpha = alpha
        
        return u_n, u_min,length, alpha #return the ion position relative to the trap centre in micrometers. Also, the minimum distance, the characteristic scale and the anisotropy
    
    def eigenparams(self):
        vector = load_workbook("./src/bnm.xlsx")
        N = self.N
        d = ('A','B','C','D','E','F','G','H','I','J','K')
        var = str(N)
        c = vector[var]
        evector = np.empty(shape=(N,N))
        evalue = np.empty(shape=(N))
        
        for i in range(N):
            evalue[i]=c[d[0]+str(i+1)].value
            for j in range(N):
                var2=d[j+1]+str(i+1)
                evector[i,j]=c[var2].value
        
        alpha_c = 2/(evalue[N-1]-1)
        self.mu = evalue
        return evector,evalue,alpha_c #eigenvectors, eigenvalues, alpha critical
    
    def axial_matrix(self):
        N = self.N
        Anm = np.empty(shape=(N,N))
         
        for i in range(N):
            for j in range(N):
                if i==j:
                    aux=np.delete(self.u_n,i)
                    Anm[i,j]=1+2*np.sum(1/np.abs(self.u_n[i]-aux)**3)
                else:
                    Anm[i,j]=-2/np.abs(self.u_n[i]-self.u_n[j])**3
                
        return Anm    
     
    def radial_matrix(self,Anm):
        
        gamma = 1/self.alpha+1/2-self.mu/2
        N = self.N
        Bnm = np.empty(shape = (N,N))
        
        for i in range(N):
            for j in range(N):
                Bnm[i,j] = (1/self.alpha+1/2)*KroneckerDelta(i,j)-1/2*Anm[i,j]
        self.gamma = gamma
        return Bnm,gamma #weight matrix for radial modes, eigenvalue (needed for calculating the vibrational frequencies)
    
    def vib_modes(self,f_z): #f_z can be the axial trap frequency as an array: noisy side_bands due to U fluctuation
        N = self.N
        freq_r = []
        freq_a = []
        lamb_dicke_r = []
        for i in range(self.N):
            freq_r.append(f_z*np.sqrt(self.gamma[i]))
            freq_a.append(f_z*np.sqrt(self.mu[i]))
        return freq_r, freq_a
    
    def lamb_dicke(self,freq_z):
        lambd2 = c/(self.HFS + c/self.lambd)
        lamb_dicke_r = []
        evector,evalue,alpha_c = self.eigenparams()
        
        if self.Laser_conf == 'counter':            #Calculates \Delta k for the counter or co-propgating config.
            delta_k = 2*pi*(1/self.lambd + 1/lambd2)*np.cos(self.angle)
        elif self.Laser_conf == 'co':
            delta_k = 2*pi*(1/self.lambd - 1/lambd2)*np.cos(self.angle)
        else: delta_k = 2*pi/self.lambd*np.cos(self.angle)
        
        freq_r,freq_a = self.vib_modes(freq_z)
        for i in range(self.N):
            lamb_dicke_r.append(delta_k*np.sqrt(hbar/(2*self.mass*am*2*pi*np.array(freq_r[i])))*evector[i])
            #lamb_dicke_r.append(delta_k*np.sqrt(hbar/(2*self.mass*am*2*pi*np.array(freq_r[i]))))

               
            
        return lamb_dicke_r #returns the lamb_dicke paramter for the vibrational modes of the radial axis
