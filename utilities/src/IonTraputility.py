from __future__ import division
import numpy as np
from scipy.constants import value, c, h, hbar, pi, k, epsilon_0
from fractions import Fraction
from scipy.interpolate import interp1d

 
am=value('atomic mass constant')
a0=value('Bohr radius')
ec=value('elementary charge')

"Calculates the stability region of a Linear Paul trap useful for trap design and evaluation of trapping stability under common experimental fluctations. This design is meant for Quantum gate implementation"

class IonTrap(object):
    "All units are in SI compliance"
    
    def __init__(self,params):
        self.f = params['f'] #Define the RF frequency [Hz]
        self.mass = params['mass'] #Define mass [a.u]
        self.V = params['V'] #Define the RF voltage [V peak]
        self.U = params['U'] #Define DC voltage [V]
        self.n = params['n'] #order of truncation for recursive equations
        self.m = params['m'] #elements of fit vectors
        self.r = params['r0'] # radious (distance) from the trap centre to the surface of the rods
        self.N = params['N'] #number of ions
        self.z0 = params['z0'] # half the distance between DC end-caps
        self.lambd = params['lambd']
        self.angle = params['angle']
        self.HFS = params['HFS']
        
        self.a = np.zeros(3)
        self.q = np.zeros(3)
        self.beta = np.zeros(3)
        self.ftrap = np.zeros(3)
        self.lamb_dicke = np.zeros(3)
        
    def all_parameter(self):
        self.a_q()
        for i in range(3):
            D0,ker,a,q,beta0 = self.Dn(self.a[i],self.q[i],0)
            self.beta[i] = self.beta_interp(ker,beta0,i)
            
        self.ftrap = self.f*self.beta
        
        return None
    
    def a_q(self): #mass,frequency,RF voltage [peak], DC voltage, radius, 
        n = 1 #particle charge in # of 'elementary charge'
        k = 1
        ax = -4*k*n*ec*self.U/(self.mass*am*(2*pi*self.f)**2*self.z0**2)
        ay = ax
        qx = 2*n*ec*self.V/(self.mass*am*(2*pi*self.f)**2*self.r**2)
        qy = -qx
        az = -2*ax
        qz = 0
        self.a = np.array([ax,ay,az])
        self.q = np.array([qx,qy,qz])
        
        return None #returns the a and q values. These are crutial to undertanding the stability of the design


    def Dn(self, a0, q0, printing, beta0=None): #, max 'q' to scan, max 'a' to scan, beta value: if this field is blank it means that the beta value is to be calculated, otherwise, it calculates the set of a-q values for given beta.
        m = self.m
        n = self.n
        gap = 1e-5
        beta = None
        
        if beta0 != None: #it takes a0 and q0 as the max number to construct the grid used to calculate the corresponding values (a and q) matching the given beta
            
            a = np.linspace(-a0,a0,2*m)
            q = np.linspace(gap,q0,m)
            A,Q = np.meshgrid(a,q)
            D = np.empty(shape=(2*n+1,m,2*m))

            for i in np.arange(-n,n+1):
                
                aux = (A-(2*i+beta0)**2)/Q
                D[i+n,:,:] = aux
                
            ker,k01,k02=self.kn(D)
            
            if (printing == 1):
                max_0 = np.max(np.abs(ker))
                self.mypltcolor(q,a,np.transpose(ker/max_0),200)
                plt.xlabel('q',fontsize=16)
                plt.ylabel('a',fontsize=16)
                plt.text(0.1,-0.6,'The frontier black-white represents the set ') 
                plt.text(0.1,-0.7,'of a-q values solution of Mathieu equation')
                plt.text(0.1,-0.9,'Beta= '+str(beta))


        elif q0 != 0:
            #if only a-q are given (scalars taken from self.a_q) then the following sentences generate a vector of beta values as a search space for the correct value
            a = a0
            q = q0
            est_beta=np.sqrt(a+q**2/2)
            #beta = np.linspace(0,1,2*m)
            beta = np.linspace(0.5*est_beta,min(1.5*est_beta,1),2*m)
            D = np.empty(shape=(2*n+1,2*m))
            #print('Radial calculation done!')
            for i in np.arange(-n,n+1):
                aux = (a0-(2*i+beta)**2)/q0
                D[i+n,:] = aux
            ker,k01,k02 = self.kn_beta(D)
                
        else: #When q=0 the result is straigthforward. So to avoid division by zero in the calculation of ker, all of them are set to None. The function self.beta_interp is in charge to calculate beta when q=0
            a = a0
            q = q0
            D = None
            ker = None
            #print('Axial calculation done')
                 
                
        return D, ker, a, q, beta
   

    def kn(self, D): #This function calculates the set of a-q parameters satisfying a given beta parameter, receives the vector of 'D' values determined in self.Dn, same order of truncation 'n'
        n = self.n
        k01 = 1/(D[-1+n,:,:]-1/(D[-2+n,:,:]-1/D[-3+n,:,:]))
        k02 = 1/(D[1+n,:,:]-1/(D[2+n,:,:]-1/D[3+n,:,:]))
            
        ker = D[n,:,:]-k01-k02
            
        return ker, k01, k02 #To understand this expression and derivation check the notes YYY in dropbox folder. The variable 'ker' is passed on to self.myinterpolation

    
    def kn_beta(self, D): #This function is only used to calculate beta for a given a-q values, receives the vector of 'D' values determined in self.Dn, same order of truncation 'n'
        n = self.n
        k01 = 1/(D[-1+n,:]-1/(D[-2+n,:]-1/D[-3+n,:]))
        k02 = 1/(D[1+n,:]-1/(D[2+n,:]-1/D[3+n,:]))
        ker = D[n,:]-k01-k02
            
        return ker, k01, k02 #To understand this expression and derivation check the notes YYY in dropbox folder. The variable 'ker' is passed on to self.myinterpolation

    
    def Dn_interm(self, beta, a, q): #Verification. order of the approximation, beta value, q vector, a vector. This verification function differs from self.Dn in the sense that it is used to verify that a particular set of 'a' and 'q' values generates the same 'beta'.
        n = self.n
        m = self.m
        gap = 1e-10
        D = np.empty(shape=(2*n+1,m))
        
        for i in np.arange(-n,n+1):
            aux = (a-(2*i+beta)**2)/q
            D[i+n,:] = aux
               
        return D # this result is passed on to self.beta
    
    
    def beta_interm(self, a, q, D): #Verification. Calculates the beta parameter for a given 'a', 'q' vectors and 'D' array. This verification function works together with self.Dnfinal
        n = self.n
        k01 = 1/(D[-1+n]-1/(D[-2+n]-1/D[-3+n]))
        k02 = 1/(D[1+n]-1/(D[2+n]-1/D[3+n]))
        beta0 = np.sqrt(a-q*(k01+k02))
        
        return beta0
    
    
    def parameter_interp(self, ker, a, q): #Fit the plot corresponding to the solution of the stability region
        
        m = self.m
        
        idx = np.where(np.abs(ker)<0.5e-4)
        agen = None
        qgen = None
        
            
        idx_q = list(idx[0])
        idx_a = list(idx[1])
        anew = a[idx_a]
        qnew = q[idx_q]
        
        f = interp1d(qnew,anew) # fitting function aproximating the a-q solution for a particular beta
        plt.plot(qnew,anew,'r')
            
            #Once the fitting function is found a new smooth sequence of a-q is generated
        qgen = np.linspace(qnew[0],qnew[len(qnew)-1],m)
        agen = f(qgen)
               
        return agen,qgen

    def beta_interp(self, ker, beta_v,i): #determines the beta value from the zero-crossing of the ker function. beta_v is a vector
        
        if self.q[i] != 0:
            
            m = self.m
            idx = np.where(np.abs(ker)<1e-3)
            #plt.plot(ker)
            beta0 = np.mean(beta_v[idx])
        else:
            beta0 = np.sqrt(self.a[i])
            
        return beta0
    
        
    def mypltcolor(self, x, y, z, f): #Only to handle large amount of data to be plotted. It simply resizes the vectors 'x', 'y' and 'z'
        #f=200 #number of points
        lx = int(len(x)/f)
        ly = int(len(y)/f)
        x1 = x[0 : len(x) : lx] #sampling the list using steps lx. Applies to 'y' and 'z'
        y1 = y[0 : len(y) : ly]
        z1 = z[0 : len(y) : ly, 0 : len(x) : lx]
        #print(x1.size,y1.size,z1.shape)
        
        plt.pcolor(x1,y1,np.log(np.abs(z1))); plt.set_cmap('hot'); plt.colorbar();    
        plt.figure()
        plt.pcolor(x1,y1,np.angle(z1)); plt.set_cmap('hot'); plt.colorbar();
        
        return None
        
    def depth(self): #Calculate the potential depth in the radial confinement.
        aux = ec*self.V**2/(4*self.mass*am*self.r**2*(2*pi*self.f)**2) #This result is already in eV
        Temp = aux*ec/k
        
        return aux, Temp

        
    def doppler_temp(self,Gamma,Omega): #linewidth and Rabi strength
        s=2*(Omega/Gamma)**2
        T= hbar*Gamma*np.sqrt(1+s)*(1+2/5)/(4*k)
        Delta = np.sqrt(1+s)/2
        Ed = hbar*Gamma*(1+2/5)/8*(np.sqrt(1+s)*2) #Energy per (radial) axis. for optimun detuning Delta. Includes secular and micromotion as if they carry the same energy.
        
        Total_E = np.array([Ed,Ed,1/2*Ed]) #three axis taking into account that the axial contribution is due to secular motion
        
        nbar = Ed/(hbar*self.ftrap)
        
        return T,nbar, Total_E,Delta #Doppler temperature limit and theoretical detuning
        
    
    def doppler_shift(self,nbar,Energy_th): #linewidth, detuning, Rabi strength
        mod_ind=[]
        k=2*pi/self.lambd
        mod_ind = 2*self.q*(self.lamb_dicke)**2*(2*nbar+1)/4 #modulation index
        Efield_emm =  self.mass*am*(2*pi*self.f)**2*mod_ind/(k*ec)    #Electric field due to excess micromotion (EMM)
        vel= mod_ind*self.f*self.lambd
        ui=ec*Efield_emm/(self.mass*am*self.ftrap**2)
        
        Energy_emm = 1/(self.mass*am)*(ec*100*Efield_emm/(2*2*pi*self.ftrap))**2
        doppler_shift = -(Energy_emm+Energy_th)/(self.mass*am*c**2)
        
        return doppler_shift, Efield_emm, vel, ui
        
    def SB_cooling(self,Gamma,Omega,det):
        alpha=2/5
        
        if (Omega < Gamma/10):
            nbar = (alpha+1)*(Gamma/2/self.ftrap)**2
            print('one')
        elif det < 2*max(self.ftrap):
            nbar = (Omega/self.ftrap)**2/8
            print('two')
        else: 
            nbar = (alpha+1)/4*det/self.ftrap
        
        return nbar
    
