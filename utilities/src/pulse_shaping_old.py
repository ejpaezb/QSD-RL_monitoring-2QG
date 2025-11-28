import numpy as np
import scipy as sp
import scipy.linalg

from tqdm import tqdm

import pickle
 
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

class IonTrap(object):
    """IonTrap class."""
    def __init__(self, N, ion_i, ion_j, omegas_path, g_ki_path, steps, nbar_c, time):
        self.N = N # Number of ions in the chain

        # Ion pair affected by the gate
        self.ion_i = ion_i
        self.ion_j = ion_j

        # Gate time
        self.gate_time = time

        # Longitudinal motional mode frequencies
        self.omegas = np.load(omegas_path)

        # Coupling factor between ith ion and kth transverse motional mode
        self.g_ki = np.load(g_ki_path)

        # Detunings list
        self.mu_list = np.linspace(0.99 * self.omegas[N - 1], 1.01 * self.omegas[0], 99)
#         self.mu_list = np.linspace(0.9 * self.omegas[N - 1], 1.1 * self.omegas[0], 200)

        self.steps = steps # Default number of pulse segments
        self.nbar_c = nbar_c # Mean phonon number in each mode
        self.rabi = np.array([0.0]) # Rabi rate

    def b(self, i, k):
        """
        b_{i,k} is the normal eigen-mode transformation matrix for ith ion and kth transverse motional mode.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).

        k : int
            Transverse motional mode's number (k >= 0).
            
        Returns
        -------
        out : b_{i,k}
        """    
        _b = np.array([np.array([1.0, 1.0]), np.array([1.0, -1.0])]) # 2 ions
        # _b = np.array([(1/np.sqrt(3.0)) * np.array([1.0, 1.0, 1.0]), (1/np.sqrt(2.0)) * np.array([-1.0, 0.0, 1.0]), (1/np.sqrt(6.0)) * np.array([1.0, -2.0, 1.0])]) # 3 ions
        return _b[k][i]

    def eta(self, k):
        """
        The Lamb-Dicke parameter for the kth motional mode.

        Parameters
        ----------
        k : int
            Transverse motional mode's number (k >= 0).
            
        Returns
        -------
        out : \eta_{k}
        """
        
        _k = np.array([0.06502, 0.06531]) # 2 ions
        # _k = np.array([0.06501954, 0.06530712, 0.06572272]) # 3 ions
        return _k[k]
        
    def g(self, i, k):
        """
        g_{i,k} is the coupling factor between ith ion and kth transverse motional mode.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        k : int
            Transverse motional mode's number (k >= 0).
            
        Returns
        -------
        out : b_{i,k}
        """
        # return b(i, k) * eta(k)
        return self.g_ki[k][i]
    
    def alpha(self, i, k, tau, mu):
        """
        The \alpha parameter as defined on pages 88 and 89 of the Manning's PhD thesis.
        It represents the coupling of ions' electronic structure to the motional modes.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        k : int
            Transverse motional mode's number (k >= 0).
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : \alpha(i, k, tau, mu)
        """
        # 1-segment pulse
        if (len(self.rabi) == 1): # analytical method
            w, z = 0, tau
            _term0 = (np.exp(1.0j*z*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(z*mu)+mu*np.sin(z*mu))-np.exp(1.0j*w*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(w*mu)+mu*np.sin(w*mu)))/(mu**2 - self.omegas[k]**2)
            return (-1.0 * 1.0j * self.rabi[0] *  self.g(i, k) * _term0)

        # multi-segment pulse
        else:
            C_mat = np.zeros([ self.steps, 1], dtype=complex)
            for p in range(1,  self.steps + 1):
                C_mat[p - 1, 0] = self.C_generator(i, k, p, tau, mu)
            return np.matmul(C_mat.T, self.rabi)[0][0]

    def alpha_manual(self, i, k, tau, mu):
        """
        The \alpha parameter as defined on pages 88 and 89 of the Manning's PhD thesis.
        It represents the coupling of ions' electronic structure to the motional modes.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        k : int
            Transverse motional mode's number (k >= 0).
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : \alpha(i, k, tau, mu)
        """
        cnt = int(tau * (self.steps / self.gate_time))
        if cnt == self.steps:
            cnt -= 1
            
        dt = self.gate_time / self.steps

        _term0 = 0
        if cnt == 0:
            _term0 += self.C_generator_manual(i=i, k=k, mu=mu, upper_bound=tau, lower_bound=0.0) * self.rabi[0]
        else:
            for idx in range(cnt):
                _term0 += self.C_generator_manual(i=i, k=k, mu=mu, upper_bound=dt*(idx+1), lower_bound=dt*idx) * self.rabi[idx]

            if tau - dt * (cnt) > 0:
                _term0 += self.C_generator_manual(i=i, k=k, mu=mu, upper_bound=tau, lower_bound=dt*cnt) * self.rabi[cnt]
            
        return _term0

    def chi(self, i, j, tau, mu):
        """
        The \chi parameter as defined on pages 88 and 89 of the Manning's PhD thesis.
        It represents the geometric phase gained by the ions after doing Molmer-Sorensen scheme.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        j : int
            Ion's number (j >= 0).
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : \chi(j, k, tau, mu) ???
        """
        # 1-segment pulse
        _chi = 0.0
        if (len(self.rabi) == 1):
            for k in range(self.N):
                _chi += 2 * self.rabi[0]**2 *  self.g(i, k) * self.g(j, k) * self.integ(k, tau, mu)
            # _term0 = 2 * self.rabi[0]**2 * g(i, 0) * g(j, 0) * integ(0, tau, mu)
            # _term1 = 2 * self.rabi[0]**2 * g(i, 1) * g(j, 1) * integ(1, tau, mu)
            
            return _chi #_term0 + _term1
        
        # multi-segment pulse
        else:
            D_mat = np.zeros([ self.steps,  self.steps], dtype=complex)
            for p in range(1,  self.steps + 1):
                for pp in range(1,  self.steps + 1):
                    D_mat[p - 1, pp - 1] = self.D_generator(i, j, p, pp, tau, mu)
            
            return (1/2)*np.matmul(self.rabi.T, np.matmul(D_mat, self.rabi))[0][0]
        
    def integ(self, k, tau, mu):
        """
        A helper function used in "def chi(...)". It reperesents the solution to the
        integral in the definition of \chi.

        Parameters
        ----------        
        k : int
            Transverse motional mode's number (k >= 0).
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : integ(...)
        """
        _term0 = -(self.omegas[k]*(np.sin(2*mu*tau)+2*mu*tau))/(4*(mu**3-mu*self.omegas[k]**2))
        _term1 = (self.omegas[k]*(mu*np.sin(mu*tau)*np.cos(self.omegas[k]*tau))-self.omegas[k]*np.cos(mu*tau)*np.sin(self.omegas[k]*tau))/(mu**2-self.omegas[k]**2)**2
        
        return _term0 + _term1


    def beta(self, k):
        """
        A helper function used in "def fidelity(...)".
        
        Parameters
        ----------        
        k : int
            Transverse motional mode's number (k >= 0).
                    
        Returns
        -------
        out : beta(k)
        """
        _tmp = (1.0 / 2.0) * np.log(1.0 + 1.0 / self.nbar_c)
        ##return np.cosh(_tmp) / np.sinh(_tmp)
        return 1.0

    def gamma_mn(self, i, tau, mu):
        """
        A helper function used in "def fidelity(...)".
        
        Parameters
        ----------        
        i : int
            Ion's number (i >= 0).
            
        tau : float 
            Interaction time.
                    
        mu : float
            Detuning from motional sideband.
        
        Returns
        -------
        out : beta(k)
        """
        _tmp = 0
        for k in range(self.N):
            _tmp += np.abs(self.alpha(i, k, tau, mu))**2 * (self.beta(k) * 2.0)

        return np.exp(-_tmp)

    def gamma_p(self, tau, mu):
        """
        A helper function used in "def fidelity(...)".
        
        Parameters
        ----------        
        tau : float 
            Interaction time.
                    
        mu : float
            Detuning from motional sideband.
        
        Returns
        -------
        out : gamma_p(tau, mu)
        """
        _tmp = 0
        for k in range(self.N):
            _tmp += np.abs(self.alpha(self.ion_i, k, tau, mu) + self.alpha(self.ion_j, k, tau, mu))**2 * (self.beta(k) * 2.0)

        return np.exp(-_tmp)

    def gamma_m(self, tau, mu):
        """
        A helper function used in "def fidelity(...)".
        
        Parameters
        ----------        
        tau : float 
            Interaction time.
                    
        mu : float
            Detuning from motional sideband.
        
        Returns
        -------
        out : gamma_m(tau, mu)
        """
        _tmp = 0
        for k in range(self.N):
            _tmp += np.abs(self.alpha(self.ion_i, k, tau, mu) - self.alpha(self.ion_j, k, tau, mu))**2 * (self.beta(k) * 2.0)

        return np.exp(-_tmp)

    def fidelity(self, tau, mu):
        """
        Calculates the gate's fidelity. For more information please refer 
        to eq. 4.8 on page 103 of Manning's PhD thesis.

        Parameters
        ----------
        tau : float 
            Interaction time.
            
        mu : float
            Detuning from motional sideband.
            
        Returns
        -------
        out : fidelity(tau, mu)
        """
        return (1.0 / 8.0) * (2.0 + 1.0j * (np.exp(-2j*np.abs(self.chi(self.ion_i, self.ion_j, tau, mu))) - np.exp(2j*np.abs(self.chi(self.ion_i, self.ion_j, tau, mu)))) * (self.gamma_mn(self.ion_i, tau, mu) + self.gamma_mn(self.ion_j, tau, mu)) + self.gamma_p(tau, mu) + self.gamma_m(tau, mu))

    def setPulse(self, pulse):
        """
        Sets the global pulse shape.

        Parameters
        ----------
        pulse : array_like 
            1-D array representing the global pulse shape.
                    
        Returns
        -------
        None
        """        
        self.rabi = np.copy(pulse)
        
    def find_rabi(self, i, j, tau, mu):
        """
        Returns the rabi rate that gives the highest gate fidelity for a one-segment constant pulse.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        j : int
            Ion's number (j >= 0).
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : find_rabi(i, j, tau, mu)
        """
        _chi = 0.0
        for k in range(self.N):
            _chi += 2 *  self.g(i, k) * self.g(j, k) * self.integ(k, tau, mu) 
        
        _bestRabi = np.sqrt(np.abs(np.pi/(4 * _chi)))
        
        self.setPulse(np.array([_bestRabi]))
        return _bestRabi

    def C_generator_manual(self, i, k, mu, upper_bound, lower_bound):
        """
        The C matrix as defined on page 100 the Manning's PhD thesis.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        k : int
            Transverse motional mode's number (k >= 0).
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : C_generator(i, k, p, tau, mu)
        """        
        w, z = lower_bound, upper_bound
        
        _term0 = (np.exp(1.0j*z*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(z*mu)+mu*np.sin(z*mu))-np.exp(1.0j*w*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(w*mu)+mu*np.sin(w*mu)))/(mu**2 - self.omegas[k]**2)

        return -1.0j *  self.g(i, k) * _term0

    def C_generator(self, i, k, p, tau, mu):
        """
        The C matrix as defined on page 100 the Manning's PhD thesis.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        k : int
            Transverse motional mode's number (k >= 0).
            
        p : int
            The pulse segmentation variable as introduced in the thesis.
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : C_generator(i, k, p, tau, mu)
        """
        dt = (tau /  self.steps)
        w, z = (p-1)*dt, p*dt
        
        _term0 = (np.exp(1.0j*z*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(z*mu)+mu*np.sin(z*mu))-np.exp(1.0j*w*self.omegas[k])*(1.0j*self.omegas[k]*np.cos(w*mu)+mu*np.sin(w*mu)))/(mu**2 - self.omegas[k]**2)

        return -1.0j *  self.g(i, k) * _term0

    def B_generator(self, p, pp, tau, mu):
        """
        The B matrix as defined on page 104 the Manning's PhD thesis.

        Parameters
        ----------
        p : int
            The pulse segmentation variable as introduced in the thesis.
            
        pp : int
            The pulse segmentation variable as introduced in the thesis.
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : B_generator(p, pp, tau, mu)
        """
        _B = 0.0
        for k in range(self.N):
            __term0 = np.conj(self.C_generator(self.ion_i, k, p, tau, mu)) * self.C_generator(self.ion_i, k, pp, tau, mu)
            __term1 = np.conj(self.C_generator(self.ion_j, k, p, tau, mu)) * self.C_generator(self.ion_j, k, pp, tau, mu)
            _B += self.beta(k) * (__term0 + __term1)

        return (1.0/4.0) * _B 
                    
    def D_generator(self, i, j, p, pp, tau, mu):
        """
        The D matrix as defined on page 101 the Manning's PhD thesis.

        Parameters
        ----------
        i : int
            Ion's number (i >= 0).
            
        j : int
            Ion's number (j >= 0).
            
        p : int
            The pulse segmentation variable as introduced in the thesis.
            
        pp : int
            The pulse segmentation variable as introduced in the thesis.
            
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        Returns
        -------
        out : D_generator(i, j, p, pp, tau, mu)
        """
        if (pp > p):
            p, pp = pp, p
        
        _D = 0.0
        for k in range(self.N):
            _D += 2 *  self.g(i, k) * self.g(j, k) * self.integ_D_gen(k, tau, mu, p, pp)
            
        return _D
            
    def integ_D_gen(self, k, tau, mu, p, pp):
        """
        A helper function used in "def D_generator(...)". It represents the solution to the
        integral in the definition of D_{p,p^{\prime}}.

        Parameters
        ----------
        k : int
            Transverse motional mode's number (k >= 0).
                
        tau : float 
            Interaction time.
            
        mu: float
            Detuning from motional sideband.
            
        p : int
            The pulse segmentation variable as introduced in the thesis.
        
        pp : int
            The pulse segmentation variable as introduced in the thesis.
            
        Returns
        -------
        out : integ_D_gen(k, tau, mu, p, pp)
        """
        dt = (tau /  self.steps)

        if (p == pp):
            w, z, m = (p-1)*dt, p*dt, (pp-1)*dt
            
            _term0 = (1/(4*mu*(mu**2-self.omegas[k]**2)**2))*(2*(w-z)*mu*(mu-self.omegas[k])*self.omegas[k]*(mu+self.omegas[k])-4*mu**2*self.omegas[k]*np.cos(z*mu)*np.cos((-m+z)*self.omegas[k])*np.sin(m*mu)-4*mu**2*self.omegas[k]*np.cos(m*mu)*np.cos((-m+w)*self.omegas[k])*np.sin(w*mu)+mu**2*self.omegas[k]*np.sin(2*w*mu)-self.omegas[k]**3*np.sin(2*w*mu)+4*mu**2*self.omegas[k]*np.cos(m*mu)*np.cos((-m+z)*self.omegas[k])*np.sin(z*mu)-mu**2*self.omegas[k]*np.sin(2*z*mu)+self.omegas[k]**3*np.sin(2*z*mu)+4*mu**3*np.sin(m*mu)*np.sin(w*mu)*np.sin((-m+w)*self.omegas[k])+4*mu*self.omegas[k]*np.cos(w*mu)*(mu*np.cos((-m+w)*self.omegas[k])*np.sin(m*mu)+self.omegas[k]*np.cos(m*mu)*np.sin((-m+w)*self.omegas[k]))-4*mu*(self.omegas[k]**2*np.cos(m*mu)*np.cos(z*mu)+mu**2*np.sin(m*mu)*np.sin(z*mu))*np.sin((-m+z)*self.omegas[k]))
            _term0 *= 2
        else:
            w, z, m, n = (p-1)*dt, p*dt, (pp-1)*dt, pp*dt
            
            _term0 = (1/(2*(mu**2-self.omegas[k]**2)**2))*(mu*((mu+self.omegas[k])*np.cos(w*(mu-self.omegas[k])+m*self.omegas[k])-(mu+self.omegas[k])*np.cos(z*mu+m*self.omegas[k]-z*self.omegas[k])-(mu-self.omegas[k])*(np.cos(m*self.omegas[k]-w*(mu+self.omegas[k]))-np.cos(m*self.omegas[k]-z*(mu+self.omegas[k]))))*np.sin(m*mu)+mu*(-(mu+self.omegas[k])*np.cos(w*(mu-self.omegas[k])+n*self.omegas[k])+(mu+self.omegas[k])*np.cos(z*mu+n*self.omegas[k]-z*self.omegas[k])+(mu-self.omegas[k])*(np.cos(n*self.omegas[k]-w*(mu+self.omegas[k]))-np.cos(n*self.omegas[k]-z*(mu+self.omegas[k]))))*np.sin(n*mu)+2*self.omegas[k]*np.cos(m*mu)*(-mu*np.cos((m-w)*self.omegas[k])*np.sin(w*mu)+mu*np.cos((m-z)*self.omegas[k])*np.sin(z*mu)-self.omegas[k]*np.cos(w*mu)*np.sin((m-w)*self.omegas[k])+self.omegas[k]*np.cos(z*mu)*np.sin((m-z)*self.omegas[k]))+self.omegas[k]*np.cos(n*mu)*((mu+self.omegas[k])*np.sin(w*(mu-self.omegas[k])+n*self.omegas[k])-(mu+self.omegas[k])*np.sin(z*(mu-self.omegas[k])+n*self.omegas[k])+(mu-self.omegas[k])*(np.sin(w*mu-n*self.omegas[k]+w*self.omegas[k])-np.sin(z*mu-n*self.omegas[k]+z*self.omegas[k]))))
        
        return _term0

    def selectOptimizedPulse(self, stp, tau, mu):
        self.steps = stp
        
        # Generating the D matrix
        D_mat = np.zeros([ self.steps,  self.steps], dtype=complex)

        for p in range(1,  self.steps + 1):
            for pp in range(1,  self.steps + 1):
                D_mat[p - 1, pp - 1] = self.D_generator(self.ion_i, self.ion_j, p, pp, tau, mu)

        # Generating the B matrix
        B_mat = np.zeros([ self.steps,  self.steps], dtype=complex)
        for p in range(1,  self.steps + 1):
            for pp in range(1,  self.steps + 1):
                B_mat[p - 1, pp - 1] = self.B_generator(p, pp, tau, mu)

                
        # Solving the generalized eigenvalue problem
        aa = (B_mat + B_mat.T)
        bb = (D_mat + D_mat.T)/2

        eigVal, eigVec = sp.linalg.eig(aa, bb)
        
        f_list = []
        # Selecting non-zero pulses
        # Calculating each pulse shape fidelity
        for i in range(self.steps):

            scale = (np.pi/4.0) / ((1.0 / 2.0) * np.matmul(eigVec[:, i].T, np.matmul(D_mat, eigVec[:, i])))            
            scale = np.sqrt(np.abs(scale))            

            # Scaling the solution to meet the geometric phase requirement
            eigVec[:, i] = eigVec[:, i] * scale
            
            pulse = eigVec[:, i].reshape(self.steps, 1)
            self.setPulse(pulse)

            max_rabi = np.max(np.abs(pulse))
            fidelity = self.fidelity(tau, mu)

            if (fidelity > 0.999):                
                f_list.append([i, max_rabi, fidelity, pulse])
                # print(max_rabi, self.fidelity(tau, mu))

        # sort the results based on the scales
        f_list.sort(key=lambda x: float(x[1]), reverse=False)
        # print(f_list)

        _bestF = f_list[0][2]
        _bestIdx = f_list[0][0]
        _optPulse = f_list[0][3]
        _maxRabi = np.max(np.abs(f_list[0][3]))
        _scaleFactor  = 1.0
            
        return _bestF, _optPulse, _maxRabi, _scaleFactor

    def generatePulse(self, time, n_segments):
        # Multi-segment pulse
        f_list = []
        f_rabi = []
        pulse_list = []

        scales = []
        chis = []

        for muTmp in tqdm(self.mu_list):
            _bestF, _optPulse, _maxRabi, _scale = self.selectOptimizedPulse(n_segments, time, muTmp)
            f_list.append(_bestF)
            f_rabi.append(_maxRabi)
            pulse_list.append(_optPulse)

            self.setPulse(_optPulse)
            scales.append(_scale)
            chis.append(self.chi(self.ion_i,self.ion_j,time,muTmp))

        # 1-segment pulse
        f_list_single = []
        f_rabi_single = []
        for muTmp in self.mu_list:    
            f_rabi_single.append(self.find_rabi(self.ion_i, self.ion_j, time, muTmp))
            f_list_single.append((self.fidelity(time, muTmp)))

        return pulse_list