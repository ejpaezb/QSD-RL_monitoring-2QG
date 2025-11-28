from qutip import *
from scipy.constants import value, c, h, hbar, pi, k, epsilon_0

import numpy as np
import src.ACstarkutility as ac
import src.IonTraputility as ion
import src.IonTraputility_decoherence as decoherence
import src.IonTraputility_crystal as crystal

am=value('atomic mass constant')
a0=value('Bohr radius')
ec=value('elementary charge')
 
class Sim(object): 
    """
    Ion trap simulator
    Parameters:
        f:           The function to be optimized (Objective Function)
        bounds:      The bounds on optimization parameters
        nPopulation: The population size
        F :          The difference amplification factor (optional)
        CR :         The cross-over probability (optional)
        seed :       The seed for the random number generator (optional)
        minimize:    If True the function "f" will be minimized else maximized (optional)
    """

    
    @staticmethod
    def __init__(params):
        Sim.t_gate = params['t_gate'] # Duration of the gate
        Sim.n_dim = params['n_dim'] # Dimension of the motional state
        Sim.delta = params['delta'] # Detuning from carrier transition
        Sim.HFS = params['HFS']
        Sim.omega = []
        Sim.stark_shift = []
        Sim.pulse_steps = params['pulse_steps']
        Sim.gamma_heating = params['gamma_heating'] # Heating rate (500 quanta/second)
        Sim.gamma_spin_dephasing = params['gamma_spin_dephasing'] # Spin dephasing rate
        Sim.gamma_spin_flip = params['gamma_spin_flip'] # Spin flip rate, both Raman and Rayleigh
        Sim.tau = params['tau'] # Motional coherence time
        Sim.V0 = params['V']
        Sim.U0 = params['U']
        Sim.rms = params['rms']
        Sim.ion1 = params['ion1']
        Sim.ion2 = params['ion2']
        
        Sim.trap = None
        Sim.Yb = None
        Sim.decoherence = None
        #Sim.opts = Options(nsteps = 1000)
        Sim.opts = Options()
        ############################################################################################################
        #### Atomic units for StarkShift calculation
        Sim.au = a0*ec
        Sim.Yb = ac.ACStark(params)
        Sim.Yb.import_leveldata('./src/yb_energy_levels.csv'); #The energy levels from the NIST database and A. Roy et al
        Sim.Yb.import_matrixelements('./src/yb_matrix_elements.csv');
        Sim.Yb.populate_transitiontable();

        lamt = Sim.Yb.transtable[3,0]
        Sim.Yb.Lpow = params['Lpow']       #Laser power [W]
        Sim.Yb.Lwaist = params['Lwaist']   #waist [m]
        Sim.Ramandet = params['Ramandet']  #detuning Hz
        f = c/lamt + Sim.Ramandet
        Sim.Yb.lam = c/f  #wavelength [m]       
        Sim.Yb.parameters_pol([0,0,1]) # Linear polirization
        
        Ef = ac.electric_field(Sim.Yb.Lwaist, Sim.Yb.Lpow)
        Sim.term0 = Sim.Yb.termlabels[0]
        Sim.shift0,Sim.terms = Sim.Yb.substate_starkshift(term=Sim.term0,F=0,mF=0, field=0) #field=0 means the Ef given by local variables power and waist
        
        Sim.matrix_ele = Sim.Yb.matrixelements[0] #Selecting the ground state as the qubit 0

        ############################################################################################################
        ### Initializing the trap parameters and decoherence
        Sim.trap = ion.IonTrap(params)
        Sim.dec  = decoherence.decoherence(params)
        Sim.trap.all_parameter()

        ############################################################################################################
        ### Initializaing the crystal
        Sim.crys = crystal.IonCrystal(params)
 
        Sim.eigen = Sim.crys.eigenparams() #'vector, values,alpha crit.'
        Sim.equi_pos = Sim.crys.equi_position(Sim.trap.ftrap) #'equlibrium positions, minimum separation, characteristic length, anisotropy'

        Sim.Anm = Sim.crys.axial_matrix()
        Sim.Bnm = Sim.crys.radial_matrix(Sim.Anm)

        Sim.freq = Sim.crys.vib_modes(Sim.trap.ftrap[2]) #returns the frequencies of the vibrationa modes both radial and axial axis, respectively

        Sim.nu_radial = Sim.freq[0] # COM frequency for radial axis
        Sim.nu_axial = Sim.freq[1] # COM frequency for axial axis
        Sim.line_width = 1e-6 * Sim.nu_radial[0]
        Sim.lamb_dicke = Sim.crys.lamb_dicke(Sim.trap.ftrap[2])
        
        
    @staticmethod
    def setPulse(pulse):
        ##global omega
        Sim.omega = np.copy(pulse)
        
        ##global pulse_steps
        Sim.pulse_steps = len(Sim.omega)
        
        ##reset the stark shift
        Sim.stark_shift = np.zeros_like(pulse)
        
    @staticmethod
    def setDelta(delta): #The input should be in Mrad/s in the same units as the stark shift given by setSatarkShift
        Sim.delta = delta
        #print(Sim.lamb_dicke)

    @staticmethod
    def setStarkShift(pulse):
        #The input must be in units of MHz. Here it is converted to Hz. This static method delivers in Hz the Stark shift of the two hyperfine levels involved. The gllbal variable stark_shift is conveniently converted to Mrad/s, this units is convenient as the solver is in units of Mrad/s and time in \mu s.
        
        pulse = np.abs(np.real(pulse))*1e6
        electricf = Sim.Yb.rabi2E(Sim.matrix_ele*Sim.au,pulse,Sim.Ramandet)

        Sim.Yb.det = 0
        ss_0 = np.zeros(len(pulse))
        #ts=np.zeros(5)
        for i in np.arange(len(pulse)):
            ss_0[i],ts = Sim.Yb.substate_starkshift(term=Sim.term0,F=0,mF=0, field=electricf[i])

        Sim.Yb.det = Sim.HFS
        ss_1 = np.zeros(len(pulse))
        #ts=np.zeros(5)
        for i in np.arange(len(pulse)):
            ss_1[i],ts = Sim.Yb.substate_starkshift(term=Sim.term0,F=1,mF=0, field=electricf[i])

        
        ss = np.zeros(len(pulse))
        ss = ss_1-ss_0
        #            ss[i],ts=Sim.Yb.substate_starkshift(term=Sim.term0,F=1,mF=0, field=electricf[i])

        ##global stark_shift
        Sim.stark_shift = (2*np.pi) * (np.copy(ss)) / 1e6 ## Hz to Mega rad/s

        return ss ## This is in Hertz 

    @staticmethod
    def setNoisySidebands(): #sets the noise to the motional frequencies based on random white noise (\nu), also sets the 'new' detuning due to the stark shift (\mu)
        noisy_mode = []
        Sim.nu_radial_noisy = []
        Sim.nu_axial_noisy = []

        for i in np.arange(1,500+1):
            Sim.trap.U = Sim.U0 + np.random.normal(0,1*Sim.rms**2)
            Sim.trap.all_parameter()
            noisy_mode.append(Sim.trap.ftrap[2])

        Sim.nu_radial_noisy, Sim.nu_axial_noisy = Sim.crys.vib_modes(np.array(noisy_mode))
        
        Sim.nu1 = 2*pi*1e-6*Sim.nu_radial_noisy[0]
        Sim.nu2 = 2*pi*1e-6*Sim.nu_radial_noisy[1]
        Sim.nu3 = 2*pi*1e-6*Sim.nu_radial_noisy[2]
        Sim.nu4 = 2*pi*1e-6*Sim.nu_radial_noisy[3]
        Sim.nu5 = 2*pi*1e-6*Sim.nu_radial_noisy[Sim.trap.N-1]
        Sim.mu = Sim.delta + np.array(Sim.stark_shift)
        
        Sim.eta_1_1 = Sim.lamb_dicke[0][Sim.ion1] #\eta_{ion}_{mode}:
        Sim.eta_2_1 = Sim.lamb_dicke[0][Sim.ion2] #lamb_dicke[mode][ion]
        Sim.eta_1_2 = Sim.lamb_dicke[1][Sim.ion1]
        Sim.eta_2_2 = Sim.lamb_dicke[1][Sim.ion2]
        Sim.eta_1_3 = Sim.lamb_dicke[2][Sim.ion1]
        Sim.eta_2_3 = Sim.lamb_dicke[2][Sim.ion2]
        Sim.eta_1_4 = Sim.lamb_dicke[3][Sim.ion1]
        Sim.eta_2_4 = Sim.lamb_dicke[3][Sim.ion2]
        Sim.eta_1_5 = Sim.lamb_dicke[4][Sim.ion1]
        Sim.eta_2_5 = Sim.lamb_dicke[4][Sim.ion2]

        
    #### Full Hamiltonian, expanding to the first order in \eta. The terms are written in the same order that they appear in the overleaf document 'Light-atom interaction Sorensen-Molmer': https://www.overleaf.com/2163481312xkrbdkmybfks        
        
    @staticmethod
    def H0_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)


    @staticmethod
    def H1_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(1.0j*Sim.nu1[idx_nu]*t)


    @staticmethod
    def H2_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(1.0j*Sim.nu2[idx_nu]*t)


    @staticmethod
    def H3_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(1.0j*Sim.nu3[idx_nu]*t)

    @staticmethod
    def H4_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(1.0j*Sim.nu4[idx_nu]*t)

    @staticmethod
    def H5_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(1.0j*Sim.nu5[idx_nu]*t)

    @staticmethod
    def H6_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(-1.0j*Sim.nu1[idx_nu]*t)


    @staticmethod
    def H7_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(-1.0j*Sim.nu2[idx_nu]*t)

    
    @staticmethod
    def H8_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(-1.0j*Sim.nu3[idx_nu]*t)

    @staticmethod
    def H9_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(-1.0j*Sim.nu4[idx_nu]*t)

    @staticmethod
    def H10_t(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        idx_nu = int(t*(499/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.cos((Sim.mu[idx])*t)*np.exp(-1.0j*Sim.nu5[idx_nu]*t)


    @staticmethod
    def gamma_motion_dephasing_coeff(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return np.sqrt(2 * decoherence.rayleigh_scatt(Sim.omega[idx], 1)[0])

    @staticmethod
    def gamma_spin_dephasing_coeff(t, args):
        idx = int(t*(Sim.pulse_steps/Sim.t_gate))
        if (idx == Sim.pulse_steps):
            idx -= 1
        return Sim.omega[idx]*np.sqrt(Sim.gamma_spin_dephasing)

#    @staticmethod
 #   def gamma_spin_flip_coeff(t, args):
  #      idx = int(t*(Sim.pulse_steps/Sim.t_gate))
   #     if (idx == Sim.pulse_steps):
    #        idx -= 1
     #   return np.sqrt(decoherence.raman_scatt(Sim.omega[idx], 1)[0])

    @staticmethod
    def integrate(solver="me", t_g=100, progress_bar=True, collapse = False):

        ##global t_gate
        Sim.t_gate = t_g

        # quantum trajectories to allow mixed states
        sm_1 = tensor(destroy(2), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{-} for the 1st ion
        sm_2 = tensor(qeye(2), destroy(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{-} for the 2nd ion

        sx_1 = tensor(sigmax(), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{x} for the 1st ion
        sx_2 = tensor(qeye(2), sigmax(), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{x} for the 2nd ion

        sy_1 = tensor(sigmay(), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{y} for the 1st ion
        sy_2 = tensor(qeye(2), sigmay(), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{y} for the 2nd ion

        sz_1 = tensor(sigmaz(), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{z} for the 1st ion
        sz_2 = tensor(qeye(2), sigmaz(), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Sigma_{z} for the 2nd ion

        a_1 = tensor(qeye(2), qeye(2), destroy(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Annihilation operator for the 1st mode
        a_2 = tensor(qeye(2), qeye(2), qeye(Sim.n_dim), destroy(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Annihilation operator for the 2nd mode
        a_3 = tensor(qeye(2), qeye(2), qeye(Sim.n_dim),qeye(Sim.n_dim), destroy(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)) # Annihilation operator for the 3rd mode
        a_4 = tensor(qeye(2), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), destroy(Sim.n_dim),qeye(Sim.n_dim)) # Annihilation operator for the 1st mode
        a_5 = tensor(qeye(2), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),destroy(Sim.n_dim)) # Annihilation operator for the 1st mode
        
        
        # Useful operators
        p_gg_gg = tensor(ket2dm(basis(2,0)), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim))
        p_ee_ee = tensor(ket2dm(basis(2,1)), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim))

        # |00>
        rho_ideal_n = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,0), basis(2,0))-1.0j*tensor(basis(2,1), basis(2,1)))), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim))
        rho_ideal_p = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,0), basis(2,0))+1.0j*tensor(basis(2,1), basis(2,1)))), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim))

        # |11>
        #rho_ideal_n = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,1), basis(2,1))-1.0j*tensor(basis(2,0), basis(2,0)))), qeye(Sim.n_dim), qeye(Sim.n_dim))
        #rho_ideal_p = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,1), basis(2,1))+1.0j*tensor(basis(2,0), basis(2,0)))), qeye(Sim.n_dim), qeye(Sim.n_dim))

        # |01>
        #rho_ideal_n = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,0), basis(2,1))-1.0j*tensor(basis(2,1), basis(2,0)))), qeye(Sim.n_dim), qeye(Sim.n_dim))
        #rho_ideal_p = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,0), basis(2,1))+1.0j*tensor(basis(2,1), basis(2,0)))), qeye(Sim.n_dim), qeye(Sim.n_dim))

        # |10>
        #rho_ideal_n = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,1), basis(2,0))-1.0j*tensor(basis(2,0), basis(2,1)))), qeye(Sim.n_dim), qeye(Sim.n_dim))
        #rho_ideal_p = tensor(ket2dm((1.0/np.sqrt(2.0))*(tensor(basis(2,1), basis(2,0))+1.0j*tensor(basis(2,0), basis(2,1)))), qeye(Sim.n_dim), qeye(Sim.n_dim))

        I = tensor(qeye(2), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim))

        N_1 = a_1.dag() * a_1
        N_2 = a_2.dag() * a_2
        N_3 = a_3.dag() * a_3
        N_4 = a_4.dag() * a_4
        N_5 = a_5.dag() * a_5

        # Phase-space q & p operators
        q_1 = a_1 + a_1.dag() # 1st motional mode
        p_1 = a_1 - a_1.dag()

        q_2 = a_2 + a_2.dag() # 2nd motional mode
        p_2 = a_2 - a_2.dag()

        q_3 = a_3 + a_3.dag() # 2nd motional mode
        p_3 = a_3 - a_3.dag()

        q_4 = a_4 + a_4.dag() # 1st motional mode
        p_4 = a_4 - a_4.dag()

        q_5 = a_5 + a_5.dag() # 2nd motional mode
        p_5 = a_5 - a_5.dag()

        # List of operators that we want to see their expectation value
        e_ops = []
        e_ops.append(p_gg_gg)
        e_ops.append(p_ee_ee)
        e_ops.append(rho_ideal_p)
        e_ops.append(rho_ideal_n)
        e_ops.append(N_1)
        e_ops.append(N_2)
        e_ops.append(N_3)
        e_ops.append(N_4)
        e_ops.append(N_5)
        e_ops.append(q_1)
        e_ops.append(p_1)
        e_ops.append(q_2)
        e_ops.append(p_2)
        e_ops.append(q_3)
        e_ops.append(p_3)
        e_ops.append(q_4)
        e_ops.append(p_4)
        e_ops.append(q_5)
        e_ops.append(p_5)


        # Carrier transition
        H0  = 1*(sx_1 + sx_2)

        # ...
        H1  = - 1*a_1.dag() * (Sim.eta_1_1 * sy_1 + Sim.eta_2_1 * sy_2) #eta_{ion}_{mode}
        
        # ...
        H2  = - 1*a_2.dag() * (Sim.eta_1_2 * sy_1 + Sim.eta_2_2 * sy_2)
        
        # ...
        H3  = - 1*a_3.dag() * (Sim.eta_1_3 * sy_1 + Sim.eta_2_3 * sy_2)

        # ...
        H4  = - 1*a_4.dag() * (Sim.eta_1_4 * sy_1 + Sim.eta_2_4 * sy_2)

        # ...
        H5  = - 1*a_5.dag() * (Sim.eta_1_5 * sy_1 + Sim.eta_2_5 * sy_2)
        
        # ...
        H6  = - 1*a_1 * (Sim.eta_1_1 * sy_1 + Sim.eta_2_1 * sy_2) #eta_{ion}_{mode}

        # ...
        H7  = - 1*a_2 * (Sim.eta_1_2 * sy_1 + Sim.eta_2_2 * sy_2)

        # ...
        H8  = - 1*a_3 * (Sim.eta_1_3 * sy_1 + Sim.eta_2_3 * sy_2)

        # ...
        H9  = - 1*a_4 * (Sim.eta_1_4 * sy_1 + Sim.eta_2_4 * sy_2)

        # ...
        H10 = - 1*a_5 * (Sim.eta_1_5 * sy_1 + Sim.eta_2_5 * sy_2)
        # ...

        H = []
        if solver == "mc":
            H.append(tensor(qeye(2), qeye(2), qeye(Sim.n_dim), qeye(Sim.n_dim),qeye(Sim.n_dim)))
        H.append([H0, Sim.H0_t])
        H.append([H1, Sim.H1_t])
        H.append([H2, Sim.H2_t])
        H.append([H3, Sim.H3_t])
        H.append([H4, Sim.H4_t])
        H.append([H5, Sim.H5_t])
        H.append([H6, Sim.H6_t])
        H.append([H7, Sim.H7_t])
        H.append([H8, Sim.H8_t])
        H.append([H9, Sim.H9_t])
        H.append([H10, Sim.H10_t])
        
        
        # Build collapse operators
     
        c_ops = []
        if collapse:
            # Motional Heating
            c_ops.append(np.sqrt(Sim.gamma_heating)*a_1) # Heating of the 1st motional mode
            c_ops.append(np.sqrt(Sim.gamma_heating)*a_1.dag()) # ...

            c_ops.append(np.sqrt(Sim.gamma_heating)*a_2) # Heating of the 2nd motional mode
            c_ops.append(np.sqrt(Sim.gamma_heating)*a_2.dag()) # ...

            # Motional Dephasing
            # c_ops.append(a_1.dag()*a_1, gamma_motion_dephasing_coeff) # Dephasing of the 1st motional mode         

            # c_ops.append(a_2.dag()*a_2, gamma_motion_dephasing_coeff) # Dephasing of the 2nd motional mode

            # Spin Dephasing
            #c_ops.append([sz_1, self.gamma_spin_dephasing_coeff]) # Dephasing of the 1st ion
            #c_ops.append([sz_2, self.gamma_spin_dephasing_coeff]) # Dephasing of the 2nd ion

            # Spin Flip
            c_ops.append(sm_1 * np.sqrt(Sim.gamma_spin_flip)) # Spin Flip of the 1st ion
            c_ops.append(sm_1.dag() * np.sqrt(Sim.gamma_spin_flip)) # ...

            c_ops.append(sm_2 * np.sqrt(Sim.gamma_spin_flip)) # Spin Flip of the 2st ion
            c_ops.append(sm_2.dag() * np.sqrt(Sim.gamma_spin_flip)) # ...

        
        print('t_gate: ', Sim.t_gate)
        tlist = np.linspace(0.0, Sim.t_gate, 1000)
        
        output = []
        if solver == "me": # me: Direct integration
            # |00>
            psi0 = (tensor(basis(2,0), basis(2,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0))) # Initial state of the system
            print('psi_0 = |00>')
            # |11>
            #psi0 = (tensor(basis(2,1), basis(2,1), basis(Sim.n_dim,0), basis(Sim.n_dim,0))) # Initial state of the system
            # |01>
            #psi0 = (tensor(basis(2,0), basis(2,1), basis(Sim.n_dim,0), basis(Sim.n_dim,0))) # Initial state of the system
            # |10>
            #psi0 = (tensor(basis(2,1), basis(2,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0))) # Initial state of the system
            #print('psi_0 = |10>')
            ## ...
            #psi0 = ket2dm(tensor(basis(2,0), basis(2,0), basis(n_dim,0), basis(n_dim,0))) # Initial state of the system
            output = mesolve(H, psi0, tlist, c_ops, e_ops, options=Sim.opts, progress_bar=progress_bar) # Numerical integration
        elif solver == "mc": # mc: Quantum trajectories method
            psi0 = tensor(basis(2,0), basis(2,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0),basis(Sim.n_dim,0), basis(Sim.n_dim,0), basis(Sim.n_dim,0)) # Initial state of the system
            output = mcsolve(H, psi0, tlist, c_ops, e_ops, options=Sim.opts, progress_bar=progress_bar, ntraj=100) # Quantum trajectories method

        return output

    