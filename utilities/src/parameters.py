import numpy as np


class parameters(object):
    
    @staticmethod
    def main():

############## ############## ############## ############## ############## ############## ##############       
        # Utilized in IonTraputility.py
        params = {'init':None}
        params.update({'t_gate' : 101}) # Duration of the gate
        params.update({'U': 3}) #DC voltage applied to end-caps
        params.update({'I': 0.5}) #Nuclear spin of Ytterbium 171
        params.update({'mass' : 171}) #define the mass in a.u
        params.update({'lambd' : 355e-9}) #Laser wavelength
        params.update({'f' : 20e6})        #Define the RF frequency [Hz]
        params.update({'n' : 3}) #Order of the truncation for recursive equations: "Dn" and "Cn"
        params.update({'m' : 1500})   #Length of the vectors     
        params.update({'r0' : 305e-6}) #trap radius in [m]. Spacing between rods
        params.update({'z0' : 1e-3}) #End caps spacing
        params.update({'V' : 397.48}) #Trap voltage
        params.update({'Laser_conf' : 'counter'}) # Laser configuration of the Raman Lasers: none: single laser, co: co-propag.
        params.update({'angle' : np.pi/4}) # incident angle of Raman lasers (regarding the trap axial axis)
        
############## ############## ############## ############## ############## ############## ##############              
        # Utilized in IonTraputiliy_simulation
        params.update({'t_gate': 100})# Duration of the gate
        params.update({'n_dim': 10})  # Dimension of the motional state
        params.update({'delta': 0.0})     # Detuning from carrier transition
        params.update({'HFS': 12.64e9})
        params.update({'pulse_steps': 10})
        params.update({'gamma_heating': 5e-4})     # Heating rate (500 quanta/second)
        params.update({'gamma_spin_dephasing': 1e-6})     # Spin dephasing rate
        params.update({'gamma_spin_flip': 2*30e-6})     # Spin flip rate, both Raman and Rayleigh
        params.update({'tau': 1e5})     # Motional coherence time
        params.update({'N': 2})     # Number of ions
        params.update({'ion1': 0})     # position of the first ion involved in the two qubit gate
        params.update({'ion2': 1})     # position of the second ion involved in the two qubit gate

############## ############## ############## ############## ############## ############## ##############       
        # Utilized in AcStarkShift
        params.update({'Lpow' : 1000e-3})   #Laser power [W]
        params.update({'Lwaist' : 30e-6 })  # waist [m]
        params.update({'Ramandet' : 33e12}) #detuning Hz
        params.update({'det' : 0})
        
############## ############## ############## ############## ############## ############## ##############       
        # Utilized in decoherence        
        params.update({'rms': 0*1.5e-3})     # Motional coherence time
        params.update({'FS' : 2*np.pi*100e12}) #Fine splitting between 2P1/2 and 2P3/2
        params.update({'rp' : 2e-6}) #radious of each patch potentials
        params.update({'Z0' : 50}) #characteristic impedance. The letter 'Z' is capital and shouldn't be confused with z0
        params.update({'Rs' : 50}) #source impedance
        params.update({'Cov' : 0.5})
        params.update({'filteratt' : 1e-3})
        params.update({'BW' : 20e6})
        params.update({'Q' : 500})
        params.update({'phi' : 1}) 
        
############## ############## ############## ############## ############## ############## ##############       
        # Utilized in the project file @CEDAR
        params.update({'det_0' : 56}) #fixed detuning of the delta_list frequency range
        params.update({'pulses' : 7}) # Number of pulses. The optimum is 2N+1
        params.update({'presc' : 0.01}) #prescision used in the DE algorith for the Rabi pulse generation 
        params.update({'maxR' : 1.0}) # Max Rabi interval of searching

             

############## ############## ############## ############## ############## ############## ##############       
        # Utilized in IonTraputility_crystal
        
        
        return params