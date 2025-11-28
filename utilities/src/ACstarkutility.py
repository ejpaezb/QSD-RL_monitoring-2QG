from __future__ import division
import numpy as np
import csv
from fractions import Fraction
from src.wigner import Wigner3j,Wigner6j

from scipy.constants import value, epsilon_0, hbar, c, pi, h,k
a0 = value('Bohr radius')
ec = value('elementary charge')

"""
Substantially Modified by Eduardo Paez 31-01-2017:
Following the model and the results given in Safronova, M. S., Johnson, W. R., & Derevianko, A. (1999). Relativistic many-body   calculations of energy levels, hyperfine constants, electric-dipole matrix elements, and static polarizabilities for alkali-metal atoms. Physical Review A, 60(6), 4476.

Here we determine the scalar, vector and tensor polarisability. Previous version only included the reduced (static) scalar polarizability and the scalar component.
"""


class ACStark:
    """Class for general calculation of AC Stark shifts in ion traps.
    
    All quantities are in SI units.
    """
    def __init__(self,params):
        self.Lpow = None #Laserpower [W], default
        self.Lwaist = None #Waist [m], default
        self.lam = None #Wavelength [m]
        self.nuc_I = params['I'] #Nuclear spin, set to Ytterbium 171 = 1/2
        self.transtable = None #Auxiliar Table for preparing atomic data 
        self.pol = np.array([0,0,1]) #Laser polarization (default set to pi)
        self.K=np.array([])         #For tensor analysis contribution
        self.det= params['det'] #in case of need, a negative detuning from the given frequency is available
        self.A= None #elipticity, given by function parameters_pol()
         
    def parameters_pol(self,polvec):
        """Set the beam polarization in spherical basis.
            Take a 3 component vector input (sig+,sig-,pi) and 
            renormalizes it"""
        e_r=-1/np.sqrt(2)*np.array([1,1j,0]) #spherical basis vector circularly right polarised light
        e_l=1/np.sqrt(2)*np.array([1,-1j,0]) # spherical basis vector circularly left polarised light
        e_z=np.array([0,0,1]) # spherical basis pi polarised light. Equivalent to z-axis
        
        polvec = np.array(polvec) #re-cast as np array
        polvec = polvec/np.linalg.norm(polvec) #normalize
        self.pol = polvec
        lightvec = polvec[0]*e_r+polvec[1]*e_l+polvec[2]*e_z
        self.A = np.real(1j*np.cross(lightvec,np.conj(lightvec))[2]) #Elipticity (-1<=A=<1)
                 
        return polvec
    
    def import_leveldata(self,filename):
        """Reads in csv file with pairs consisting of a term label and energy in
            inverse cm"""
        termlabels = []
        termqnumbers = []
        energies = []
        with open(filename,'r') as file:
            csv_read = csv.reader(file,delimiter=',')
            for line in csv_read:
                
                termlabels.append(line[0])
                energies.append(float(line[1]))
                termqnumbers.append(parse_term(line[0]))
        file.close()
        
        #sort by energy 
        energies,termqnumbers,termlabels = zip(*sorted(zip(energies,
                                                    termqnumbers,termlabels)))
        self.termlabels = termlabels
        self.termqnumbers = termqnumbers
        self.termenergies = energies 
        
        return termlabels,energies 
        
    def import_matrixelements(self,filename):
        """Reads in csv file with rows consisting of a lower term label, 
           upper term label, and reduced matrix element in [e a_0]"""
        lowerterms = []
        upperterms = []
        matrixelements = []
        with open(filename,'r') as file:
            csv_read = csv.reader(file,delimiter=',')
            for line in csv_read:
                lowerterms.append(line[0])
                upperterms.append(line[1])
                matrixelements.append(float(line[2]))
        file.close()
        lowerterms,upperterms = check_transition_order(lowerterms,upperterms,
                                             self.termlabels,self.termenergies)
        self.lowertransitionterms = lowerterms
        self.uppertransitionterms = upperterms
        self.matrixelements = matrixelements
        return lowerterms,upperterms,matrixelements
   
    def populate_transitiontable(self):
        """Populate a data table which has term energies on diagonal,
           matrix element of upper right, and transition wavelength on
           lower left"""
        numterm = len(self.termlabels)
        self.numterms = numterm
        transtable = np.zeros((numterm,numterm))
        for (i,j) in np.ndindex(numterm,numterm):
            if i==j:
                transtable[i,j] = self.termenergies[i]
            elif i>j:
                transtable[i,j] = wavelength(self.termenergies[i]
                                                    ,self.termenergies[j])
            elif j>i:
                for (term1,term2,matrixel) in zip(self.lowertransitionterms,
                                self.uppertransitionterms, self.matrixelements):
                    if term1==self.termlabels[i] and term2==self.termlabels[j]:
                        transtable[i,j] = matrixel
        self.transtable = transtable
        
        return transtable
        
    def check_parameters(self):
        """Check if various parameters have been defined, if not, raise 
           exception"""
        variables = ('.Lpow','.Lwaist','.lam','.nuc_I')
        for x,name in zip([self.Lpow,self.Lwaist,self.lam,self.nuc_I],['.Lpow','.Lwaist','.lam','.nuc_I']):
            if x == None:
                errmsg = 'No value given for \'self'+name+'\''
                raise ValueError(errmsg)
        if self.transtable.any() == None:
            self.populate_transitiontable();
            print ('Transition Table was empty: populated with current data.')
                
    def substate_starkshift(self,term,F,mF,field=None):
        """Calculates the Stark shift for a specific (F,mF) state. Returns both 
            the total shift, and the contributions by both term
            and individual state (sorted by size of shift).
            Output:
            totalshift : cumulative stark shift for state, in [Joules]
            reltermcontributions : dictionary with rel term contributions
            relstatecontributions : dictionary with rel state contributions
            """
        self.check_parameters()
        check_state(term,self.nuc_I,F,mF)
        istate = statedic(term,F,mF) #build dictionary for initial state
        index = self.termlabels.index(term) #indicates the index for the corresponding term of interest. 
        I=self.nuc_I
        
            
        Fi = istate['F']
        Ji = istate['J']
        mFi = istate['mF']
        
        #building the vector tensor
        if np.any(field) == None:
            field = electric_field(self.Lwaist,self.Lpow)        
        
        Efield= -(field/2)**2/h

        K=self.rank_vector(self.pol, self.A, Ji) # Determines the kind of contribution for the polarizability
        
        for K0 in K:
            rankfactor=self.rank_factor(istate, I, K0)
            geofactor = self.geometry_factor(self.A, self.pol)
            termcontribution=self.state_contribution(istate,I,K0,index)
            termcontribution['values']=termcontribution['values']*rankfactor
            
        
        totalshift = Efield*sum(termcontribution['values'])
        
        reltermcontribution = sort_and_normalize(termcontribution)

        return totalshift, termcontribution
        
    
    def substate_starkshift_uK(self,term,F,mF):
        '''Calls 'substate_starkshift' method but discards discard all output 
            except the overall shift, which is then scaled to uK'''
        shift , discard1, discard2 = self.substate_starkshift(term,F,mF)
        return shift/k*1e6
        
    def substate_starkshift_MHz(self,term,F,mF):
        '''Calls 'substate_starkshift' method but discards discard all output 
            except the overall shift, which is then scaled to MHz'''
        shift , discard1, discard2 = self.substate_starkshift(term,F,mF)
        return shift/hbar/2/pi*1e-6
    
    def manifold_starkshift(self,term,F):
        ''' Return Stark shift [in MHz] for all substates in a given manifold
            output:
            mF_list: array of mF values
            shifts : array of corresponding shift in MHz
            labels : list of state labels'''
        mF_list = np.arange(-F,F+1).tolist()
        shifts = []
        labels = []
        for mF in mF_list:
            labels.append(statedic(term,F,mF)['label'])
            shifts.append(self.substate_starkshift_MHz(term,F,mF))
        return mF_list,shifts,labels
        
    def term_starkshift(self,term):
        '''Return Stark shift [in MHz] for all substates in a given term
            output:
            F_list : list of F values
            mF_list: lists of mF values
            shifts : lists of corresponding shift in MHz
            labels : list of state labels'''
        [S,L,J] = parse_term(term)
        I=self.nuc_I
        F_list = np.arange(np.abs(J-I),J+I+1).tolist()
        mF_list = []
        shift_list = []
        labels = []
        for F in F_list:
            m,shifts,label = self.manifold_starkshift(term,F)
            mF_list.append(m)
            shift_list.append(shifts)
            labels.append(label)
        return F_list,mF_list,shift_list,labels
    
    def rank_factor(self,istate,I,K):
        '''Calculates coeficient for the tensorial part in the second order reduced polarizability.
            Inputs: dictionary for initial and final state, nuclear spin, and polarization component. 
            Returned coefficient is an amplitude.'''
            
        Fi = istate['F']
        Ji = istate['J']
        mFi = istate['mF']
        factor1 = 0
        factor2 = 0
        factor3 = 0
        
        if K==0:
            factor1 = 2
            factor2 = 1/(3)*1/(2*Ji+1)
            factor3 = 1
        elif K==1:
            factor1 = mFi/Fi
            factor2 = (-1)**(Ji+Fi+I+1)*Wigner6j(Fi,Ji,I,Ji,Fi,1)*np.sqrt(Fi*(2*Fi+1)*(2*Ji+1)*(Ji+1)/(Ji*(Fi+1)) )
            factor3 = np.sqrt(24*Ji/((Ji+1)*(2*Ji+1)))
        elif K==2:
            factor1 = (3*mFi-Fi*(Fi+1))/(Fi*(2*Fi-1))
            factor2 = (-1)**(Ji+Fi+I)*Wigner6j(Fi,Ji,I,Ji,Fi,2)*np.sqrt(Fi*(2*Fi-1)*(2*Fi+1)/((2*Fi+3)*(Fi+1)))*np.sqrt((2*Ji+3)*(2*Ji+1)*(Ji+1)/(Ji*(2*Ji-1)))
            factor3 = np.sqrt(40*Ji*(2*Ji-1)/(3*(Ji+1)*(2*Ji+3)*(2*Ji+1)))
        
        return factor1*factor2*factor3

        
    def state_contribution(self,istate,I,k,index):
        '''

              Inputs:
            istate : dictionary type with quantum numbers of initial state
            targetterm : parseable string labeling target term
            I : nuclear spin
            k: rank tensor contribution
            index : the index of the initial state

              Outputs:
            total : float, sum of all contributions for the term
            contributions : dictionary with state labels and contributions 
            '''
        au = a0*ec #atomic units for reduced matrix elements
        planck=hbar*(2*pi)
        Fi = istate['F']
        Ji = istate['J']
        mFi = istate['mF']
        K0=k
        #print(K0,'rank of contribution')
        termcontribution = {'labels':[],'values':np.array([])}
        #print('upper')

                        #sum over all couplings to terms of higher energy:
        freqL = lam_to_freq(self.lam) #current laser wavelength
       
        for j in np.arange(index+1,self.numterms):

            matrixelement =  self.transtable[index,j]
            
            fterm=self.termlabels[j]
            [Sf,Lf,Jf] = parse_term(fterm)


               #for non-zero matrix element, proceed to calc shift
            if matrixelement != 0 and j!=index:
                wavelength = self.transtable[j,index] 
                freqtrans = lam_to_freq(wavelength) + self.det
                freq_factor = 1/h*(1/(freqtrans+freqL)+(-1)**K0/(freqtrans-freqL))
                
                if K0==0:
                    prefactor=1
                else:
                    prefactor = (-1)**(Ji+Jf)*Wigner6j(Ji,K0,Ji,1,Jf,1)               
                
                term_contribution = (matrixelement*au)**2*freq_factor*prefactor
                termcontribution['labels'].append(self.termlabels[j])
                termcontribution['values']=np.append(termcontribution['values'],term_contribution)
                #print(c/freqtrans*1e9,freqtrans*1e-12,freqL*1e-12,term_contribution)

            else:
                termcontribution['labels'].append(self.termlabels[j])
                termcontribution['values']=np.append(termcontribution['values'],0)

                
                        #sum over all couplings to terms of lower energy:
        #print('lower')
  
        for i in np.arange(0,index):
            matrixelement =  self.transtable[i,index]
  
            fterm=self.termlabels[i]
            [Sf,Lf,Jf] = parse_term(fterm)

                #for non-zero matrix element, proceed to calc the energy level contribution from the table transtable
            if matrixelement != 0:
                wavelength = self.transtable[index,i]
                freqtrans = lam_to_freq(wavelength)+self.det
                freq_factor = 1/h*(1/(freqtrans+freqL)+(-1)**K0/(freqtrans-freqL))
                
                if K0==0:
                    prefactor=1
                else:
                    prefactor = (-1)**(Ji+Jf)*Wigner6j(Ji,K0,Ji,1,Jf,1)
                
                term_contribution = (matrixelement*au)**2*freq_factor*prefactor
                termcontribution['labels'].append(self.termlabels[i])
                termcontribution['values']=np.append(termcontribution['values'],term_contribution)
                #print(term_contribution)

            else:
                termcontribution['labels'].append(self.termlabels[i])
                termcontribution['values']=np.append(termcontribution['values'],0)
            #print(term_contribution)
        return termcontribution       


    def geometry_factor(self,A,pol):
        return 1

    def rank_vector(self,pol,A,J):
        K=[]
        rank_contribution= []
        if pol[2] != 0:
            K.append(0)
            rank_contribution.append('Scalar Contribution')

        if np.abs(A-1)<0.1:
            
            K.append(1)
            rank_contribution.append('Vector Contribution')
            
        if J>1/2:
            K.append(2)
            rank_contribution.append('Tensor Contribution')

        return K

    def rabi2E(self,matrix,rabi,det): #Rabi in rad/s, detuning in Hz, matrix element should include atomic units (au)
        rabi = rabi/(2*pi) #rad/s to Hz
        E = 2*2*pi*hbar*np.sqrt(2*det*rabi)/matrix
        
        return E

    
####### End ACStark Class        
        
####### Auxiliary functions #######

###Methods for handling data structures
def parse_term(term):
    """Takes a term label (sting) of form '5_2p3/2' and returns (S,L,J)
        as per the usual term label conventions. I.e for above example:
        S=(2-1)/2=1/2
        L=p=1
        J=3/2"""
    tri=None
    term=term.split('_')[1]
    for char in term:
        if char.islower() or char.isupper():
            part1 = term.split(char)[0]
            part2 = char.lower()
            part3 = term.split(char)[1]
            break
    S = (float(part1)-1)/2 # total spin
    Ldic = {'s': 0,'p': 1,'d': 2,'f': 3,'g': 4}
    L = Ldic[part2]
    J = float(Fraction(part3))
    
    return [S,L,J]

def statedic(term,F,mF):
    '''Given the term label, F, mF, return a dictionary with all quantum number 
       describing the state'''
    [S,L,J] = parse_term(term)
    label = term + '_%s_%s'%(Fraction(F),Fraction(mF))
    if S==1:
        print('Singlet state')
    elif S==3:
        print('Triplet state')
    
    return {'term':term,'label':label,'S':S,'L':L,'J':J,'F':F,'mF':mF}
    
def check_transition_order(lowerterms,upperterms,termlabels,energies):
    """Here we ensure that the order of the labels in the matrix elements is correct (from lower energy to upper). Returns corrected lists"""
    newlowerlist = []
    newupperlist = []
    for lower,upper in zip(lowerterms,upperterms):
        lower_energy = energies[termlabels.index(lower)]
        upper_energy = energies[termlabels.index(upper)]
        if upper_energy>lower_energy:
            newlowerlist.append(lower)
            newupperlist.append(upper)
        elif upper_energy<lower_energy:
            newlowerlist.append(upper)
            newupperlist.append(lower)
            print ('Swapped terms', lower,'and',upper)
    return newlowerlist,newupperlist


def sort_and_normalize(contributions):
    '''Take a dictionary which has a list of labels and array of values. Sorts and 
        renormalizes value relative to largest element. This method is ugly but 
        don't see a simple way around since I have mixed dictionary,list,and 
        numpy array types all together.'''
    labels = contributions['labels']
    values = contributions['values']
    refvalues = np.abs(values)
    values = values/np.max(refvalues)
    refvalues,labels,values = zip(*sorted(zip(refvalues.tolist(),
                                          labels,values.tolist()),reverse=True))
    return {'labels':labels,'values':values}
    
    
###Methods for error handling
def check_state(term,I,F,mF):
    '''Check that the term, I, F, and mF are a quantum mechanically valid 
        combination'''
    [S,L,J] = parse_term(term)
    Fmax = J + I
    Fmin = np.abs(J-I)
    if F>Fmax or F<Fmin:
        err1 = 'F = %0.1f out of bounds, must lie between %0.1f '%(F,Fmin)
        err2 = 'and %0.1f for a '%(Fmax)+term+' term with J=%0.1f,I=%0.1f'%(J,I)
        raise ValueError(err1+err2)
    if np.abs(mF)>F:
        err1 = 'mF = %0.1f out of bounds, must lie between -%0.1f and'%(mF,F)
        err2 = ' %0.1f for an F=%0.1f manifold'%(F,F)
        raise ValueError(err1+err2)
    if F not in np.arange(Fmin,Fmax+1):
        raise ValueError('%0.1f is invalid F value'%(F))
    if mF not in np.arange(-F,F+1):
        raise ValueError('%0.1f is invalid mF value'%(mF))
        
###General Methods for basic calculations
def wavelength(energy1,energy2):
    """Takes two energies [cm^-1] and returns transition wavelength [m]"""
    deltaE = np.abs(energy2-energy1)
    return 1e-2/deltaE
    
def lam_to_freq(wavelength):
    '''Take wavelength [m] and returns frequency [Hz]'''
    return c/wavelength
    
def laser_intensity(waist,power):
    '''Takes beam waist and power, returns peak  intensity'''
    return  2*power/(pi*waist**2)
    
def electric_field(waist,power):
    '''Takes beam waist and power, returns electric field amplitude'''
    inten = laser_intensity(waist,power)
    return np.sqrt(2*inten/(c*epsilon_0))
    
def pol_projection(q,polvec):
    '''For a normalized polarization vector, returns the projection factor onto
        the q-component. (amplitude squared)'''
    if q==-1:
        return (polvec[0])**2 
    elif q==0:
        return (polvec[1])**2
    elif q==1:
        return (polvec[2])**2
    else:
        raise ValueError('%f invalid q value.'%(q))


