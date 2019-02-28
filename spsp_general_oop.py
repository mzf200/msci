
"""
Created on Fri Jan 11 15:37:02 2019
@author: avinashmocherla
"""
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.nan)


#%%

class trapped_ions:
    """
    Returns the dressed state Hamiltonian. In the form: 
        H = H_mot + H_spinspin + H_int
    
    Args:
        J: Number of ions in trap 
        basis: 0 (bare) or 1(dressed) is the working basis 
        nu: trap frequency 
        nphonon: maximum number of phonons in system
        ldp: Lamb-Dicke Parameter
        rabi: driving field rabi frequency     
        spinspin: bool, include spinspin hamiltonian
    """
    
    def __init__(self,J, basis = 0, nu=1, nphonon=1, ldps=[1], rabi=1, spinspin = True):

        self.__J
        self.__basis = basis
        self.__nu = nu
        self.__nphonon = nphonon
        self.__ldps = ldps
        self.__rabi = rabi
        self.__spinspin = spinspin
        self.__eigs = None
    
        self.__H = qt.Qobj( np.zeros(shape=(4**J,4**J)), dims = [[4 for j in range(J)],[4 for j in range(J)]],)
        """Dressed Basis States"""
        self.__basis_m1 = qt.tensor(qt.basis(4,0))
    
        self.__basis_0 = qt.tensor(qt.basis(4,1))

        self.__basis_p1 = qt.tensor(qt.basis(4,2))

        self.__basis_0tag = qt.tensor(qt.basis(4,3))
        
        self.__dark = 1/(np.sqrt(2))*(self.__basis_p1 - self.__basis_m1)
        
    
        self.__sigma_z = [[-1,0,0,0],
                          [0,0,0,0],
                          [0,0,1,0],
                          [0,0,0,0]]
    
        self.__sigma_z = qt.Qobj(self.__sigma_z)
        """Motional mode operators"""
        self.__a = qt.tensor(qt.qeye(4), qt.qeye(4))

        self.__a_dag = qt.tensor(qt.qeye(4), qt.qeye(4))

        """Spin Spin Operators"""
        self.__sigmaz_1 = qt.tensor(self.__sigma_z, qt.qeye(4))
    
        self.__sigmaz_2 = qt.tensor(qt.qeye(4), self.__sigma_z)
    
        """More Spin-Spin Operators""" 
        
        
        c1 = self.__basis_p1 * self.__basis_0.trans()
        c2 = self.__basis_0 * self.__basis_m1.trans()
        t1 = (c1 + c2)  + (c1.trans() + c2.trans())
        H_int = [t1 for j in range(J)]
        self.__H += 0.5* self.__rabi*qt.tensor(H_int)
# =============================================================================
#         if J == 1: 
#              c1 = self.__basis_p1 * self.__basis_0.trans()
#              c2 = self.__basis_0 * self.__basis_m1.trans()
#              t1 = (c1 + c2)  + (c1.trans() + c2.trans())
#              print(t1)
#              H_int = t1
#              self.__H += 0.5 * self.__rabi * H_int
#     
#         elif J == 2: 
#              c1 = self.__basis_p1 * self.__basis_0.trans()
#              c2 = self.__basis_0 * self.__basis_m1.trans()
#              t1 = (c1 + c2)  + (c1.trans() + c2.trans())
#         
#              H_int = qt.tensor(t1,t1)
#             
#              self.__H += 0.5 * rabi * H_int
# =============================================================================
         

        if spinspin == True:
            for j in range(J-1):
                for k in range(j+1,J):
                    spsp = [qt.qeye(4) for l in J]
                    spsp[j],spsp[k] = self.__sigma_z, self.__sigma_z
                    self.__H +=  self.__nu * self.__ldps[j]*self.__ldps[k] * (spsp)
        else: pass 
  
#%%
    def time_evol(self, psi0, psif, t0 = 0, T = 10, dt = 1e-2, plot=True):
        '''
        Simulate the unitary time evolution for a given Hamiltonian over a certain
        timespan. Then return output data and plot if desired.
        
        Arguments:
            
            H - The Hamiltonian to be applied in the time evolution unitary.
            psi0 - Initial state of system.
            psif - Final (desired) state of system against which to check fidelity
            t0 - Start time of simulation
            T - End time of simulation
            dt - Time increment
            
           
            display_progress - Choose whether to show progress of simulation whilst
                               it runs.
            plot - Boolean to set whether a graph is plotted.
        
        '''
        times = np.arange(t0,T,dt)
        
        
        optns = qt.Options()
        
        optns.nsteps = 10000
       
        results = qt.sesolve(self.__H, psi0, times, options = optns,progress_bar=True)
        states = results.states
        
        overlap = []
        
        for i in states:
                ov = i.overlap(psif)
                overlap.append(ov*np.conj(ov))
        
        if plot==True:
    
            plt.plot(times, overlap)
            plt.xlabel('Time/(hbar*nu)')
            plt.ylabel('|<psi(t)|D>|^2')
       #     plt.ylim([0.96,1.04])
            plt.show()
            
    
        return(overlap)
        
#%%
        
    def eig_states(self):
        if self.__eigs != None:
            return self.__eigs
        elif self.__eigs == None:
            return(self.__H.eigenstates())
        
#%%
    def gen_state(self, states = []):
        if type(states) != list:
            st = []
            for j in range(J):
                if states == 'plus':
                    st.append(self.__basis_p1)
                if states == 'minus'or 'm1' or '-':
                    st.append(self.__basis_m1)
                if states == 'zero' or '0':
                    st.append(self.__basis_0)
                if states == 'zero tag' or 'zero_tag' or '0tag':
                    st.append(self.__basis_0tag)
                if states == 'D' or 'Dark' or 'dark':
                    st.append(self.__dark)
            st = qt.tensor(st)
            return st
        
        else:
            st = []
            for i in states:
                if i == 'plus' or 'p1' or '+':
                    st.append(self.__basis_p1)
                if i == 'minus' or 'm1' or '-':
                    st.append(self.__basis_m1)
                if i == 'zero' or '0':
                    st.append(self.__basis_0)
                if i == 'zero tag' or 'zero_tag' or '0tag':
                    st.append(self.__basis_0tag)
                if i == 'D' or 'Dark' or 'dark':
                    st.append(self.__dark)
            st = qt.tensor(st)
            return st
                    
                    
    
#%%  

def one_qubit_eigenbasis_check():
    """
    Find eigenstates of H_int for 1 qubit.
    Evolves each eigenstate to check if stationary, as expected.
    """
    trap = trapped_ions(1, spinspin = False)
    eigs = trap.eig_states()
    print(eigs)
    eigs = eigs[1:]
    for i in range(len(eigs[0])):
        psi0 = eigs[0][i]
        psif = psi0
        overlap = trap.time_evol(psi0, psif)
        
#%%
def two_qubit_eigenbasis_check():
    """
    Find eigenstates of H_int for 2 qubit.
    Evolves each eigenstate to check if stationary, as expected.`
    """
    trap = trapped_ions(2, spinspin = False)
    eigs = trap.eig_states()
    eigs = eigs[1:]
    for i in range(len(eigs[0])):
        psi0 = eigs[0][i]
        psif = psi0
        overlap = trap.time_evol(psi0, psif)
        
#%%
def test1(step):
    """
    Measures the effect of the spin spin interaction on the evolution 
    varying strength interaction 
    
    Args:
        step = number of times to increase spin-spin strength
    """
    
    ldp_strength = np.linspace(0,5,3)
    for i in range(len(ldp_strength)):
        trap = trapped_ions(2, ldp = ldp_strength[i], spinspin = True)
        eigs = trap.eigenstates()
        eigs = eigs[1:]
        print(eigs[0])
        
        for j in range(len(eigs[0])):
            print(i,j)
            psi0 = eigs[0][j]
            psif = psi0 
            overlap = trap.time_evol(psi0, psif)
            
#%%


nu_tag = 1
nphonon = 1
nu = 2*np.pi*459.34E3
ldps=[0.0041,-0.0041]
rabi = 12.6E9/nu
J=2

trap = trapped_ions(J, basis = 0, nu=nu_tag, nphonon= nphonon,
                          ldp=ldps[1], rabi=rabi, spinspin = True)

D = trap.gen_state('D')

psi0 = D#qt.tensor(basis_m1,basis_m1)
psif = psi0         



fid = trap.time_evol(psi0, psif, t0 = 0, T = 100, dt = 1E-4, plot=True)
            
            
            
            
            
            
            
            


