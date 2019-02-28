import numpy as np 
import qutip as qt
import matplotlib.pyplot as plt
from qutip import Bloch

from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

import colorednoise as cn
beta = 2 # the exponent
samples = 2**18 # number of samples to generate
y = cn.powerlaw_psd_gaussian(beta, samples)


omega_plus = 1 
omega_minus = 1
nphonon = 5 
ld = 1
nu = 1

sigma_x = qt.sigmax()
sigma_y = qt.sigmay()
sigma_z = qt.sigmaz()
identity = qt.qeye(2)

#%%
    
def H_laser_coeff(t,args):
    return np.cos(2*t)

def H1_coeff(t,args):    
    return  0.1 * (np.cos(0.001* t))

def pink_noise_coeff(t, args):
    return 0.1 * np.abs(y[int(t)])
 
    
def pulse(initial_state, duration, rabi, error = False):
    t0 = omega_plus + omega_minus #transition frequency
    
    H0 = ((1/2) * t0 * sigma_z)  
    
    H_laser = rabi * sigma_x
    
    if error == False:  
    
        H_final = [H0, [H_laser, H_laser_coeff]]
    
    if error == True: 
        H1 = (1/2) * sigma_z
        
        H_final = [H0, [H1, pink_noise_coeff], [H_laser, H_laser_coeff]]
  
   
   
    optns = qt.Options()
    
    optns.nsteps = 1000
    
    states = qt.sesolve(H_final, initial_state, np.linspace(0,duration,100), options = optns).states

   
    final_state = states[-1]

   
    return(final_state, states)



def free_ev(initial_state, duration, error = False):

    t0 = omega_plus + omega_minus #transition frequency
    
#    H_couple = nu * ld * (a + a_dag) * sigma_z
    
    H0 = ((1/2) * t0 * sigma_z)  
    
    if error == False:
    
        H_final = H0
        
    if error == True: 
     
        H1 = (1/2) * sigma_z
        
        H_final = [H0, [H1, pink_noise_coeff]]
    
    optns = qt.Options()
    
    optns.nsteps = 100000
    
    states = qt.mesolve(H_final, rho0 = initial_state, tlist= np.linspace(0,duration,100),  options = optns).states
    
    
    

    final_state = states[-1]
    
    return(final_state, states)

def y_pulse(initial_state):
    
    return(1j * sigma_y * initial_state)
    
#%%
 
def decoherence(rabi):
    """
    performs a free induction decay on a superposition qubit state 
    as a function of time to measure coherence 
    """
    initial_state = (1/np.sqrt(2)) * (qt.basis(2,0) + 1j * qt.basis(2,1))
    
#    p_initial = initial_state * initial_state.dag()
    
    duration =  10 * 2*np.pi / rabi
    
    state = free_ev(initial_state, duration, error = True)    
    
    state1 = pulse(state[0], (np.pi / rabi),0.01, error = False)
    
    print(state1)
    
#    print(state[1])
    plt.plot(state1[0])
    
#    W = np.zeros(len(state[1]))
 
#    for i in range(len(state[1])):    
#        p_pm = state[1][i].full()[0][1]
#        p_pm = np.sqrt(np.conj(p_pm)*p_pm)
#       
#        W[i] = p_pm
        
        
#    W = W/np.abs((-1j / 2))
    

    
    plt.show()
    
#    return(state[1])
    
#    duration = (2*np.pi / 2 * rabi) 
#      
#    state = pulse(initial_state, duration, rabi, error = False)
#    
#    state1 = free_ev(state[0], duration * 2)
#    
#    state2 = pulse(state1[0], duration, rabi, error = False)
#    
#    probabilities =[]
#    for i in state[1]:
#          overlap = final_state.overlap(i)
#          metric = np.dot(overlap, overlap.conj())
#          probabilities.append(metric)
#    
#    plt.plot(probabilities)
#    plt.show()
    
    return(state[1])
    
def coherence(rabi):
    
    initial_state = qt.basis(2,0) #Set initial State of qubit to 0
    
    final_state = qt.basis(2,0)
    
    duration = 10*(2*np.pi / rabi) 
      
    state = pulse(initial_state, duration, rabi, error = True)
   
    
    
    probabilities =[]
    for i in state[1]:
          overlap = final_state.overlap(i)
          metric = np.dot(overlap, overlap.conj())
          probabilities.append(metric)
    
    plt.plot(probabilities)
    plt.show()
    
#    return(state[1])
 #%%   
def decouple(rabi, method):
    """
    Perform a single-qubit rotation of pi radians (from 0 to 1)
    with dynamical decoupling included. 
    
    If the dynamical decoupling procedure is perfectly successful,
    then the overlap after the final pi/2 rotation with the 1 state 
    equals 1. If there has been a lot of uncorrected dephasing 
    (the qubit state has not been aligned perfectly along the y axis), 
    then this will reduce to 0.5). 
    
    Args: 
        rabi: rabi frequency of applied pi pulse
        method: dynamical decoupling procedure to use 
        (varies the number and duration of pulses) (default = SE)
    
    """
    initial_state = 1/sqrt(2) * (qt.basis(2,0) + 1j*qt.basis(2,1)) #Set initial State of qubit to 0
    
    final_state = qt.basis(2,1) #Set final comparison state to test how good procedure was
    
#    state = pulse(initial_state, -np.pi/(2*rabi), rabi) # Apply first pi/2 pulse to create superposition state
    
    if method == 'None': 
        duration = (2*np.pi / rabi) / 2 
        
        state1 = free_ev(initial_state, duration, error=True)
        
        state_final = free_ev(state1[0], duration, error=True)
        
        states = [*[initial_state for i in range(10)], *state1[1], *[state1[0] for i in range(10)],  *state_final[1]]
        
    if method == "SE":
    
        duration = (2*np.pi / rabi) / 2 #split the free evolution time 
       
        state1 = free_ev(initial_state, duration, error = True)
        
        state2 = y_pulse(state1[0])
        
      
        print(state2)
        state_final = free_ev(state2, duration, error = True)
        print(state_final[1])
        
        states = [*[initial_state for i in range(10)], *state1[1], *[state1[0] for i in range(10)], *[state2 for i in range(10)], *state_final[1] ]
        
    
    final = pulse(state_final[0], -np.pi/(2*rabi), rabi, error = False)

    overlap = final_state.overlap(final[0])
  
    metric = np.dot(overlap, overlap.conj())

    for i in final[1]:
        states.append(i)
   
    return(states, metric)

#%%
def low_freq_noise_results(rabi):
    """
    Compares None, SE and CPMG decoupling for low frequency noise (AC errors)
    """
    no_correction = decouple(rabi, 'None')[1]
    SE = decouple(rabi, 'SE')[1]

    return(no_correction, SE) 
    
def high_freq_noise_results(rabi):
    """
    Compares None, SE and CPMG decoupling for high frequency noise
    """
    no_correction = decouple(rabi, 'None')[1]
    SE = decouple(rabi, 'SE')[1]
    
    return(no_correction, SE) 
    
def pink_noise_results():
    """
    Compares None, SE and CPMG decoupling for 1/f magnetic noise
    """
    pass

#%%
def animate(states):

    fig = figure()
    ax = Axes3D(fig,azim=-40,elev=30)
    sphere = Bloch(axes=ax)
    
    def animate(i):
        sphere.clear()
    
        sphere.add_states([states[i]])
        sphere.make_sphere()
        return ax
    
    def init():
        sphere.vector_color = ['r']
        return ax
    
    ani = animation.FuncAnimation(fig, animate, np.arange(len(states)),
                                init_func=init, blit=False, repeat=True)
    ani.save('decoherence.mp4', fps=20)


