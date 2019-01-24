# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:34:04 2019

@author: Max
"""

import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

def time_evol(psi0,H,t0,T,dt,display_progress=True,plot=False):
    '''
    Simulate the unitary time evolution for a given Hamiltonian over a certain
    timespan. Then return output data and plot if desired.
    
    Arguments:
        psi0 - Initial state of system.
        H - The Hamiltonian to be applied in the time evolution unitary.
        dt - The increment between calculated points in the simulation.
        t0 - Start time of simulation
        T - End time of simulation
        display_progress - Choose whether to show progress of simulation whilst
                           it runs.
        plot - Boolean to set whether a graph is plotted.
    
    '''
    times = np.arange(t0,T,dt)
    
    
    return(qt.mesolve(H,psi0,times,progress_bar=display_progress))
    
    
    
ham = qt.sigmax()

up = qt.basis(2,0)
down = qt.basis(2,1)

t0=0
T=10
dt=0.01

res = time_evol(up,ham,t0,T,dt,False)
states = res.states

up_states = []
down_states = []

for i in range(len(states)):
    up_states.append(states[i][0])
    down_states.append(states[i][1])
    
up_states,down_states = np.array(up_states),np.array(down_states)

up_prob,down_prob = np.absolute(up_states).flatten(),np.absolute(down_states).flatten()
print(up_prob)
plt.plot(res.times,up_prob)
plt.show()

    
    