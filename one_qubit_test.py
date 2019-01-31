#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:37:02 2019

@author: avinashmocherla
"""
import qutip as qt
import numpy as np


np.set_printoptions(threshold=np.nan)
def dressed_Hamiltonian(rabi):
    """
    Returns the dressed state Hamiltonian. In the form: 
        H = H_mot + H_spinspin + H_int
    
    Args:
        J: Number of ions in trap 
        nu: trap frequency 
        nphonon: maximum number of phonons in system
        ldp: Lamb-Dicke Parameter
        rabi: driving field rabi frequency
        
    """
    

    """Dressed Basis States"""
    basis_m1 = qt.tensor(qt.basis(4,0))
    
    basis_0 = qt.tensor(qt.basis(4,1))

    basis_p1 = qt.tensor(qt.basis(4,2))

    basis_0tag = qt.tensor(qt.basis(4,3))
    

    

#    c1 = qt.tensor((basis_p1 * basis_0.trans()), qt.qeye(4),qt.qeye(nphonon) )
# 
#    
#    c2 = qt.tensor((basis_0 * basis_m1.trans()), qt.qeye(4), qt.qeye(nphonon))
#
#    
#    c3 = qt.tensor(qt.qeye(4), (basis_p1 * basis_0.trans()), qt.qeye(nphonon))
#
#    
#    c4 = qt.tensor(qt.qeye(4), (basis_0 * basis_m1.trans()), qt.qeye(nphonon))
    
    
    c1 = qt.tensor((basis_p1 * basis_0.trans()))
 
    
    c2 = qt.tensor((basis_0 * basis_m1.trans()))


#    t1 = (c1 + c2) * (ldp*(a_dag - a)).expm()   + (c1.trans() + c1.trans()) * (ldp*(a_dag - a)).expm()
#    
#    t2 = (c3 + c4) * (ldp*(a_dag - a)).expm() + (c3.trans() + c4.trans()) * (ldp*(a_dag - a)).expm()
    
    t1 = (c1 + c2)  + (c1.trans() + c2.trans())
    

    H_int =  t1

    return H_int



def time_evol(H, psi0, psif, t0, T, dt, display_progress=True, plot=True):
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
   
    results = qt.sesolve(H,psi0,times,options = optns)
    states = results.states
    
    fidelities = []
    
    for i in states:
            fidelities.append(qt.fidelity(i,psif))
    
    if plot==True:

        plt.plot(times, fidelities)
        plt.show()
        

   # print(optns)
    return(fidelities)
#    
##%%  
#        
#"""System Params"""    
#rabi= np.sqrt(2)*2*np.pi*45.4e3
#
#"""basis definitions"""   
#basis_m1 =   qt.basis(4,0)
#basis_0 =    qt.basis(4,1)
#basis_p1 =   qt.basis(4,2)
#basis_0tag = qt.basis(4,3)
#
#"""dressed state basis definitions """
#up = 0.5 * basis_p1 + 0.5 * basis_m1 + (1/np.sqrt(2))* basis_0
#down = 0.5 * basis_p1 + 0.5 * basis_m1 - (1/np.sqrt(2))* basis_0
#D = (1/np.sqrt(2)) * basis_p1 - (1/np.sqrt(2))*basis_m1
#
#    
#
#
#"""TimeEvol Params"""    
#t0=0
#T=0.4
#dt = 1e-3
#psi0 = up
#psif = psi0
#
#
#H = dressed_Hamiltonian(rabi)
#
#fidelities = time_evol(H, psi0, psif, t0,T,dt, True, True)
#

