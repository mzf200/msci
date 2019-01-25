#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:37:02 2019
@author: avinashmocherla
"""
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import copy as cp 

# q = qt.Qobj(([0],[1]))

np.set_printoptions(threshold=np.nan)

#%%

def dressed_Hamiltonian(J, nu, nphonon, ldps, rabi):
   
    
    """
    Returns the dressed state Hamiltonian. In the form: 
        H = H_mot + H_spinspin + H_int
    
    Args:
        J: Number of ions in trap 
        nu: trap frequency 
        nphonon: maximum number of phonons in system
        ldps: Lamb-Dicke Parameters
        rabi: driving field rabi frequency
        
    """
    if len(ldps) != J:
      raise Exception('There must be a Lamb-Dicke Parameter for each ion.')
      
    """Initialisation of Basis Vectors"""  
      
    # Define basis vectors of dressed system   
    basis_m1 =   qt.basis(4,0)
    basis_0 =    qt.basis(4,1)
    basis_p1 =   qt.basis(4,2)
    basis_0tag = qt.basis(4,3)


    # Initialise list of ion operators 
    eye_list = [qt.qeye(4) for j in range(J)]
    
    # Define sigma-z operator in dressed basis 
    sigmaz_fl = [[-1,0,0,0],
                 [0,0,0,0],
                 [0,0,1,0],
                 [0,0,0,0]]
    
    sigmaz_fl = qt.Qobj(sigmaz_fl)
  
    # Initialise Hamiltonian 
    H_int_J = 0 
    
    eye_temp = cp.deepcopy(eye_list)
    eye_temp.append(qt.create(nphonon))
    a_dag = qt.tensor(eye_temp)
    
    eye_temp = cp.deepcopy(eye_list)
    eye_temp.append(qt.destroy(nphonon))
    a =     qt.tensor(eye_temp)

    """ Motional Hamiltonian """      

    
    
    
    
    
    
    
    H_mot = nu * a_dag * a
    
    
    """ Internal  Hamiltonian """  
    for j in range(J):
        pol,pol_dag = ldps[j]*(a_dag-a),-ldps[j]*(a_dag-a)
        pol = pol.expm()
        pol_dag = pol_dag.expm()
        
        eye_temp1 = cp.deepcopy(eye_list)
        eye_temp2 = cp.deepcopy(eye_list)
        eye_temp3 = cp.deepcopy(eye_list)
        eye_temp4 = cp.deepcopy(eye_list)
        
        cross_term_1 = basis_p1 * basis_0.trans()
        cross_term_2 = basis_0 * basis_m1.trans()
        cross_term_3 = basis_p1 * basis_0.trans()
        cross_term_4 = basis_m1* basis_0.trans()       
        
        eye_temp1[j] = cross_term_1
        eye_temp2[j] = cross_term_2
        eye_temp3[j] = cross_term_3
        eye_temp4[j] = cross_term_4     
        
      
        
        
        d1 = qt.tensor(eye_temp1)
        d2 = qt.tensor(eye_temp2)
        d3 = qt.tensor(eye_temp3)
        d4 = qt.tensor(eye_temp4)
        
        d1 = qt.tensor(d1, qt.qeye(nphonon))
        d2 = qt.tensor(d2, qt.qeye(nphonon))
        d3 = qt.tensor(d3, qt.qeye(nphonon))
        d4 = qt.tensor(d4, qt.qeye(nphonon))
        
        H_int_j = (rabi/2) *(d1*pol +
        d2*pol +
        d3*pol_dag +
        d4*pol_dag)
        
        H_int_J += H_int_j
             

                                              
     
    """ Spin-Spin Hamiltonian """                           
    for j in range(J-1):
       for k in range(j+1,J): #iteration to avoid double summation or summation over the same index.
           #print(j,k)
           H_temp = eye_list
           H_temp[j] = sigmaz_fl
           H_temp[k] = sigmaz_fl
           H_temp = qt.tensor(H_temp)*ldps[j]*ldps[k]
          # print(H_temp)
           
           if j==0 and k==1:
               H_spsp = H_temp
           else:
               H_spsp += H_temp
           #print(H_spsp) 
																	
		
    H_spinspin = nu * qt.tensor(H_spsp,qt.qeye(nphonon))
    
    eig_up = [0.5*basis_p1 + 0.5*basis_m1 + (1./(np.sqrt(2)))*basis_0\
              for j in range(J)]
  #  print(eig_up)

    eig_dn = [0.5*basis_p1 + 0.5*basis_m1 - (1./(np.sqrt(2)))*basis_0\
              for j in range(J)]
    
    eig_up.append(qt.basis(nphonon,0))
    eig_dn.append(qt.basis(nphonon,0))
    
    eig_up = qt.tensor(eig_up)
    eig_dn = qt.tensor(eig_dn)
    print(eig_up)
                                             
    return (H_int_J + H_spinspin + H_mot)
	
# J, nu, nphonon, ldps, rabi
#print(dressed_Hamiltonian(J=2,nu=1,nphonon=3,ldps=[1,1],rabi=1))
#%%

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
    
#%%  
        
"""System Params"""    
J = 2
nu = 2*np.pi*459.341e3
nphonon = 3
ldps=[0.0041,-0.0041]
rabi= np.sqrt(2)*2*np.pi*45.4e3

"""basis definitions"""   
basis_m1 =   qt.basis(4,0)
basis_0 =    qt.basis(4,1)
basis_p1 =   qt.basis(4,2)
basis_0tag = qt.basis(4,3)

eig_up = [0.5*basis_p1 + 0.5*basis_m1 + (1./(np.sqrt(2)))*basis_0\
              for j in range(2)]
eig_up.append(qt.basis(nphonon,0))
eig_up = qt.tensor(eig_up)
    
D = [(1/np.sqrt(2)) * basis_p1 - (1/np.sqrt(2))*basis_m1 for j in range(J)]
D.append(qt.basis(nphonon,0))
D = qt.tensor(D)
"""TimeEvol Params"""    
t0=0
T=1
dt = 1e-4
psi0 = qt.tensor(qt.basis(4,0), qt.basis(4,0), qt.basis(nphonon,0))
psif = D


H = dressed_Hamiltonian(J,nu,nphonon,ldps, rabi)

fidelities = time_evol(H, psi0, psif, t0,T,dt, True, True)































				