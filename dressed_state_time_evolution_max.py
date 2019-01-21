#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:37:02 2019

@author: avinashmocherla
"""
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
# q = qt.Qobj(([0],[1]))



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
      
    basis_m1 =   qt.basis(4,0)
    basis_0 =    qt.basis(4,1)
    basis_p1 =   qt.basis(4,2)
    basis_0tag = qt.basis(4,3)

    sigmaz_fl = [[-1,0,0,0],
                 [0,0,0,0],
                 [0,0,1,0],
                 [0,0,0,0]]
    
    sigmaz_fl = qt.Qobj(sigmaz_fl)
  
    H_int_J = [0 for j in range(J)]
    
    a_dag = qt.create(nphonon)
    a =     qt.destroy(nphonon)

#     sigmaz = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(nphonon))
#     sigmax_2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(nphonon))
    
    H_mot = nu * a * a_dag
    
    for j in range(J):
		H_int_j = (rabi/2) * np.exp(ldp[j]*(a_dag - a)) * (qt.tensor(basis_p1 * basis_0.trans(), qt.qeye(4)) +
		qt.tensor(basis_0 * basis_m1.trans(), qt.qeye(4)) + (rabi/2) * np.exp(ldp[j]*(a_dag - a)) * np.exp(ldp*(a_dag - a)) *
											(qt.tensor(basis_p1 * basis_0.trans(), qt.qeye(4))  +
											 qt.tensor(basis_0 * basis_m1.trans(), qt.qeye(4))) +
		(rabi/2) * np.exp(ldp[j]*(a_dag - a)) * np.exp(ldp*(a_dag - a)) *
											(qt.tensor(basis_p1 * basis_0.trans(), qt.qeye(4))  +
											 qt.tensor(basis_0 * basis_m1.trans(), qt.qeye(4))) +
		(rabi/2) * np.exp(ldp[j]*(a_dag - a)) * np.exp(ldp*(a_dag - a)) *
											(qt.tensor(basis_p1 * basis_0.trans(), qt.qeye(4)) +
											 qt.tensor(basis_0 * basis_m1.trans(), qt.qeye(4)))
		H_int_J[j]=H_int_j
    
    H_int_J = qt.tensor(H_int_J)
                                              
    eye_list = [qt.eye(4) for j in J]
    H_spsp = 0                                      
	for j in range(J):
        for k in range(j+1,J): #iteration to avoid double summation or summation over the same index.
            H_temp = ldps[j]*ldps[k]*eye_list
            H_temp[j] = sigmaz_fl
            H_temp[k] = sigmaz_fl
            H_temp = qt.tensor(H_temp)
            H_spsp += H_temp
																					
		
    H_spinspin = nu * H_spsp
                                                 
    return(H_int_J, H_mot, H_spinspin)
	
print (dressed_Hamiltonian(3,1,5,[1,1,2],2))
						
							 
							 
