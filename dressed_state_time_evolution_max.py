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
      
    basis_m1 =   qt.basis(4,0)
    basis_0 =    qt.basis(4,1)
    basis_p1 =   qt.basis(4,2)
    basis_0tag = qt.basis(4,3)

    eye_list = [qt.qeye(4) for j in range(J)]

    sigmaz_fl = [[-1,0,0,0],
                 [0,0,0,0],
                 [0,0,1,0],
                 [0,0,0,0]]
    
    sigmaz_fl = qt.Qobj(sigmaz_fl)
  
    H_int_J = [0 for j in range(J)]
    
    eye_temp = cp.deepcopy(eye_list)
    eye_temp.append(qt.create(nphonon))
    a_dag = qt.tensor(eye_temp)
    
    eye_temp = cp.deepcopy(eye_list)
    eye_temp.append(qt.destroy(nphonon))
    a =     qt.tensor(eye_temp)

#     sigmaz = qt.tensor(qt.sigmaz(), qt.qeye(2), qt.qeye(nphonon))
#     sigmax_2 = qt.tensor(qt.qeye(2), qt.sigmaz(), qt.qeye(nphonon))
    
    H_mot = nu * a_dag * a
    
    for j in range(J):
        pol,pol_dag = ldps[j]*(a_dag-a),-ldps[j]*(a_dag-a)
        pol = pol.expm()
        pol_dag = pol_dag.expm()
        eye_temp = cp.deepcopy(eye_list)
        
        d1 = qt.tensor(basis_p1 * basis_0.trans(),qt.qeye(nphonon))
        d2 = qt.tensor(basis_0 * basis_m1.trans(),qt.qeye(nphonon))
        d3 = qt.tensor(basis_p1 * basis_0.trans(),qt.qeye(nphonon))
        d4 = qt.tensor(basis_m1* basis_0.trans(),qt.qeye(nphonon))
        
        H_int_j = (rabi/2) *(d1*pol +
        d2*pol +
        d3*pol_dag +
        d4*pol_dag)
             
        H_int_J[j]=H_int_j
    
    H_int_J = qt.tensor(H_int_J)
                                              

  #  H_spsp = 0                                      
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
                                                 
    return(H_int_J, H_mot, H_spinspin)
	
#%%

    
    
# J, nu, nphonon, ldps, rabi
print (dressed_Hamiltonian(J=2,nu=1,nphonon=3,ldps=[1,1],rabi=1))
						
							 
							 
