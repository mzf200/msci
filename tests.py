import qutip as qt
import copy as cp 
J=3
nu=1
nphonon=5
ldps=[1,1,2]
rabi=2

basis_m1 =   qt.basis(4,0)
basis_0 =    qt.basis(4,1)
basis_p1 =   qt.basis(4,2)
basis_0tag = qt.basis(4,3)


sigmaz_fl = [[-1,0,0,0],
             [0,0,0,0],
             [0,0,1,0],
             [0,0,0,0]]

sigmaz_fl = qt.Qobj(sigmaz_fl)

# =============================================================================
# a_dag = qt.create(7)
# a =  qt.destroy(7)
# pol = (a_dag-a)
# H_test = qt.tensor(basis_0 *basis_m1.trans(), pol)
# =============================================================================

#print(H_test)



eye_list = [qt.qeye(4) for j in range(J)]
#print(eye_list)
  #  H_spsp = 0                                      
for j in range(J-1):
   for k in range(j+1,J): #iteration to avoid double summation or summation over the same index.
       print(j,k)
       H_temp = cp.deepcopy(eye_list)
       print('1=',H_temp)
       H_temp[j] = sigmaz_fl
       H_temp[k] = sigmaz_fl
       print('2=',H_temp)

       H_temp = H_temp
       print('3=',H_temp)
      # print('H_temp=',H_temp)
       H_temp = qt.tensor(H_temp)*ldps[j]*ldps[k]
#       print(H_temp)
           
       if k ==j+1:
           H_spsp = H_temp
       else:
           for l,m in enumerate(H_spsp):
               print('l=',H_temp,'m=',m) 
               H_spsp += H_temp 
               
       print(H_spsp) 																	

# =============================================================================
# eye_list = [qt.qeye(4) for j in range(J)]	
# print('eye0',eye_list)
# H_temp = cp.deepcopy(eye_list)
# H_temp[0] = sigmaz_fl
# H_temp[1] = sigmaz_fl
# print(H_temp)
# H_temp[0] += eye_list[0]
# print(H_temp)
# print('eye', eye_list)
# =============================================================================


#H_spinspin = nu * H_spsp
#print(H_spinspin)       

       
       
       
       
       
       