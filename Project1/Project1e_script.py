

#Comparing the cpu time for the special tridiagonal algorithm
#and the general tridiagonal algorithm implemented in the 
#function TriDiagSpes in Diagonal_Solvers.py, with the solver based on LU 
#decomposition in LU_Solver.py. The comparison is performed for increasing number 
#of grid points n.

import numpy as np
from Diagonal_Solvers import *
from LU_Solver import *



#An array with increasing number of grid points n 
#n_array = np.array([5,10])
n_array = np.array([10, 100, 1000, 10**4])

#An array where the last element really puts the LU solver to the test with
#10**5 gridpoints, resulting in a (10**5)*(10**5) matrix: 

#n_array = np.array([10, 10**5])

cpu_LU = np.zeros(np.size(n_array))
cpu_special = np.zeros(np.size(n_array))
cpu_general = np.zeros(np.size(n_array))


#Finding elapsed cpu time and number of flops for each n in n_array
#for both algorithms
i=0
for n in n_array:
    #step size
    h = 1.0/(n+1)
 
    x = np.linspace(h, n*h, n)
    f = np.zeros(n)

    #Our second derivative u'' = f:
    f = h**2*(100*np.exp(-10*x))

    #Calling the function with the special tridiagonal algorithm
    [v,flop_s,cpu_s] = TriDiagSpes(n,f)

    
    a = -1
    b = 2
    c = -1
    #Calling the function with the general tridiagonal algoritm
    
    [v,flop_g,cpu_g] = TriDiag(n,a,b,c,f)

    #Calling the function with the algorithm based on LU decomposition.
    #Need to create the tridiagonal matrix A first
    A = np.zeros((n,n))
    A[0,0] = b
    A[0,1] = c
    A[n-1,n-1] = b
    A[n-1,n-2] = a


    for j in range(1,n-1):
        A[j,j] = b
        A[j, j -1] = a
        A[j,j+1] = c

    [v,cpu_lu] = LU_Solver(A,f,n)
    
    cpu_LU[i] = cpu_lu
    cpu_special[i] = cpu_s
    cpu_general[i] = cpu_g
    #print(A)
    i = i+1

#Table
print('') 
print('Table with  elapsed cpu time for the algorithms:')
print('')
print('                 elapsed cpu time     ')
print('    n      general    special   LU-factorization')               
for i in range(0,np.size(n_array)):
    print("""  %.0E   %.3E  %.3E  %.3E"""
    %(n_array[i], cpu_general[i],cpu_special[i],cpu_LU[i]))
	

#Running the script
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python Project1e_script.py

Table with  elapsed cpu time for the algorithms:

                 elapsed cpu time
    n      general    special   LU-factorization
  1E+01   2.182E-05  2.053E-05  2.665E-03
  1E+02   2.011E-04  1.347E-04  7.962E-04
  1E+03   2.116E-03  1.444E-03  4.597E-02
  1E+04   2.283E-02  1.587E-02  9.671E+00
"""


