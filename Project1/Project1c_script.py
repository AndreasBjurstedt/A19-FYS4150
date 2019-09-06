
#Comparing the cpu time for the special tridiagonal algorithm
#and the general tridiagonal algorithm, which functionality is explained
#in Diagonal_Solvers.py, for increasing number of grid points n
#Also giving the number of flops for both algorithms-

import numpy as np
from Diagonal_Solvers import *

#An array with increasing number of grid points n 
n_array = np.array([10, 100, 1000, 10**4, 10**5, 10**6])
flop_special = np.zeros(np.size(n_array))
flop_general = np.zeros(np.size(n_array))
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

    #Calling the function with the general tridiagonal algoritm
    a = -1
    b = 2
    c = -1
    [v,flop_g,cpu_g] = TriDiag(n,a,b,c,f)

    
    flop_special[i] = flop_s
    flop_general[i] = flop_g
    cpu_special[i] = cpu_s
    cpu_general[i] = cpu_g

    i = i+1

#Table
print('') 
print('Table with  elapsed cpu time and flops for both algorithms:')
print('')
print('            elapsed cpu time                flops  ')
print('    n      general    special       general       special')               
for i in range(0,np.size(n_array)):
    print("""  %.0E   %.3E  %.3E  %10.0f    %10.0f"""
    %(n_array[i], cpu_general[i],cpu_special[i],flop_general[i], flop_special[i]))

#Running the script
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python Project1c_script.py

Table with  elapsed cpu time and flops for both algorithms:

            elapsed cpu time                flops
    n      general    special       general       special
  1E+01   2.181E-05  2.021E-05          73            37
  1E+02   1.950E-04  1.293E-04         793           397
  1E+03   2.104E-03  1.382E-03        7993          3997
  1E+04   2.068E-02  1.522E-02       79993         39997
  1E+05   2.108E-01  1.627E-01      799993        399997
  1E+06   2.097E+00  1.459E+00     7999993       3999997
"""





