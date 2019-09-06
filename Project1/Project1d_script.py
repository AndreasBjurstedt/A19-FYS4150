
#Using the special tridiagonal algorithm implemented in the 
#function TriDiagSpes, which functionality is explained
#in Diagonal_Solvers.py, in order to compute the relative error for 
#increasing number of grid points n

import numpy as np
import matplotlib.pyplot as plt
from Diagonal_Solvers import *

#An array with increasing number of grid points n 
n_array = np.array([10, 100, 1000, 10000, 10**5, 10**6, 10**7])
rel_error =np.zeros(np.size(n_array))

#Finding the relative error for each n in n_array
i=0
for n in n_array:
    #step size
    h = 1.0/(n+1)

    x = np.linspace(h, n*h, n)
    f = np.zeros(n)
	
    #Our second derivative u'' = f:
    f = h**2*(100*np.exp(-10*x))
	
    #Our exact solution in the interval [0+h, 1-h] of x:	
    #xx = np.linspace(0+h,1-h,n)
    u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

    #Calling the function with the special tridiagonal algorithm
    [v,flop,cpu_s] = TriDiagSpes(n,f)
	
	
    relerr = np.abs((v - u)/u)
    relerr = np.max(relerr)	
    relerr = np.log10(relerr) 
    rel_error[i] = relerr
    i = i+1

#Table of relative errors
print('Table of e(n) = log10(abs(relative error)):')
print('    n        e(n)')               
for i in range(0,np.size(rel_error)):
    print("""  %.0E   %.4E"""%(n_array[i],rel_error[i]))
	

#Plot of the relative error.

plt.figure()
ax = plt.subplot(111) 
ax.set_xscale("log") 
plt.plot(n_array,rel_error,'bo-')
#plt.axis([0, max(xx), min(u) - 0.5, max(u) + 0.5])   # [ xx_min, xx_max, u_cordmin, u_max]
plt.xlabel('Grid points n', fontsize = 14)
plt.ylabel('log10 of maximum relative error', fontsize = 14)
plt.title("""Maximum relative error when using 10**n grid points""", fontsize = 14) 
ax.tick_params(labelsize = 14)
plt.show()


#Running the script:
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python Project1d_script.py
Table of e(n) = log10(abs(relative error)):
    n        e(n)
  1E+01   -1.1797E+00
  1E+02   -3.0880E+00
  1E+03   -5.0801E+00
  1E+04   -7.0793E+00
  1E+05   -9.0776E+00
  1E+06   -1.0145E+01
  1E+07   -9.0902E+00
"""

















