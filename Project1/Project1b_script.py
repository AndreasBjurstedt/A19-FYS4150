
#Using the general tridiagonal algorithm implemented in the 
#function TriDiag(n,a,b,c,f), which functionality is explained
#in Diagonal_Solvers.py, in order to create plots and to
#find the algorithm's number of Floating Point operations for a specific
#number of grid points n  

import numpy as np
import matplotlib.pyplot as plt
from Diagonal_Solvers import *

#Number of grid points 
n = 1000

#step size
h = 1.0/(n+1)

x = np.linspace(h, n*h, n)
f = np.zeros(n)

#Our second derivative u'' = f:
f = h**2*(100*np.exp(-10*x))

#Calling the function with the general tridiagonal algorithm
a = -1
b = 2
c = -1
[v,flop_g,total_g] = TriDiag(n,a,b,c,f)
v = np.hstack((0.0,v,0.0))
#-----------------------------------------------

#Our exact solution:	
xx = np.linspace(0,1,1000)
u = np.zeros(1000)
#for elements in xx:
u = 1 - (1 - np.exp(-10))*xx - np.exp(-10*xx)
#-----------------------------------------------

print('Floating Point Operations for n = %d: %d'%(n,flop_g))
print('Elapsed CPU time for n = %d: %G s'%(n,total_g))

#For large n, the cpu requires a lot of time to create the plot.
if n < 10000:
    x = np.hstack((0.0,x,1.0))
    plt.figure()
    ax = plt.subplot(111)
    #plt.subplot(3,1,1)  
    plt.plot(xx, u, 'r')
    plt.plot(x,v,'b')
    #plt.axis([0, max(xx), min(u) - 0.5, max(u) + 0.5])   # [ xx_min, xx_max, u_cordmin, u_max]
    plt.xlabel('x',fontsize = 14)
    plt.ylabel('u(x)',fontsize = 14)
    plt.title("""Numerical solution (blue) and exact analytical \n solution (red)""", fontsize = 14) 
    ax.tick_params(labelsize = 14)
    plt.show()


#Running the script for n = 10
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python Project1b_script.py
Floating Point Operations for n = 10: 73
Elapsed CPU time for n = 10: 2.3739E-05 s
"""


