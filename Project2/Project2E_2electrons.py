
import numpy as np
#import scipy.linalg as scl
from Jacobi_Eigen import *


#Matrix dimension NxN 
N =160

#An array with varying strength of wr
wr_array = np.array([0.01,0.25,0.5,1,5])
#An array with varying strength of wr squared
wr2_array = np.array([0.01**2,0.25**2,0.5**2,1,5**2])

#Matrix for storing 4 lowest eigenvalues for all the wr values in wr_array
eig_A =np.zeros((np.size(wr_array),4))

j=0

rhoO = 0.0
rhoN = 10.0

for wr2 in wr2_array:
    #Distance between the discrete points
    h = (rhoN -rhoO)/(N+1)
    #Matrix element e on the diagonal just above and below the
    #main diagonal A[i,i] in the tridiagonal matrix A
    e = -1.0/(h**2)

    #The tridiagonal matrix A with element d on the main diagonal
    A = np.zeros((N,N))

    for i in range(N):
        d = 2.0/(h**2) + (wr2)*(rhoO + (i+1)*h)**2 + 1.0/(rhoO + (i+1)*h)
        A[i,i] = d

    for i in range(N-1):
        A[i+1,i] = e
        A[i,i+1] = e
    
    #print(h)
    #print(A)

    [A,st,cpu_j] = Eigen_Jacobi(A)
    a = np.diag(A)
    a = np.sort(a)
    eig_A[j,0] = a[0]
    eig_A[j,1] = a[1]
    eig_A[j,2] = a[2]
    eig_A[j,3] = a[3]
    j = j+1

print("""The four lowest scaled two-electron energies with varying parameter wr,""") 
print("""when using a %d x %d matrix in approximating them"""%(N,N))
print("""     wr         E1          E2          E3         E4""")
for i in range(0,np.size(wr_array)):
    print("""  %.1E   %.4E  %.4E  %.4E %.4E"""%(wr_array[i],eig_A[i,0],eig_A[i,1],eig_A[i,2],eig_A[i,3]))


#Running the script
"""
(base) M:\FYS4150-H19\Prosjekt2\Programmer>python Project2E_2electrons.py
The four lowest scaled two-electron energies with varying parameter wr,
when using a 160 x 160 matrix in approximating them
     wr         E1          E2          E3         E4
  1.0E-02   3.1163E-01  6.8179E-01  1.2228E+00 1.9470E+00
  2.5E-01   1.2499E+00  2.1898E+00  3.1498E+00 4.1231E+00
  5.0E-01   2.2298E+00  4.1331E+00  6.0705E+00 8.0244E+00
  1.0E+00   4.0566E+00  7.9039E+00  1.1805E+01 1.5730E+01
  5.0E+00   1.7417E+01  3.6923E+01  5.6489E+01 7.6027E+01
 (base) M:\FYS4150-H19\Prosjekt2\Programmer> 
"""

