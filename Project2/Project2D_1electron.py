
import numpy as np
#import scipy.linalg as scl
from Jacobi_Eigen import *


#An array with increasing matrix dimension NxN
#N_array = np.array([5,10]) 
N_array = np.array([10,20,40,80,160])

#Matrix for storing 4 lowest eigenvalues for all the N values in N_array
eig_A =np.zeros((np.size(N_array),4))

j=0

rhoO = 0
rhoN = 5.0
for N in N_array:
    #Distance between the discrete points
    h = (rhoN -rhoO)/(N+1)
    #Matrix element e on the diagonal just above and below the
    #main diagonal A[i,i] in the tridiagonal matrix A
    e = -1.0/(h**2)

    #The tridiagonal matrix A with element d on the main diagonal
    A = np.zeros((N,N))

    for i in range(N):
        d = 2.0/(h**2) + (rhoO + (i+1)*h)**2
        A[i,i] = d

    for i in range(N-1):
        A[i+1,i] = e
        A[i,i+1] = e
    
    #Finding the eigenvalues of A
    [A,st,cpu_j] = Eigen_Jacobi(A)
    a = np.diag(A)
    a = np.sort(a)
    #print(a)
    eig_A[j,0] = a[0]
    eig_A[j,1] = a[1]
    eig_A[j,2] = a[2]
    eig_A[j,3] = a[3]
    j = j+1


print("""The four lowest scaled one-electron energies when using an NxN matrix""")
print("""in approximating them""")
print("""     N         E1          E2          E3         E4""")
for i in range(0,np.size(N_array)):
    print("""  %.1E   %.4E  %.4E  %.4E %.4E"""%(N_array[i],eig_A[i,0],eig_A[i,1],eig_A[i,2],eig_A[i,3]))

#Running the script
"""
(base) M:\FYS4150-H19\Prosjekt2\Programmer>python Project2D_1electron.py
The four lowest scaled one-electron energies when using an NxN matrix
in approximating them
     N         E1          E2          E3         E4
  1.0E+01   2.9338E+00  6.6598E+00  1.0142E+01 1.3345E+01
  2.0E+01   2.9822E+00  6.9102E+00  1.0780E+01 1.4592E+01
  4.0E+01   2.9953E+00  6.9767E+00  1.0943E+01 1.4900E+01
  8.0E+01   2.9988E+00  6.9940E+00  1.0986E+01 1.4979E+01
  1.6E+02   2.9997E+00  6.9985E+00  1.0997E+01 1.4999E+01
(base) M:\H2018\FYS4150-H19\Prosjekt2\Programmer>
"""    






