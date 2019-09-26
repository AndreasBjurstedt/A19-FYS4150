
import numpy as np
import scipy.linalg as scl
from Jacobi_Eigen import *
from time import perf_counter

#------------We make a table with eigenvalues:--------------------------------

#number of discrete points, excluding the two end points
N = 6
#Distance between the discrete points
h = 1.0/(N+1)
#Matrix elements in the tridiagonal matrix
d = 2.0/(h**2)
a = -1.0/(h**2)

#The tridiagonal matrix
A = np.zeros((N,N))

for i in range(N):
    A[i,i] = d

for i in range(N-1):
    A[i+1,i] = a
    A[i,i+1] = a


#Analytical expression for the eigenvalues of the tridiagonal matrix
eig_a = np.zeros(N)
for j in range(N):
    eig_a[j] = d + 2*a*np.cos((j+1)*np.pi/(N+1))
    

#Scipy function for calculating the eigenvalues of the tridiagonal matrix
eig_p = scl.eigvals(A)
 

#Finding the eigenvalues of the tridiagonal matrix by using Jacobi's method
[A,st,cpu_j] = Eigen_Jacobi(A)


#Table with eigenvalues

    
#Sorting the numerical eigenvalues from lowest value up to highest value
eig_p = np.sort(eig_p)
a = np.diag(A)
a = np.sort(a)

print("""Table with eigenvalues for a %d x %d tridiagonal matrix."""%(N,N))
print("""  Analytical       Python          Jacobis method""")               
for i in range(0,N):
    print("""  %.8E   %.8E  %.8E"""%(eig_a[i],eig_p[i].real,a[i]))

#-------------------------------------------------------------------

#---------Number of rotation transformations and elapsed cpu time--------

#An array with increasing matrix dimension NxN 
N_array = np.array([10, 20, 40, 80, 160])

#Arrays for elapsed cpu times and number of rotation transformations
#and rotation transformations vs N
cpu_pa =np.zeros(np.size(N_array))
cpu_ja =np.zeros(np.size(N_array))
st_array =np.zeros(np.size(N_array))


j=0

for N in N_array:
    #Distance between the discrete points
    h = 1.0/(N+1)
    #Matrix elements in the tridiagonal matrix
    d = 2.0/(h**2)
    a = -1.0/(h**2)

    #The tridiagonal matrix
    A = np.zeros((N,N))

    for i in range(N):
        A[i,i] = d

    for i in range(N-1):
        A[i+1,i] = a
        A[i,i+1] = a


    #Elapsed CPU time for scipy function calculating the eigenvalues.
    start = perf_counter()
    eig_p = scl.eigvals(A)
    slutt = perf_counter()
    total = slutt - start
    cpu_pa[j] = total

    #Elapsed CPU time and number of rotation transformations for Jacobis method finding
    #the eigenvalues.	
    [A,st,cpu_j] = Eigen_Jacobi(A)	
    cpu_ja[j] = cpu_j
    st_array[j] = st
    j = j+1
	
    
    
#Rotation transformation factor as N increases
rotfac = np.zeros(np.size(N_array))
for i in range(1,np.size(N_array)):
    rotfac[i] = float(st_array[i])/st_array[i-1]
	 
   	
print(""" """)
print("""Table with elapsed cpu time for Scipy eigenvalue function and""")
print("""Jacobis method for finding eigenvalues. Also includes number of rotation""")
print("""transformations  rotrans used in Jacobis method for finding the eigenvalues,""")
print("""and the factor between number of rotation transformations as N doubles""")
print("""   N        cpu_python  cpu_Jacobi  rotrans     rotrans(N)/rotrans(N/2)""")               
for i in range(0,np.size(N_array)):
    print("""  %.1E   %.4E  %.4E  %.4E  %.3E"""%(N_array[i],cpu_pa[i],cpu_ja[i],st_array[i],rotfac[i]))


#----------Running the script-----------------------------------------
"""
(base) M:\FYS4150-H19\Prosjekt2\Programmer>python Project2b.py
Table with eigenvalues for a 6 x 6 tridiagonal matrix.
  Analytical       Python          Jacobis method
  9.70505095E+00   9.70505095E+00  9.70505095E+00
  3.68979994E+01   3.68979994E+01  3.68979994E+01
  7.61929485E+01   7.61929485E+01  7.61929485E+01
  1.19807052E+02   1.19807052E+02  1.19807052E+02
  1.59102001E+02   1.59102001E+02  1.59102001E+02
  1.86294949E+02   1.86294949E+02  1.86294949E+02
  
Table with elapsed cpu time for Scipy eigenvalue function and
Jacobis method for finding eigenvalues. Also includes number of rotation
transformations  rotrans used in Jacobis method for finding the eigenvalues,
and the factor between number of rotation transformations as N doubles
   N        cpu_python  cpu_Jacobi  rotrans     rotrans(N)/rotrans(N/2)
  1.0E+01   3.1824E-04  1.0363E-02  1.5800E+02  0.000E+00
  2.0E+01   1.9537E-04  1.4329E-01  6.7900E+02  4.297E+00
  4.0E+01   4.1287E-04  2.1992E+00  2.8400E+03  4.183E+00
  8.0E+01   4.4992E-03  3.4348E+01  1.1589E+04  4.081E+00
  1.6E+02   1.0029E-02  5.4968E+02  4.7307E+04  4.082E+00

(base) M:\FYS4150-H19\Prosjekt2\Programmer>
"""
















