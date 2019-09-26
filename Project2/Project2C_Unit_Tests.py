
import numpy as np
import scipy.linalg as scl
from Jacobi_Eigen import *

#Test 1: An algortihm searching for the largest non-diagonal element
#in the symmetrical matrix A is neccesary in order to find the eigenvalues
#of this matrix with the Jacobi-method. This algorithm is part of the
#function "Jacobi_Eigen.py", and it is tested here.  

#A symmetric test matrix A:
A = np.matrix([[1,2,3,4],[2,7,9,1],[3,9,6,3], [4,1,3,8]])
A=1.0*A
#Finding  the abs(maximum value) element in the symmetric matrix A

a = np.shape(A)
n = a[0]
Amax = 0.0
        
for i in range(n-1):
    for j in range(i+1,n):
        if abs(A[i,j]) >= Amax:
            Amax = abs(A[i,j])
            k = i
            l = j

print("""Largest matrix value: A[%d,%d] = A[%d,%d] = %d."""%(k,l,l,k,Amax))

#Test2: Checking that the function "Eigen_Jacobi(A)" finds the correct eigenvalues
#of the matrix A. The result is compared with the eigenvalues which the eigvals function
#from the Scipy libary calculates:


#Scipy function for calculating the eigenvalues of the symmetric matrix
eig_p = scl.eigvals(A)
 

#Finding the eigenvalues of the tridiagonal matrix by using Jacobi's method
[A,st,total] = Eigen_Jacobi(A)

#Sorting the numerical eigenvalues from lowest value up to highest value
eig_p = np.sort(eig_p)
a = np.diag(A)
a = np.sort(a)

print(""" """)
print("""Eigenvalues for the %d x %d symmetric matrix A."""%(n,n))
print("""   Python function   Jacobis method""")               
for i in range(0,n):
    print("""  %15.8E   %15.8E """%(eig_p[i].real,a[i]))


#Running the script
"""
(base) M:\FYS4150-H19\Prosjekt2\Programmer>python Project2C_Unit_Tests.py
[abjurste@fleshkind Programmer]$ python Project2C_Unit_Tests.py 
Largest matrix value: A[1,2] = A[2,1] = 9.
 
Eigenvalues for the 4 x 4 symmetric matrix A.
   Python function   Jacobis method
  -2.77703310E+00   -2.77703310E+00 
  -1.06939477E+00   -1.06939477E+00 
   8.06906286E+00    8.06906286E+00 
   1.77773650E+01    1.77773650E+01 
(base) M:\FYS4150-H19\Prosjekt2\Programmer>
"""



