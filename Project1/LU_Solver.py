

import numpy as np
import scipy.linalg as scl
from time import perf_counter


def LU_Solver(A,f,n):
    """
    Solves a system of n linear equations with n unknowns with
    the aid of LU decomposition of the matrix A : A*v = f => (L*U)*v = f
    => U*v = (inverse(L))*f. Then y = inverse(L)*f and we find v with
    backward substitution of y = U*v   
    """
    y = np.zeros(n)
    v = np.zeros(n)

    #Starting timer
    start = perf_counter()

    #LU decomposition of quadratic matrix A
    P,L,U = scl.lu(A)

    #Inverting L
    #LL = scl.inv(L)
    #Matrix multiplication between matrix LL and vector f
    #y = np.dot(LL,f)
	
    #Forward substitution to find vector y = L**(-1)*f
    for i in range(n):
        y[i] = f[i] - np.dot(L[i,:],y) 	

    #Backward substitution to find solution vector v
    for i in range(n-1,-1,-1):
        v[i] = (float(y[i] - np.dot(U[i,:],v)))/U[i,i]

    #Elapsed cpu time during the run of the algoritm.
    slutt = perf_counter()
    total = slutt - start
    return v, total


#Testing of the LU solver
"""
a = -1
b = 2
c = -1
f = np.asarray([1,2,3,4])
f = 1.0*f

n = 4
A = np.zeros((n,n))
A[0,0] = b
A[0,1] = c
A[n-1,n-1] = b
A[n-1,n-2] = a

for i in range(1,n-1):
    A[i,i] = b
    A[i, i -1] = a
    A[i,i+1] = c

print('Testing of the LU solver:')
[v,Total_LU] = LU_Solver(A,f,n)
print('v = %s'%v)
print('Elapsed CPU time: %G s'%Total_LU)
"""

#The result of the testing
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python LU_Solver.py
Testing of the LU solver:
v = [4. 7. 8. 6.]
Elapsed CPU time: 0.00305115 s
"""





