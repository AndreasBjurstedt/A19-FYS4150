

import numpy as np
#import time
from time import perf_counter

def TriDiag(n,a,b,c,f):
    """
    Solves a system of n linear equations with n unknowns,
    when the matrix A in the nxn matrix system A*v = f is 
    tri-diagonal: 
    The array b of length n is the numbers on the diagonal of A. 
    If all the numbers are the same, b can be given as that single number.
    The array c of length n-1 is the numbers just above the diagonal of A. 
    If all the numbers are the same, c can be given as that single number.
    The array a of length n-1 is the numbers just below the diagonal of A. 
    If all the numbers are the same, a can be given as that single number.
    The rest of the numbers in A are equal to zero.
    """
    flop = 0
    #Checks out if b is a number. Assumes
    #then that the user has given a and c respectively as numbers too.
    if isinstance(b,(float,int)):
	    #Including an extra number in the array 'b' in order to make the array
        #compatible with the algorithm below.
        a = a*np.ones(n-1)
        b = b*np.ones(n)
        c = c*np.ones(n-1)
    else:
        #List 'a' converted to array 'a'. Also make
        #sure that array 'a' is a 'float' array, not an 'int' array.
        a = np.asarray(a)
        a = 1.0*a
        #a = np.hstack((0.0,a))
        #List 'b' converted to array 'b'. Also making sure that array 'b' is a
        #'float' array by multiplying it with 1.0.
        b = np.asarray(b)
        b = 1.0*b
  
        c = np.asarray(c)
        c = 1.0*c

    f = np.asarray(f)
    f = 1.0*f	
    v = np.zeros(n)

    #Starting timer
    start = perf_counter()
    #Forward substitution step:
    for i in range(1,n):
        d = (float(a[i-1])/b[i-1])
        b[i] = b[i] - d*c[i-1]
        f[i] = f[i] - d*f[i-1]
        flop = flop + 5

    #Backward substitution step
    v[n-1] =(f[n-1])/b[n-1]
    flop = flop + 1
    for i in range(n-1,0,-1):
        v[i-1] =(f[i-1] - c[i-1]*v[i])/b[i-1]
        flop = flop + 3
    #Elapsed cpu time while it is running the algorithm.
    slutt = perf_counter()
    total = slutt - start
    return v, flop, total




def TriDiagSpes(n,f):
    """
    Solves a system of n linear equations with n unknowns,
    when the matrix A in the nxn matrix system A*v = f is 
    tri-diagonal: 
    All the numbers of the diagonal of A is equal to 2.
    The numbers just above and below the diagonal is equal to -1.
    The rest of the numbers in A are equal to zero.
    """
    flop = 0
    a = -1*np.ones(n-1)
    b = 2*np.ones(n)
    c = -1*np.ones(n-1)
   
    #List 'f' converted to array 'f'. Also making sure that array 'f' is a
    #'float' array by multiplying it with 1.0.
    f = np.asarray(f)
    f = 1.0*f	
	
    v = np.zeros(n)

    #Starting timer
    start = perf_counter()
    #c0 = time.time()
    #Forward substitution step:
    for i in range(1,n):
        b[i] = float((i+2))/(i+1)
        f[i] = f[i] + ((float(i))/(i+1))*f[i-1]
        flop = flop + 2

    #Backward substitution step
    v[n-1] = (f[n-1])/b[n-1]
    flop = flop + 1
    for i in range(n-1,0,-1):
        v[i-1] = (float(i)/(i+1))*(f[i-1] + v[i])
        flop = flop + 2
    #Elapsed cpu time while it is running the algorithm.
    slutt = perf_counter()
    total = slutt - start
    return v,flop,total




#Testing of the two tridiagonal solvers
"""
a = -1
b = 2
c = -1
n = 4
f = [1,2,3,4]

print('--------------------------')
print('Testing of the general tridiagonal solver:')
[v1,flop_g, total_g] = TriDiag(n,a,b,c,f)
print('v = %s'%v1)
print('Floating Point Operations for n = %d: %d'%(n,flop_g))
print('Elapsed CPU time: %G s'%total_g)


print('--------------------------')
print('Testing of the special tridiagonal solver:')
[v2,flop_s, total_s] = TriDiagSpes(n,f)
print('v = %s'%v2)
print('Floating Point Operations for n = %d: %d'%(n,flop_s))
print('Elapsed CPU time: %G s'%total_s)
"""

#The result of the testing
"""
(base) M:\FYS4150-H19\Prosjekt1\Programmer>python Diagonal_Solvers.py
--------------------------
Testing of the general tridiagonal solver:
v = [4. 7. 8. 6.]
Floating Point Operations for n = 4: 25
Elapsed CPU time: 1.7965E-05 s
--------------------------
Testing of the special tridiagonal solver:
v = [4. 7. 8. 6.]
Floating Point Operations for n = 4: 13
Elapsed CPU time: 1.7002E-05 s
"""

