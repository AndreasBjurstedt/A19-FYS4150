

import numpy as np
#import scipy.linalg as scl
from time import perf_counter



def Eigen_Jacobi(A):
    """
    The function Eigen_Jacobi(A) is diagonalizing the symmetric nxn matrix A using 
    Jacobi's method. The matrix eigenvalues then appear on the diagonalized matrix
    diagonal. The function returns the diagonal matrix with the eigenvalues, the number
    of times the matrix A had to be transformed by a rotation transformation in order
    to be diagonalized, and the elapsed CPU time used from the start to the end of
    diagonalizing process. 
    """
    
    a = np.shape(A)
    n = a[0]
    
    #Initial value of the norm of the non-diagonal elements in the symmetric
    #matrix A n order to get in to while loop.	
    snond = 1.0
    #Number of rotation transformations
    st = 0

    #Starting timer
    start = perf_counter()

    while snond > 1.0*10**(-8):  # and st < n**3:
        snond = 0.0
        Amax = 0.0
        #Number of rotation transformations increased by one more
        st = st+1
 
    #Finding  the abs(maximum value) non-diagonal element in the symmetric matrix A
    #and the norm of the non-diagonal elements in A.	
        
        for i in range(n-1):
            for j in range(i+1,n):
                snond = snond + 2*A[i,j]**2
                if abs(A[i,j]) >= Amax:
                    Amax = abs(A[i,j])
                    k = i
                    l = j
        snond = np.sqrt(snond)      
    #The rotation transfomation preserving the eigenvalues of matrix A
    #First task: Calculating c = cos(theta) and s = sin(theta) for the orthogonal
    #rotation matrix 
    
        if A[k,l] != 0:
            tau = (1.0*A[l,l] - A[k,k])/(2*A[k,l])
           
            if tau > 0:
                t = -tau + np.sqrt(tau**2 + 1)
            else:
                t = -tau - np.sqrt(tau**2 + 1)

            c = 1.0/(np.sqrt(1 + t**2))
            s = t*c
        else:
            c = 1
            s = 0

    #The rotation transfomation preserving the eigenvalues of matrix A
    #Second task: Changing matrix elements with kk, ll, kl and lk indices

        a_kk = A[k,k]
        a_ll = A[l,l]
        A[k,k] = a_kk*c*c - 2*A[k,l]*c*s + a_ll*s*s
        A[l,l] = a_ll*c*c + 2*A[k,l]*c*s + a_kk*s*s
        A[k,l] = 0.0
        A[l,k] = 0.0
        

    #The similarity transfomations preserving the eigenvalues of matrix A
    #Third task: Changing matrix elements with ik, ki, il and li indices,
    #where i!=k and i!=l
    
        for i in range(n):
            if i!=k and i!=l:
                a_ik = A[i,k]
                a_il = A[i,l]
                A[i,k] = a_ik*c - a_il*s
                A[i,l] = a_il*c + a_ik*s
                A[k,i] = A[i,k]
                A[l,i] = A[i,l]
            
         
    
    
    #Elapsed cpu time for finding the diagonal matrix.
    slutt = perf_counter()
    total = slutt - start
    #cpu_j = time.time() - c0

    


    return A, st, total


#For testing of the function: see Project 2C_Unit_tests.py


 

