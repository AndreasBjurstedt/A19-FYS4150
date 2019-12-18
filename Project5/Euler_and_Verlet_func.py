
import numpy as np
#import time
from time import perf_counter

def ForwardEuler(f,x0,vx0,y0,vy0,t0,T,n):
    """ 
    Uses the Forward Euler method to find the path of a planet around the Sun 
    in the solar system. The function f gives the planet centripetal acceleration in x 
    and y direction. Initial planet position components are given by x0 and y0. 
    Its initial speed components are given by vx0 and vy0. The calculations are
    preformed from t[0] = t0 to t[n] = T with time step dt = T/n.
    """

    
    #Time vector t[0],t[1],t[2],......,t[n]
    t = np.zeros(n+1)
    
    #Vector for p for x and y component of position
    #Vector v for x and y component of speed
    p = np.zeros([n+1,2])
    v = np.zeros([n+1,2])
    
    #Initial values of position and speed in x and y direction 
    #when t = t[0], and for time
    p[0,:] = np.array([float(x0),float(y0)])
    v[0,:] = np.array([float(vx0),float(vy0)])
    t[0] = float(t0)
    
    #Time step between time-grid points
    dt = (T - t[0])/float(n)

    
    #Starting timer
    #The Forward Euler method
    start = perf_counter()
    #c0 = time.time()
    for k in range(n):
        t[k+1] = t[k] + dt
        p[k+1] = p[k] + dt*v[k]
        v[k+1] = v[k] + dt*f(p,v,k)
    #cpu time while it is running the Forward Euler algorithm.
    slutt = perf_counter()
    cpu_fe = slutt - start
    #cpu_fe = time.time() - c0
    return t,p,v, cpu_fe
    

def VelocityVerlet(f,x0,vx0,y0,vy0,t0,T,n):

    
    #Time vector t[0],t[1],t[2],......,t[n]
    t = np.zeros(n+1)
    
    #Vector for p for x and y component of position
    #Vector v for x and y component of speed
    p = np.zeros([n+1,2])
    v = np.zeros([n+1,2])
    
    #Initial values for position and speed when t = t[0],
    #and for time
    p[0,:] = np.array([x0,y0])
    v[0,:] = np.array([vx0,vy0])
    t[0] = t0
    
    #Time step between time-grid points
    dt = (T - t[0])/float(n)
    
    #Starting timer
    #The Velocity Verlet method
    #c1 = time.time()
    start = perf_counter()
    for k in range(n):
        fpk = f(p,v,k)
        t[k+1] = t[k] + dt
        p[k+1] = p[k] + dt*v[k] + 0.5*dt**2*fpk
        v[k+1] = v[k] + 0.5*dt*(f(p,v,k+1) + fpk)       
    #cpu time while it is running the Velocity Verlet algorithm.
    slutt = perf_counter()
    cpu_vv = slutt - start
    #cpu_vv = time.time() - c1
    return t,p,v, cpu_vv












