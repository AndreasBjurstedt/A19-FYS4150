

import numpy as np
import matplotlib.pyplot as plt
from Euler_and_Verlet_func import *


#The centripetal acceleration needed to find the velocity components in the 
#Forward Euler method and the position and velocity components in the Velocity Verlet
#method. Unit AU/(year)**2 (AU - Astronomical Unit)

#Factor bb in r**bb in f(p,k) in interval [3,4].
#Physically correct factor: bb = 2

def f(p,k):    
    f = np.array([-4*np.pi**2*p[k,0],-4*np.pi**2*p[k,1]])
    r = np.sqrt(p[k,0]**2 + p[k,1]**2)
    f = f/float(r**(bb+1))
    return f

#Initial time t0
t0 = 0
#Initial position of planet Earth: x0 = 1 AU and y0 = 0 AU
x0 = 1
y0 = 0
#Initial velocity of planet Earth: vx0 = 0 AU/year and vy0 = 2*pi AU/year
vx0 = 0
#vy0 = 2*np.pi

bb = 2.0    

##-----Test 1, escape velocity
#vy0 = 2*np.pi*np.sqrt(2-0.5)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
#T = 4
##-----number of time steps in the analysis: dt = T/n
#n = 20000
#tittel = """Initial velocity of Planet earth:\n v_x0 = 0 and v_y0 = 2*pi*(2-0.5)**0.5"""
#posisjon = """lower left"""

##-----Test 2, escape velocity
#vy0 = 2*np.pi*np.sqrt(2-0.1)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
#T = 40
##-----number of time steps in the analysis: dt = T/n
#n = 50000
#tittel = """Initial velocity of Planet earth:\n v_x0 = 0 and v_y0 = 2*pi*(2-0.1)**0.5"""
#posisjon = """lower left"""

##-----Test 3, escape velocity
#vy0 = 2*np.pi*np.sqrt(2-0.01)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
#T = 1500
##-----number of time steps in the analysis: dt = T/n
#n = 100000
#tittel = """Initial velocity of Planet earth:\n v_x0 = 0 and v_y0 = 2*pi*(2-00.1)**0.5"""
#posisjon = """lower left"""


##-----Test 4: Analytical escape velocity
#vy0 = 2*np.pi*np.sqrt(2)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
#T = 1500
##-----number of time steps in the analysis: dt = T/n
#n = 100000
#tittel = """Initial velocity of Planet earth:\n v_x0 = 0 and v_y0 = 2*pi*(2)**0.5"""
#posisjon = """lower left"""


##-----Test5,  Factor bb = 2.99 and slight increase in y component of initial velocity.
#bb = 2.99 
#vy0 = 2*np.pi*np.sqrt(1.001)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
#T = 100
##-----number of time steps in the analysis: dt = T/n
#n = 100000
#tittel = """Orbit of Planet earth with beta = 2.99, v_x0= 0 and\n v_y0 = 2*pi*(1.001)**0.5"""
#posisjon = """center"""


##-----Test6,  Factor bb = 3 and slight increase in y component of initial velocity. 
bb = 3.0
vy0 = 2*np.pi*np.sqrt(1.001)
##-----Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
T = 100
##-----number of time steps in the analysis: dt = T/n
n = 100000
tittel = """Orbit of Planet earth with beta = 3, v_x0= 0 and\n v_y0 = 2*pi*(1.001)**0.5"""
posisjon = """lower right"""


plott1 = 'ja'
#Plot of earth's position if plott = 'ja' 
if plott1 == 'ja':
    
    [t,pvv,vvv,cpu_vv] = VelocityVerlet(f,x0,vx0,y0,vy0,t0,T,n)


    #A perfect circle 2*np.pi
    ns = 1000
    Ts = 2*np.pi
    ts = np.zeros(ns+1)
    xs = np.zeros(ns+1)
    ys = np.zeros(ns+1)
    dts = (Ts)/float(ns)

    xs[0] = 1.0
    ys[0] = 0.0
    for i in range(ns):
            ts[i+1] = ts[i] + dts
            xs[i+1] = np.cos(ts[i+1])
            ys[i+1] = np.sin(ts[i+1])
        

    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(pvv[:,0], pvv[:,1], color = 'b',label = 'Earth orbit')
    ax.plot(xs, ys, color = 'r', label = 'Circular orbit')
    #Increase margins on axes
    ax.set_xmargin(0.1)
    ax.axis('equal')
    #plt.axis('equal')
    ax.set_xlabel('x(t) [AU]', fontsize = 15)
    ax.set_ylabel('y(t) [AU]', fontsize = 15)
    ax.set_title(tittel, fontsize = 16) 
    ax.legend(loc= posisjon, fontsize = 14)
    ax.tick_params(labelsize = 14)
    plt.show()

