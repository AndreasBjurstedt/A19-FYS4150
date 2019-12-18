

import numpy as np
import matplotlib.pyplot as plt
from Euler_and_Verlet_func import *


#The centripetal acceleration needed to find the velocity components in the 
#Forward Euler method and the position and velocity components in the Velocity Verlet
#method. Unit AU/(year)**2 (AU - Astronomical Unit)
def f(p,v,k):
    f = np.array([-4*np.pi**2*p[k,0],-4*np.pi**2*p[k,1]])
    r = np.sqrt(p[k,0]**2 + p[k,1]**2)
    f = f/float(r**3)
    return f

#Initial time t0
t0 = 0
#Initial position of planet Earth: x0 = 1 AU and y0 = 0 AU
x0 = 1
y0 = 0
#Initial velocity of planet Earth: vx0 = 0 AU/year and vy0 = 2*pi AU/year
vx0 = 0
vy0 = 2*np.pi    
#*np.sqrt(2-0.1)


##########Plot of Earth's orbit with Forward Eueler and Velocity Verlet method########

plott1 = 'ja'
#Plot of earth's position if plott = 'ja' 
if plott1 == 'ja':

    #Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
    T = 2
    #number of time steps in the analysis: dt = T/n
    n = 20000

    [t,peuler,veuler, cpu_fe] = ForwardEuler(f,x0,vx0,y0,vy0,t0,T,n)
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
    ax.plot(pvv[:,0], pvv[:,1], color = 'b',label = 'Velocity Verlet')
    ax.plot(peuler[:,0], peuler[:,1], color = 'g', label = 'Forward Euler')
    ax.plot(xs, ys, color = 'r', label = 'Circular orbit')
    #Increase margins on axes
    ax.set_xmargin(0.1)
    ax.axis('equal')
    #plt.axis('equal')
    ax.set_xlabel('x(t) [AU]', fontsize = 15)
    ax.set_ylabel('y(t) [AU]', fontsize = 15)
    ax.set_title('Planet Earth orbiting 2 times around the Sun', fontsize = 16) 
    ax.legend(loc='center', fontsize = 14)
    ax.tick_params(labelsize = 14)
    plt.show()


###Plot of maximum relative difference in kinetic energy, potential energy 
###and angular momentum of planet Earth during one orbit around the Sun (one year)

plott2 = 'nei'
#If plott2 = 'ja':
#Plot of maximum relative difference in kinetic energy, potential energy
#and angular momentum of planet Eart during one orvbit around the Sun.
#Also elapsed CPU time while running Forward Euler and while running
#Velocity Verlet. 
if plott2 == 'ja':

    #Time span of numerical analysis: t[0] = 0 to t[n] = T [year]
    T = 1
    #number of time steps n in the analysis, where: dt = T/n
    narray = np.array([10, 100, 1000, 10000, 10**5, 10**6])
    dE_kEuler = np.zeros(np.size(narray))
    dE_kvv =np.zeros(np.size(narray))
    dE_pEuler = np.zeros(np.size(narray))
    dE_ppvv =np.zeros(np.size(narray))
    dE_amEuler = np.zeros(np.size(narray))
    dE_ampvv =np.zeros(np.size(narray))
    cpu_FE = np.zeros(np.size(narray))
    cpu_VV = np.zeros(np.size(narray))
    cpu_factor = np.zeros(np.size(narray))

    for k in range (np.size(narray)):
        n = narray[k]
        [t,peuler,veuler,cpu_fe] = ForwardEuler(f,x0,vx0,y0,vy0,t0,T,n)
        [t,pvv,vvv,cpu_vv] = VelocityVerlet(f,x0,vx0,y0,vy0,t0,T,n)
 

        #Calculation of velocity squared, Forward Euler method [(AU/year)**2]
        veuler2 = veuler[:,0]**2 + veuler[:,1]**2
        veuler2_max = np.max(veuler2)
        veuler2_min = np.min(veuler2)

    
        #Calculation of velocity squared, Velocity Verlet method [(AU/year)**2]
        #Approximate velocity squared
        vvv2 = vvv[:,0]**2 + vvv[:,1]**2
        vvv2_max = np.max(vvv2)
        vvv2_min = np.min(vvv2)

        #Maximum relative difference in Earth's kinetic energy
        dE_kEuler[k] = (veuler2_max - veuler2_min)/np.mean(veuler2)
        dE_kvv[k] = (vvv2_max - vvv2_min)/np.mean(vvv2)

        #Calculation of position r, Forward Euler method [(AU/year)**2]
        ppeuler = np.sqrt(peuler[:,0]**2 + peuler[:,1]**2)
        ppeuler_max = np.max(ppeuler)
        ppeuler_min = np.min(ppeuler)

    
        #Calculation of position r, Velocity Verlet method [(AU/year)**2]
        ppvv = np.sqrt(pvv[:,0]**2 + pvv[:,1]**2)
        ppvv_max = np.max(ppvv)
        ppvv_min = np.min(ppvv)

        #Maximum relative difference in Earth's potential energy
        dE_pEuler[k] = (1.0/ppeuler_min - 1.0/ppeuler_max)/np.mean(1.0/ppeuler)
        dE_ppvv[k] = (1.0/ppvv_min - 1.0/ppvv_max)/np.mean(1.0/ppvv)

        #Calculation of angular momentum per unit mass, Forward Euler method
        ameuler = peuler[:,0]*veuler[:,1] - peuler[:,1]*veuler[:,0]
        ameuler_max = np.max(ameuler)
        ameuler_min = np.min(ameuler)

        #Calculation of angular momentum per unit mass, Velocity Verlet method
        amvv = pvv[:,0]*vvv[:,1] - pvv[:,1]*vvv[:,0]
        amvv_max = np.max(amvv)
        amvv_min = np.min(amvv)

        #Maximum relative difference in Earth's angular momentum
        dE_amEuler[k] = (ameuler_max - ameuler_min)/np.mean(ameuler)
        dE_ampvv[k] = (amvv_max - amvv_min)/np.mean(amvv)
        #print np.mean(amvv)
        #print amvv_max
        #print amvv_min

        #Elapsed CPU time while running Forvard Euler algorithm and when
        #running the the Velocity Verlet algorithm
        cpu_FE[k] = cpu_fe
        cpu_VV[k] = cpu_vv
        cpu_factor[k] = cpu_VV[k]/float(cpu_FE[k])

    #Table
    print('') 
    print("""Elapsed CPU time while running each of the two algorithms:""")
    print('')
    print('    n      Forward Euler  Velocity Verlet  VelVer/FEuler')               
    for i in range(0,np.size(narray)):
        print("""  %.0E     %.3E       %.3E        %.3f"""
        %(narray[i], cpu_FE[i],cpu_VV[i],cpu_factor[i]))




 

    #Plot of maximum relative difference in kinetic energy
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.loglog(narray, dE_kvv, 'bo-',label = 'Velocity Verlet')
    ax.plot(narray, dE_kEuler, 'go-', label = 'Forward Euler')
    #Increase margins on axes
    ax.set_xmargin(0.1)
    #ax.axis('equal')
    #plt.axis('equal')
    ax.set_xlabel('Number of time steps dt', fontsize = 15)
    ax.set_ylabel('Relative difference', fontsize = 15)
    ax.set_title('Maximum relative difference in kinetic energy', fontsize = 16) 
    #ax.set_title('Relative difference between \n maximum and minimum potential energy', fontsize = 16)
    ax.legend(loc='lower left', fontsize = 14)
    ax.tick_params(labelsize = 14)
    plt.grid(True)
    plt.show()

    #Plot of maximum relative difference in potential energy
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.loglog(narray, dE_ppvv, 'bo-',label = 'Velocity Verlet')
    ax.plot(narray, dE_pEuler, 'go-', label = 'Forward Euler')
    #Increase margins on axes
    ax.set_xmargin(0.1)
    #ax.axis('equal')
    #plt.axis('equal')
    ax.set_xlabel('Number of time steps dt', fontsize = 15)
    ax.set_ylabel('Relative difference', fontsize = 15)
    ax.set_title('Maximum relative difference in potential energy', fontsize = 16) 
    ax.legend(loc='lower left', fontsize = 14)
    ax.tick_params(labelsize = 14)
    plt.grid(True)
    plt.show()

    #Plot of maximum relative difference in angular momentum 
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.loglog(narray, dE_ampvv, 'bo-',label = 'Velocity Verlet')
    ax.plot(narray, dE_amEuler, 'go-', label = 'Forward Euler')
    #Increase margins on axes
    ax.set_xmargin(0.1)
    #ax.axis('equal')
    #plt.axis('equal')
    ax.set_xlabel('Number of time steps dt', fontsize = 15)
    ax.set_ylabel('Relative difference', fontsize = 15)
    ax.set_title('Maximum relative difference in angular momentum', fontsize = 16) 
    ax.legend(loc='center left', fontsize = 14)
    ax.tick_params(labelsize = 14)
    plt.grid(True)
    plt.show()


#Table from running the script

#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5c.py
#
#Elapsed CPU time while running each of the two algorithms:
#
#    n      Forward Euler  Velocity Verlet  VelVer/FEuler
#  1E+01     1.344E-04       2.197E-04        1.635
#  1E+02     9.788E-04       1.746E-03        1.784
#  1E+03     1.060E-02       1.786E-02        1.684
#  1E+04     9.866E-02       1.756E-01        1.779
#  1E+05     9.840E-01       1.762E+00        1.790
#  1E+06     9.906E+00       1.757E+01        1.773
#
#(base) M:\H2018\FYS4150-H19\Prosjekt5\Programmer>









