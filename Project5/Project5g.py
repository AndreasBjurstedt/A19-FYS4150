
from time import perf_counter
import numpy as np
import matplotlib.pyplot as plt

#The program uses the Velocity Verlet method to simulate the
#perihelion percession of Mercury over TT years divided into n
#time steps. In order to avoid problems with insuficient compter
#memory, one year at the time is simulated. When one year is simulated,
#the progtam finds the positions of Mercury closest to the Sun that current year.
#The initial conditions for next year to be simulated are extracted from the last time
#step in the current year simulated.


#The Velocity Verlet method. 
def VelocityVerlet(f,x0,vx0,y0,vy0,t0,dt,n):

    
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


#The acceleration needed to find the velocity components in the 
#Velocity Verlet methodUnit AU/(year)**2 (AU - Astronomical Unit)
#The accelaration is without the general relativistic correction,
#that is acceleration derived from Newtons classical law of gravitation.
def f_N(p,v,k):
    f = np.array([-4*np.pi**2*p[k,0],-4*np.pi**2*p[k,1]])
    r = np.sqrt(p[k,0]**2 + p[k,1]**2)
    f = f/float(r**3)
    return f


#The acceleration needed to find
#the position and velocity components in the Velocity Verlet
#method. Unit AU/(year)**2 (AU - Astronomical Unit).
#The general relativistic correction is included.
#c is the speed of light with unit [AU/Year]
def f_E(p,v,k):
    c = 63241.077084266275                                         
    l = np.abs(p[k,0]*v[k,1]-p[k,1]*v[k,0])
    r = np.sqrt(p[k,0]**2 + p[k,1]**2)
    f = np.array([-4*np.pi**2*(1.0+3.0*l**2/(r**2*c**2))*p[k,0],-4*np.pi**2*(1.0+3.0*l**2/(r**2*c**2))*p[k,1]])
    f = f/float(r**3)
    return f

#Initial time t0
t0 = 0

#Initial position of Mercury: x0 = 0.3075 AU and y0 = 0 AU
x0_N = 0.3075
y0_N = 0
x0_E = 0.3075
y0_E = 0
#Initial velocity of Mercury: vx0 = 0 AU/year and vy0 = 12.44 AU/year
vx0_N = 0
vy0_N = 12.44 
vx0_E = 0
vy0_E = 12.44 
   
#Numer of years TT,total number of time steps n ans time step length
TT = 100
n = 2*(10**6)
dt = (TT - t0)/float(n*TT)

#Initializon of lists
perihel_N = []
#kvalue_N = []
perihel_E = []
#kvalue_E = []
tt_N = []
tt_E = []
x_N = []
y_N = []
x_E = []
y_E = []

#Initializon of year counter to be printed on screen
#during simulation
teller = 0

#Simulation of one year at the time, without relativistic correction 
#(index N for Newton) and with relativistic correction (index E for Einstein)
for i in range(TT):
    T = 1 + i

    [t,pvv_N,vvv_N,cpu_vv_N] = VelocityVerlet(f_N,x0_N,vx0_N,y0_N,vy0_N,t0,dt,n)

    [t,pvv_E,vvv_E,cpu_vv_E] = VelocityVerlet(f_E,x0_E,vx0_E,y0_E,vy0_E,t0,dt,n)
    
#Initial conditions for simulation of next year    
    t0 = t[n]
    x0_N = pvv_N[n,0]
    y0_N = pvv_N[n,1]
    x0_E = pvv_E[n,0]
    y0_E = pvv_E[n,1]
    vx0_N = vvv_N[n,0]
    vy0_N = vvv_N[n,1]
    vx0_E = vvv_E[n,0]
    vy0_E = vvv_E[n,1]
    
    
#Distances between Sun and Mercury (no relativistic correction)
    
    rr_N = np.sqrt(pvv_N[:,0]**2 + pvv_N[:,1]**2)
    #rmax_N = np.max(rr_N)
    #rmin_N = np.min(rr_N)
    
#Distances between Sun and Mercury (relativistic correction included)
    
    rr_E = np.sqrt(pvv_E[:,0]**2 + pvv_E[:,1]**2)
    #rmax_E = np.max(rr_E)
    #rmin_E = np.min(rr_E)

    
#Finding the positions and corresponding time steps where Mecury is
#closest to the Sun (no relativistic correction)   
    for k in range(1,np.size(rr_N)-1):

        if rr_N[k] < rr_N[k-1] and rr_N[k] < rr_N[k+1]:
            perihel_N.append(rr_N[k])
            tt_N.append(t[k])
            x_N.append(pvv_N[k,0])
            y_N.append(pvv_N[k,1])
            #kvalue_N.append(k)
            
    
#perihel_N = np.asarray(perihel_N)
#kvalue_N = np.asarray(kvalue_N) 

    
    
#--------------------------------   
    
#Finding the positions and corresponding time steps where Mecury is
#closest to the Sun (relativistic correction included)  
    for k in range(1,np.size(rr_E)-1):
        if rr_E[k] < rr_E[k-1] and rr_E[k] < rr_E[k+1]:
            perihel_E.append(rr_E[k])
            tt_E.append(t[k])
            x_E.append(pvv_E[k,0])
            y_E.append(pvv_E[k,1])  
            #kvalue_E.append(k)
            
#Printing curret year just simulated on screen.   
    teller = teller + 1
    print(teller)    

#--------------------------------------------

#Making arrays of lists with results from TT years of
#simulation
tt_N = np.asarray(tt_N)
x_N = np.asarray(x_N)
y_N = np.asarray(y_N)
perihel_N = np.asarray(perihel_N)
tt_E = np.asarray(tt_E)
x_E = np.asarray(x_E)
y_E = np.asarray(y_E) 
perihel_E = np.asarray(perihel_E)
#print(np.size(tt_N))








#x_N = np.zeros(np.size(kvalue_N))
#y_N = np.zeros(np.size(kvalue_N))
#theta_N = np.zeros(np.size(kvalue_N))

#for k in range(np.size(kvalue_N)):
#    x_N[k] = pvv_N[kvalue_N[k],0]
#    y_N[k] = pvv_N[kvalue_N[k],1]
#    theta_N[k] = np.arctan(y_N[k]/x_N[k])
    
#------------------------------------------

#x_E = np.zeros(np.size(kvalue_E))
#y_E = np.zeros(np.size(kvalue_E))
#theta_E = np.zeros(np.size(kvalue_E))

#for k in range(np.size(kvalue_E)):
#    x_E[k] = pvv_E[kvalue_E[k],0]
#    y_E[k] = pvv_E[kvalue_E[k],1]
#    theta_E[k] = np.arctan(y_E[k]/x_E[k])
    



#print(np.size(tt_N))
#print(np.size(tt_E))
#print(np.size(x_N))
#print(np.size(y_N))
#print(np.size(x_E))
#print(np.size(y_E))


#Writing final results to file
outfile = open('data5g_N.txt','w')
for i in range(np.size(x_N)):
     outfile.write("""%2.12f %2.12f %2.12f""" % (tt_N[i], x_N[i], y_N[i]))
     outfile.write('\n')
outfile.close()

outfile = open('data5g_E.txt','w')
for i in range(np.size(x_E)): 
     outfile.write("""%2.12f %2.12f %2.12f""" % (tt_E[i], x_E[i], y_E[i]))
     outfile.write('\n')
outfile.close()

    
#Plotting Mercurys orbit after TT rears of simulation
#(both with and without relativistic correction) as
#a check that they look okay.    
plott1 = 'ja'
#Plot of earth's position if plott = 'ja' 
if plott1 == 'ja':       


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(pvv_N[:,0], pvv_N[:,1],'b',label = 'Newton')
    ax.plot(pvv_E[:,0], pvv_E[:,1],'r',label = 'Einstein')
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







 

   






