
##Using the classes in Project5_classes.py to calculate the orbits of planet Earth
#and Jupiter when the Sun is fixed in Origo.
#Calculating maximum relative difference in total energy
#and maximum relative difference in angular momentum for the system

import numpy as np
import matplotlib.pyplot as plt
from Project5_classes import *

#Mass factor for planet Jupiter
mf = 1000

#Time span of numerical analysis: t[0] = 0 to t[n] = T = .. years
T = 36
#number of time steps in the analysis: dt = T/n
n = 3*10**6
#number of points in the plots. For best result: n/pp should be a
#number without remaider. If we for example choose pp = 5000, choose
#n = frac*5000, where frac 2,3,4, ... 
pp = 5000

#Mass of Sun, Earth, Jupiter [kg]
M_sun = 2*(10**30)
M_E = 6*(10**24)
M_J = 1.9*(10**27)*mf


#Relative mass of Sun, Earth, Jupiter
M = 1.0*np.array([M_sun, M_E, M_J])
M = M/float(M_sun)



#Initial position (in Astronomical units AU) of Earth and Jupiter relative to  
#Solar System Barycenter at 
# 2019-Dec-31 00:00. For now, we let the Sun's initial position be 
#x0 = y0 = 0, where it stays throughout the numerical analysis!

#x0 = np.array([0, -1.528711301502854*0.1, 5.149322540479010*0.1])
#y0 = np.array([0, 9.793926255065635*0.1, -5.194692596109877])


#Initial position (in Astronomical units AU) of Earth and Jupiter relative to  
#Solar System Barycenter at 
# 2019-Nov-29 00:00. For now, we let the Sun's initial position be 
#x0 = y0 = 0, where it stays throughout the numerical analysis!

x0 = np.array([0, 3.948527228009325*0.1, 2.771209156933313*0.1])
y0 = np.array([0, 9.100160380472437*0.1, -5.224508231691265])


#-----------------------------------------------------------------------------------


#Initial velocity (AU/year) of Earth and Jupiter relative to Solar System Barycenter
# 2019-Dec-31 00:00. For now, we let the initial velocity of the Sun be vx = vy = 0,
#and its velocity remains zero throughout the numerical analysis!

#vx0 = np.array([0, -1.729804323045028*(10**(-2)), 7.417650642179744*(10**(-3))])
#vx0 = vx0*365.25
#vy0 = np.array([0, -2.680908671655300*(10**(-3)),  1.104560540984019*(10**(-3))])
#vy0 = vy0*365.25


#Initial velocity (AU/year) of Earth and Jupiter relative to Solar System Barycenter
# 2019-Nov-29 00:00. For now, we let the initial velocity of the Sun be vx = vy = 0,
#and its velocity remains zero throughout the numerical analysis!

vx0 = np.array([0, -1.603566066362973*(10**(-2)), 7.443826199722129*(10**(-3))])
vx0 = vx0*365.25
vy0 = np.array([0, 6.880628606437826*(10**(-3)),  7.587138148929383*(10**(-4))])
vy0 = vy0*365.25



#Sun Object/Instance in the Planet class
Sun = Planet(M[0],x0[0],vx0[0],y0[0],vy0[0],n)
#Earth object in the Planet class
Earth = Planet(M[1],x0[1],vx0[1],y0[1],vy0[1],n)
#Jupiter Object in the Planet class
Jupiter = Planet(M[2],x0[2],vx0[2],y0[2],vy0[2],n)

#List of the objects in our system
Planets = [Sun, Earth, Jupiter]

#Object/Instance in the System class, representing our Sun, Earth, Jupiter system
Alle = System(Planets,T,n)





#Calculation of the total energy of the planet system:
#(Total energy = kinetic energy + potential energy)
#The sun is at rest in origo and therefore has no kinetic or potential energy 

E_TOT = 0
for j in range(1,len(Planets)):
    vel2 =  Planets[j].vel[:,0]**2 + Planets[j].vel[:,1]**2
    E_kin = 0.5*Planets[j].M*vel2
    pos2 = np.sqrt(Planets[j].pos[:,0]**2 + Planets[j].pos[:,1]**2)
    E_pot = -4*np.pi**2*Planets[j].M/pos2
    #Array with total energy for each object in the planet system, except the sun
    E_tot = E_kin+E_pot
    #Array with total energy for the whole planet system
    E_TOT = E_TOT + E_tot

#Max and min total energy
E_MAX = np.max(E_TOT)
E_MIN = np.min(E_TOT)
#print E_MAX
#print E_MIN


#Maximum relative difference in the total energy of the system
dE_TOT = (E_MAX - E_MIN)/(np.abs(np.mean(E_TOT)))


#Calculation of angular momentum for our  planet system
#The sun is at rest in origo and has no angular momentum
ANGM = 0
for j in range(1,len(Planets)):
    #Array with angular momentum for each object in the planet system
    #except the sun
    angm1 = Planets[j].pos[:,0]*Planets[j].vel[:,1]
    angm2 = Planets[j].pos[:,1]*Planets[j].vel[:,0]
    angm = Planets[j].M*(angm1 - angm2)
    #Array with total angular momentum for the whole planet system
    ANGM = ANGM + angm

#Max and min angular momentum
ANGM_max = np.max(ANGM)
ANGM_min = np.min(ANGM)

#Maximum relative difference in Earth's angular momentum
dE_ANGM = (ANGM_max - ANGM_min)/np.mean(ANGM)

#----------------------------------------------------------------------
#In order to speed up the making of the orbit-plots, not all calculated 
#planet positions need to be plotted, we choose for example pp = 5000 points 
#for each of the planet orbits.

x_E   = np.zeros(pp+1)   #Earth orbit. Sun fixed in x=y=0
y_E   = np.zeros(pp+1)
x_J   = np.zeros(pp+1)   #Jupiter orbit. Sun fixed in x=y=0
y_J   = np.zeros(pp+1)

frac = (n+1)//pp
if frac*pp > n:
    pp = pp+1
    frac = (n+1)//pp
for i in range(np.size(x_E)):
    x_E[i] = Alle.Planets[1].pos[frac*i,0]
    y_E[i] = Alle.Planets[1].pos[frac*i,1]
    x_J[i] = Alle.Planets[2].pos[frac*i,0]
    y_J[i] = Alle.Planets[2].pos[frac*i,1]


fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_E, y_E, color = 'r')
ax.plot(x_J, y_J, color = 'b')
#Increase margins on axes
ax.set_xmargin(0.1)
ax.axis('equal')
#plt.axis('equal')
ax.set_xlabel('x(t) [AU]', fontsize = 15)
ax.set_ylabel('y(t) [AU]', fontsize = 15)
ax.set_title('Planet Earth and Jupiter', fontsize = 16) 
#ax.legend(loc='center', fontsize = 14)
ax.tick_params(labelsize = 14)
plt.show()

#-------------------------------------------------------------------


print("""Maximum relative difference in total energy for the system with""")
print("""n = %d and Jupiter mass = %.1E : %.4E""" %(n, M_J, dE_TOT)) 
print(""" """)
print("""Maximum relative difference in angular momentum for the system with""")
print("""n = %d and Jupiter mass = %.1E : %.4E""" %(n, M_J, dE_ANGM)) 
print(""" """)




#---------------------------------------------------

#--Running the script with Jupiter mass M_J = 1.9*(10**27)*1000.
#--and initial conditions from  2019-Nov-29 00:00.

#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5e.py
#Maximum relative difference in total energy for the system with
#n = 3000000 and Jupiter mass = 1.9E+30 : 2.3971E-02

#Maximum relative difference in angular momentum for the system with
#n = 3000000 and Jupiter mass = 1.9E+30 : 1.5515E-13



#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5e.py





#--Running the script with Jupiter mass M_J = 1.9*(10**27)*950.
#--and initial conditions from  2019-Nov-29 00:00.

#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5e.py
#Maximum relative difference in total energy for the system with
#n = 1000000 and Jupiter mass = 1.8E+30 : 6.7459E-04

#Maximum relative difference in angular momentum for the system with
#n = 1000000 and Jupiter mass = 1.8E+30 : 9.3579E-14



#--Running the script with Jupiter mass M_J = 1.9*(10**27)*1000.
#--and initial conditions from  2019-Dec-31 00:00.

#(base) M:FYS4150-H19\Prosjekt5\Programmer>python Project5e.py
#Maximum relative difference in total energy for the system with
#n = 3000000 and Jupiter mass = 1.9E+30 : 2.4992E-02
#
#Maximum relative difference in angular momentum for the system with
#n = 3000000 and Jupiter mass = 1.9E+30 : 2.5085E-13


#(base) M:\FYS4150-H19\Prosjekt5\Programmer>




#Testprints
#print Solarsystem[1].M
#print Earth.y0
#print Earth.vel[0,0]
#print Alle.Planets[0].M
#print Alle.leng
#print Alle.Planets[1].pos



