
#Using the classes in Project5_classes.py to calculate the orbits of planet Earth
#and Jupiter when: 
#      1) The Sun is fixed in Origo
#      2) The Sun, Planet Earth and Jupiter are orbiting their common mass center. 
#These orbits are compared in a plot. Also a plot of the Sun's orbit around the 
#common mass center. 
#
#Also included for analysis with common mass center in origo:
#      Maximum relative difference in total energy for the system
#      Maximum relative difference in angular momentum for the system
#      Maximum relative difference in total energy for the system when excluding the
#      contribution from the Sun.

import numpy as np
import matplotlib.pyplot as plt
from Project5_classes import *

#Time span of numerical analysis: t[0] = 0 to t[n] = T = 12 years
T = 12
#number of time steps in the analysis: dt = T/n
n = 10000
#number of points in the plots. For best result: n/pp should be a
#number without remaider. If we for example choose pp = 5000, choose
#n = frac*5000, where frac 2,3,4, ... 
pp = 5000

#Mass of Sun, Earth, Jupiter [kg]
M_sun = 2*(10**30)
M_E = 6*(10**24)
M_J = 1.9*(10**27)

#Relative mass of Sun, Earth, Jupiter
M = 1.0*np.array([M_sun, M_E, M_J])
M = M/float(M_sun)

#Initial position (in Astronomical units AU) of Earth and Jupiter relative to  
#Solar System Barycenter (Mass center) at 2019-Nov-29 00:00.
#The initial position of the Sun is not calculated yet and is set to x = y = 0
#for now

x0 = np.array([0, 3.948527228009325*0.1, 2.771209156933313*0.1])
y0 = np.array([0, 9.100160380472437*0.1, -5.224508231691265])

vx0 = np.array([0, -1.603566066362973*(10**(-2)), 7.443826199722129*(10**(-3))])
vx0 = vx0*365.25
vy0 = np.array([0, 6.880628606437826*(10**(-3)),  7.587138148929383*(10**(-4))])
vy0 = vy0*365.25






#Function which calculates the orbit of the Earth and Jupiter with
#the Velocity Verlet method, when the position of the Sun is in
#x = y = 0 and its velocity is vx = vy = 0. This is essentially what we did in
#Project5e.py.

def Project5e(T,n,M,x0,y0,vx0,vy0):

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
    
    return Alle

#Function which calculates the orbit of the Sun, the Earth and Jupiter
#around their common mass center

def Project5f(T,n,M,x0,y0,vx0,vy0):

    #The initial position of the Sun relative to  the mass center of the system
    Mxx = 0
    Myy = 0
    for i in range(1,np.size(x0)):
        Mx = M[i]*x0[i]
        My = M[i]*y0[i]
        Mxx = Mxx + Mx
        Myy = Myy + My
    x0[0] = -Mxx/float(M[0])
    y0[0] = -Myy/float(M[0])

    #The initial velocity of the Sun relative to the mass center of the system
    #(The momentum (bevegelsesmengde) of the system is zero.)
    Vxx = 0
    Vyy = 0
    for i in range(1,np.size(x0)):
        Vx = M[i]*vx0[i]
        Vy = M[i]*vy0[i]
        Vxx = Vxx + Vx
        Vyy = Vyy + Vy
    vx0[0] = -Vxx/float(M[0])
    vy0[0] = -Vyy/float(M[0])

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

    return Alle


#Calculating the orbit of the Earth and of Jupiter when the
#position of the Sun is fixed in x = y = 0
Alle_e= Project5e(T,n,M,x0,y0,vx0,vy0)

#Calculating the orbit of the Sun, the Earth and Jupiter 
#around their common mass center
Alle_f = Project5f(T,n,M,x0,y0,vx0,vy0)


#------------------------------------------------------------------------

#Calculation of the total energy of the planet system and the sun
#in orbit around their common mass center:
#(Total energy = kinetic energy + potential energy)

E_TOT = 0
for j in range(0,len(Alle_f.Planets)):
    vel2 =  Alle_f.Planets[j].vel[:,0]**2 + Alle_f.Planets[j].vel[:,1]**2
    E_kin = 0.5*Alle_f.Planets[j].M*vel2
    pos2 = np.sqrt(Alle_f.Planets[j].pos[:,0]**2 + Alle_f.Planets[j].pos[:,1]**2)
    E_pot = -4*np.pi**2*Alle_f.Planets[j].M/pos2
    #Array with total energy for each object in the planet system, including the sun
    E_tot = E_kin+E_pot
    #Array with total energy for the whole planet system including the sun
    E_TOT = E_TOT + E_tot


#E_TOT = E_TOT[np.size(E_TOT)/2:np.size(E_TOT)]
#Max and min total enrgy
E_MAX = np.max(E_TOT)
E_MIN = np.min(E_TOT)

#The position of max and min total energy in the total energy array.
E_MAXP = np.argmax(E_TOT)
E_MINP = np.argmin(E_TOT)

#Maximum relative difference in the total energy of the system
dE_TOT1 = np.abs((E_MAX - E_MIN)/np.mean(E_TOT))

#-----------------------------------------------------------------------------

#Calculation of the total energy of the planet system 
#in orbit around their common mass center:
#(Total energy = kinetic energy + potential energy)
#The total energy due to the motion of the sun is excluded.
E_TOT = 0
for j in range(1,len(Alle_f.Planets)):
    vel2 =  Alle_f.Planets[j].vel[:,0]**2 + Alle_f.Planets[j].vel[:,1]**2
    E_kin = 0.5*Alle_f.Planets[j].M*vel2
    pos2 = np.sqrt(Alle_f.Planets[j].pos[:,0]**2 + Alle_f.Planets[j].pos[:,1]**2)
    E_pot = -4*np.pi**2*Alle_f.Planets[j].M/pos2
    #Array with total energy for each object in the planet system, the sun excluded
    E_tot = E_kin+E_pot
    #Array with total energy for the whole planet system including, the sun excluded
    E_TOT = E_TOT + E_tot

#Max and min total enrgy
E_MAX = np.max(E_TOT)
E_MIN = np.min(E_TOT)
#The position of max and min total energy in the total energy array.
E_MAXP = np.argmax(E_TOT)
E_MINP = np.argmin(E_TOT)

#Maximum relative difference in the total energy of the system
dE_TOT2 = np.abs((E_MAX - E_MIN)/np.mean(E_TOT))

#-------------------------------------------------------------------------------

#Calculation of angular momentum for our  planet system, including the sun
#in orbit around the common mass center of the system.
ANGM = 0
for j in range(len(Alle_f.Planets)):
    #Array with angular momentum for each object in the planet system
    #except the sun
    angm1 = Alle_f.Planets[j].pos[:,0]*Alle_f.Planets[j].vel[:,1]
    angm2 = Alle_f.Planets[j].pos[:,1]*Alle_f.Planets[j].vel[:,0]
    angm = Alle_f.Planets[j].M*(angm1 - angm2)
    #Array with total angular momentum for the whole planet system
    ANGM = ANGM + angm

#Max and min angular momentum
ANGM_max = np.max(ANGM)
ANGM_min = np.min(ANGM)

#Maximum relative difference in the system's angular momentum
dE_ANGM = (ANGM_max - ANGM_min)/np.mean(ANGM)


#---------------------------------------------------------------------
#Calculating the maximum relative difference in the position of the
#Earth, comparing its orbit when the Sun is fixed in origo with
#the orbit when the Sun, the Earth and Jupiter are orbiting around 
#their common mass center. The average distance between the Sun and planet
#Earth is 1 Au.

#The calculated positions of planet Earth, the sun fixed in origo
x_Ee = Alle_e.Planets[1].pos[:,0] 
y_Ee = Alle_e.Planets[1].pos[:,1]

#The calculated positions of planet Earth when it orbits around the common
#mass center.
x_Ef = Alle_f.Planets[1].pos[:,0] 
y_Ef = Alle_f.Planets[1].pos[:,1]

r_Ee = np.sqrt(x_Ee**2 + y_Ee**2)
r_Ef = np.sqrt(x_Ef**2 + y_Ef**2)

diff_E = np.abs(r_Ee - r_Ef)
M_diff_E = np.max(diff_E)

#print M_diff_E

#----------------------------------------------------------------------

#Calculating the maximum relative difference in the position of 
#Jupiter, comparing its orbit when the Sun is fixed in origo with
#the orbit when the Sun, the Earth and Jupiter are orbiting around 
#their common mass center. The average distance between the Sun and
#Jupiter is 5.2 Au.

#The calculated positions of Jupiter, the sun fixed in origo
x_Je = Alle_e.Planets[2].pos[:,0] 
y_Je = Alle_e.Planets[2].pos[:,1]

#The calculated positions of Jupiter when it orbits around the common
#mass center.
x_Jf = Alle_f.Planets[2].pos[:,0] 
y_Jf = Alle_f.Planets[2].pos[:,1]

r_Je = np.sqrt(x_Je**2 + y_Je**2)
r_Jf = np.sqrt(x_Jf**2 + y_Jf**2)

diff_J = np.abs(r_Je - r_Jf)
M_diff_J = np.max(diff_J)/5.2
 
#print M_diff_J

#-----------------------------------------------------------------------------------

#In order to speed up the making of the orbit-plots, not all calculated 
#planet positions need to be plotted, we choose for example pp = 5000 points 
#for each of the planet orbits.

x_suf = np.zeros(pp+1)    #Sun in orbit around common mass center
y_suf = np.zeros(pp+1)
x_Ee   = np.zeros(pp+1)   #Earth orbit. Sun fixed in x=y=0
y_Ee   = np.zeros(pp+1)
x_Ef  = np.zeros(pp+1)    #Earth orbit. Sun also in orbit around
y_Ef  = np.zeros(pp+1)    #common mass center
x_Je   = np.zeros(pp+1)   #Jupiter orbit. Sun fixed in x=y=0
y_Je   = np.zeros(pp+1)
x_Jf  = np.zeros(pp+1)    #Jupiter orbit. Sun also in orbit around
y_Jf  = np.zeros(pp+1)    #common mass center


frac = (n+1)//pp
if frac*pp > n:
    pp = pp+1
    frac = (n+1)//pp
for i in range(np.size(x_Ee)):
    x_suf[i] = Alle_f.Planets[0].pos[frac*i,0] 
    y_suf[i] = Alle_f.Planets[0].pos[frac*i,1]
    x_Ee[i] = Alle_e.Planets[1].pos[frac*i,0]
    y_Ee[i] = Alle_e.Planets[1].pos[frac*i,1]
    x_Ef[i] = Alle_f.Planets[1].pos[frac*i,0]
    y_Ef[i] = Alle_f.Planets[1].pos[frac*i,1]
    x_Je[i] = Alle_e.Planets[2].pos[frac*i,0]
    y_Je[i] = Alle_e.Planets[2].pos[frac*i,1]
    x_Jf[i] = Alle_f.Planets[2].pos[frac*i,0]
    y_Jf[i] = Alle_f.Planets[2].pos[frac*i,1]
    
#Plot of orbit of planet Earth and Jupiter. --- lines are for the orbits
#when the Sun is fixed in origo of our coordinate system.
#Whole lines are for the orbits when the common mass center are in origo
#of our coordinate system
fig = plt.figure()
ax = plt.subplot(111)
#ax.plot(x_Ef, y_Ef, color = 'b')
ax.plot(x_Ef, y_Ef, 'b', label = 'Common Mass center in origo')
ax.plot(x_Ee, y_Ee, 'g--',label = 'The Sun in origo')
ax.plot(x_Jf, y_Jf, 'b')
ax.plot(x_Je, y_Je, 'g--')

#Increase margins on axes
#ax.set_xmargin(0.1)
#ax.axis('equal')
plt.axis('equal')
plt.axis([min(x_Jf)-0.5, max(x_Jf)+0.5, min(y_Jf)-0.5, max(y_Jf)+2.2])
ax.set_xlabel('x(t) [AU]', fontsize = 15)
ax.set_ylabel('y(t) [AU]', fontsize = 15)
ax.set_title('The orbit of planet Earth and Jupiter', fontsize = 16) 
ax.legend(loc='upper right', fontsize = 14)
ax.tick_params(labelsize = 14)
plt.show()


#Plot of the orbit of the Sun around the common mass center.
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_suf, y_suf, 'b')
#Increase margins on axes
#ax.set_xmargin(0.1)
#ax.axis('equal')
plt.axis('equal')
plt.axis([min(x_suf)-0.0, max(x_suf)+0.0, min(y_suf)-0.0005, max(y_suf)+0.0005])
ax.set_xlabel('x(t) [AU]', fontsize = 15)
ax.set_ylabel('y(t) [AU]', fontsize = 15)
ax.set_title('The orbit of the Sun around common Mass center', fontsize = 16) 
ax.tick_params(labelsize = 14)
plt.show()


print(""" """)
print("""Number of time steps in analysis: n = %d."""%(n))
print("""Number of years in analysis: T = %d."""%(T))

print(""" """)
print("""Maximum relative difference in total energy for the system with""")
print("""the Sun included:%.4E."""%(dE_TOT1)) 

print(""" """)
print("""Maximum relative difference in total energy for the system with""")
print("""the Sun excluded: %.4E""" %(dE_TOT2)) 
print(""" """)
print("""Maximum relative difference in angular momentum for the system with""")
print("""the Sun included: %.4E""" %(dE_ANGM)) 


#-------------------------------------------------------------------
#Running the script with the Sun, planet Earth and Jupiter orbiting around their
#common mass center


#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5f.py
#
#Number of time steps in analysis: n = 10000.
#Number of years in analysis: T = 12.
#
#Maximum relative difference in total energy for the system with
#the Sun included:9.9464E-02.
#
#Maximum relative difference in total energy for the system with
#the Sun excluded: 7.0705E-04
#
#Maximum relative difference in angular momentum for the system with
#the Sun included: 1.4414E-14
#
#(base) M:\FYS4150-H19\Prosjekt5\Programmer>














