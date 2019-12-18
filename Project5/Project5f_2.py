
#Using the classes in Project5_classes.py to calculate the planet orbits.
#These orbits are compared in a plot. 
#
#Also included for analysis with common mass center in origo:
#      Maximum relative difference in total energy for the system
#      Maximum relative difference in angular momentum for the system
#      Maximum relative difference in total energy for the system when excluding the
#      contribution from the Sun.

import numpy as np
import matplotlib.pyplot as plt
from Project5_classes import *


#Time span of numerical analysis: t[0] = 0 to t[n] = T = 165 years
T = 165
#number of time steps in the analysis: dt = T/n
n = 150000
#number of points in the plots. For best result: n/pp should be a
#number without remaider. If we for example choose pp = 5000, choose
#n = frac*5000, where frac 2,3,4, ... 
pp = 5000

#Mass of Sun, Earth, Jupiter and the other planets in the solarsystem [kg]
M_sun = 2*(10**30)           #The Sun
M_Me = 3.3*(10**23)          #Mercury
M_Ve = 4.9*(10**24)          #Venus
M_E = 6*(10**24)             #Earth
M_Ma = 6.6*(10**23)          #Mars
M_J = 1.9*(10**27)           #Jupiter
M_Sa = 5.5*(10**26)          #Saturn
M_U = 8.8*(10**25)           #Uranus
M_N =1.03*(10**26)           #Neptune

#Relative mass of Sun, Earth, Jupiter and the other planets
M_Me  = M_Me/float(M_sun)
M_Ve  = M_Ve/float(M_sun)
M_E   = M_E/float(M_sun)
M_Ma  = M_Ma/float(M_sun)
M_J   = M_J/float(M_sun)
M_Sa  = M_Sa/float(M_sun)
M_U   = M_U/float(M_sun)
M_N   = M_N/float(M_sun)
M_sun = M_sun/float(M_sun)



#Initial position (in Astronomical units AU) of Earth, Jupiter and the other planets
#relative to the Solar System Barycenter (Mass center) at 2019-Nov-29 00:00.
x0_Me = -3.089137495084154*0.1   #Mercury
x0_Ve = 4.814605067455450*0.1    #Venus
x0_E =  3.948527228009325*0.1    #Earth
x0_Ma = -1.543208932754952       #Mars
x0_J =  2.771209156933313*0.1    #Jupiter
x0_Sa = 3.632628879585697        #Saturn
x0_U = 1.629688404988837*10      #Uranus
x0_N = 2.921750763559268*10      #Neptune


y0_Me = 1.744886373010318*0.1       #Mercury 
y0_Ve = -5.345470370402726*0.1      #Venus
y0_E =  9.100160380472437*0.1       #Earth
y0_Ma = -5.042083307102040*0.1      #Mars
y0_J = -5.224508231691265           #Jupiter
y0_Sa = -9.348288811274543          #Saturn
y0_U = 1.128605542338266*10         #Uranus
y0_N = -6.461552366128481           #Neptune

#Initial velocity (AU/year) of Earth,Jupiter and the other planets relative to Solar 
#System Barycenter (Mass center) 2019-Nov-29 00:00.

vx0_Me = -1.928258980107407*0.01        #Mercury
vx0_Ve = 1.493272115673404*0.01         #Venus
vx0_E =  -1.603566066362973*0.01        #Earth
vx0_Ma = 4.926983218616618*10**(-3)     #Mars
vx0_J = 7.443826199722129*10**(-3)      #Jupiter
vx0_Sa = 4.891570166847385*10**(-3)     #Saturn
vx0_U = -2.268171563746997*10**(-3)     #Uranus
vx0_N = 6.568727514842341*10**(-4)      #Neptune

vx0_Me = vx0_Me*365.25
vx0_Ve = vx0_Ve*365.25
vx0_E  = vx0_E*365.25
vx0_Ma = vx0_Ma*365.25
vx0_J  = vx0_J*365.25
vx0_Sa = vx0_Sa*365.25
vx0_U  = vx0_U*365.25
vx0_N  = vx0_N*365.25


vy0_Me = -2.350312105925493*0.01         #Mercury 
vy0_Ve = 1.341199462523215*0.01          #Venus
vy0_E = 6.880628606437826*10**(-3)       #Earth
vy0_Ma = -1.208451455394788*0.01         #Mars
vy0_J =   7.587138148929383*10**(-4)     #Jupiter
vy0_Sa = 2.004764720827292*10**(-3)      #Saturn
vy0_U = 3.050128900966884*10**(-3)       #Uranus
vy0_N = 3.084041276878159*10**(-3)       #Neptune

vy0_Me = vy0_Me*365.25
vy0_Ve = vy0_Ve*365.25
vy0_E  = vy0_E*365.25
vy0_Ma = vy0_Ma*365.25
vy0_J  = vy0_J*365.25
vy0_Sa = vy0_Sa*365.25
vy0_U  = vy0_U*365.25
vy0_N  = vy0_N*365.25


#Mercury object in the Planet class
Mercury = Planet(M_Me,x0_Me,vx0_Me,y0_Me,vy0_Me,n)
#Venus Object in the Planet class
Venus = Planet(M_Ve,x0_Ve,vx0_Ve,y0_Ve,vy0_Ve,n)
#Earth Object/Instance in the Planet class
Earth = Planet(M_E,x0_E,vx0_E,y0_E,vy0_E,n)
#Mars Object/Instance in the Planet class
Mars = Planet(M_Ma,x0_Ma,vx0_Ma,y0_Ma,vy0_Ma,n)
#Jupiter object in the Planet class
Jupiter = Planet(M_J,x0_J,vx0_J,y0_J,vy0_J,n)
#Saturn Object in the Planet class
Saturn = Planet(M_Sa,x0_Sa,vx0_Sa,y0_Sa,vy0_Sa,n)
#Uranus object in the Planet class
Uranus = Planet(M_U,x0_U,vx0_U,y0_U,vy0_U,n)
#Neptune Object in the Planet class
Neptune = Planet(M_N,x0_N,vx0_N,y0_N,vy0_N,n)


#Array with relative masses of the Sun, Earth, Jupiter and the other
#planets
M = 1.0*np.array([M_sun, M_Me, M_Ve, M_E, M_Ma, M_J, M_Sa, M_U, M_N])

#Array with initial positions (in Astronomical units AU) of Earth, Jupiter and 
#the other planets relative to the Solar System Barycenter (Mass center) 
#at 2018-Nov-18 00:00.
x0 = np.array([0.0, x0_Me, x0_Ve, x0_E, x0_Ma, x0_J, x0_Sa, x0_U, x0_N])
y0 = np.array([0.0, y0_Me, y0_Ve, y0_E, y0_Ma, y0_J, y0_Sa, y0_U, y0_N])



#The initial position of the Sun relative to the Mass center of the system
Mxx = 0
Myy = 0
for i in range(1,np.size(x0)):
    Mx = M[i]*x0[i]
    My = M[i]*y0[i]
    Mxx = Mxx + Mx
    Myy = Myy + My
x0[0] = -Mxx/float(M[0])
y0[0] = -Myy/float(M[0])

#print x0
#print y0

############Test of mass center: that it is in x = y = 0.########
Mxx = 0
Myy = 0
Mx = 0
My = 0

for i in range(np.size(x0)):
    Mx = M[i]*x0[i]
    My = M[i]*y0[i]
    Mxx = Mxx + Mx
    Myy = Myy + My

#print Mxx
#print Mxx
###################

#Array with initial velocities (AU/year) of Earth, Jupiter and 
#the other planets relative to the Solar System Barycenter (Mass center) 

vx0 = np.array([0.0, vx0_Me, vx0_Ve, vx0_E, vx0_Ma, vx0_J, vx0_Sa, vx0_U, vx0_N])
vy0 = np.array([0.0, vy0_Me, vy0_Ve, vy0_E, vy0_Ma, vy0_J, vy0_Sa, vy0_U, vy0_N])

#vx0 = np.array([0.0, vx0_Ve, vx0_E, vx0_J])
#vy0 = np.array([0.0, vy0_Ve, vy0_E, vy0_J])

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

###########Testing that the systems momentum is zero##########
Vxx = 0
Vyy = 0
Vx = 0
Vy = 0

for i in range(np.size(x0)):
    Vx = M[i]*vx0[i]
    Vy = M[i]*vy0[i]
    Vxx = Vxx + Vx
    Vyy = Vyy + Vy

#print Vxx
#print Vyy
###################

#Sun Object/Instance in the Planet class
Sun = Planet(M_sun,x0[0],vx0[0],y0[0],vy0[0],n)

#List of the objects in our system
Planets = [Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune]

#Object/Instance in the System class, representing our Sun, Earth, Jupiter system
Alle = System(Planets,T,n)


#------------------------------------------------------------------------

#Calculation of the total energy of the planet system and the sun
#in orbit around their common mass center:
#(Total energy = kinetic energy + potential energy)

E_TOT = 0
for j in range(0,len(Alle.Planets)):
    vel2 =  Alle.Planets[j].vel[:,0]**2 + Alle.Planets[j].vel[:,1]**2
    E_kin = 0.5*Alle.Planets[j].M*vel2
    pos2 = np.sqrt(Alle.Planets[j].pos[:,0]**2 + Alle.Planets[j].pos[:,1]**2)
    E_pot = -4*np.pi**2*Alle.Planets[j].M/pos2
    #Array with total energy for each object in the planet system, including the sun
    E_tot = E_kin+E_pot
    #Array with total energy for the whole planet system including the sun
    E_TOT = E_TOT + E_tot

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
for j in range(1,len(Alle.Planets)):
    vel2 =  Alle.Planets[j].vel[:,0]**2 + Alle.Planets[j].vel[:,1]**2
    E_kin = 0.5*Alle.Planets[j].M*vel2
    pos2 = np.sqrt(Alle.Planets[j].pos[:,0]**2 + Alle.Planets[j].pos[:,1]**2)
    E_pot = -4*np.pi**2*Alle.Planets[j].M/pos2
    #Array with total energy for each object in the planet system, the sun excluded
    E_tot = E_kin+E_pot
    #Array with total energy for the whole planet system including, the sun excluded
    E_TOT = E_TOT + E_tot

#Max and min total energy
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
for j in range(len(Alle.Planets)):
    #Array with angular momentum for each object in the planet system
    #except the sun
    angm1 = Alle.Planets[j].pos[:,0]*Alle.Planets[j].vel[:,1]
    angm2 = Alle.Planets[j].pos[:,1]*Alle.Planets[j].vel[:,0]
    angm = Alle.Planets[j].M*(angm1 - angm2)
    #Array with total angular momentum for the whole planet system
    ANGM = ANGM + angm

#Max and min angular momentum
ANGM_max = np.max(ANGM)
ANGM_min = np.min(ANGM)

#Maximum relative difference in the system's angular momentum
dE_ANGM = (ANGM_max - ANGM_min)/np.mean(ANGM)


#---------------------------------------------------------------------


#In order to speed up the making of the orbit-plots, not all calculated 
#planet positions need to be plotted, we choose for example pp = 5000 points 
#for each of the planet orbits.

#aa = np.zeros(pp+1)
#x_sun = np.zeros(pp+1)    #Sun
#y_sun = np.zeros(pp+1)
x_Me  = np.zeros(pp+1)    #Mercury
y_Me  = np.zeros(pp+1) 
x_Ve  = np.zeros(pp+1)    #Venus
y_Ve  = np.zeros(pp+1)
x_E   = np.zeros(pp+1)    #Earth
y_E   = np.zeros(pp+1)
x_Ma  = np.zeros(pp+1)    #Mars
y_Ma  = np.zeros(pp+1)
x_J   = np.zeros(pp+1)    #Jupiter
y_J   = np.zeros(pp+1)
x_Sa  = np.zeros(pp+1)    #Saturn
y_Sa  = np.zeros(pp+1)
x_U   = np.zeros(pp+1)    #Uranus
y_U   = np.zeros(pp+1)
x_N   = np.zeros(pp+1)    #Neptune
y_N   = np.zeros(pp+1)


frac = (n+1)//pp
if frac*pp > n:
    pp = pp+1
    frac = (n+1)//pp
for i in range(np.size(x_Me)):
    x_Me[i] = Alle.Planets[1].pos[frac*i,0]
    y_Me[i] = Alle.Planets[1].pos[frac*i,1]
    x_Ve[i] = Alle.Planets[2].pos[frac*i,0]
    y_Ve[i] = Alle.Planets[2].pos[frac*i,1]
    x_E[i] = Alle.Planets[3].pos[frac*i,0]
    y_E[i] = Alle.Planets[3].pos[frac*i,1]
    x_Ma[i] = Alle.Planets[4].pos[frac*i,0]
    y_Ma[i] = Alle.Planets[4].pos[frac*i,1]
    x_J[i] = Alle.Planets[5].pos[frac*i,0]
    y_J[i] = Alle.Planets[5].pos[frac*i,1]
    x_Sa[i] = Alle.Planets[6].pos[frac*i,0]
    y_Sa[i] = Alle.Planets[6].pos[frac*i,1]
    x_U[i] = Alle.Planets[7].pos[frac*i,0]
    y_U[i] = Alle.Planets[7].pos[frac*i,1]
    x_N[i] = Alle.Planets[8].pos[frac*i,0]
    y_N[i] = Alle.Planets[8].pos[frac*i,1]


#Plot of the planet orbits
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_Me,y_Me,color = 'b')
ax.plot(x_Ve,y_Ve,color = 'b')
ax.plot(x_E,y_E,color = 'r')
ax.plot(x_Ma,y_Ma,color = 'b')
ax.plot(x_J,y_J,color = 'b')
ax.plot(x_Sa,y_Sa,color = 'b')
ax.plot(x_U,y_U,color = 'b')
ax.plot(x_N,y_N,color = 'b')

##ax.plot(Alle.Planets[1].pos[:,0], Alle.Planets[1].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[2].pos[:,0], Alle.Planets[2].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[3].pos[:,0], Alle.Planets[3].pos[:,1], color = 'r')
##ax.plot(Alle.Planets[4].pos[:,0], Alle.Planets[4].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[5].pos[:,0], Alle.Planets[5].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[6].pos[:,0], Alle.Planets[6].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[7].pos[:,0], Alle.Planets[7].pos[:,1], color = 'b')
##ax.plot(Alle.Planets[8].pos[:,0], Alle.Planets[8].pos[:,1], color = 'b')
#Increase margins on axes
#ax.set_xmargin(0.1)
#ax.axis('equal')
plt.axis('equal')
plt.axis([min(x_N)-1, max(x_N)+1, min(y_N)-3, max(y_N)+3])
ax.set_xlabel('x(t) [AU]', fontsize = 15)
ax.set_ylabel('y(t) [AU]', fontsize = 15)
ax.set_title('The planetary orbits around the Mass center of the system', fontsize = 16) 
#ax.legend(loc='center', fontsize = 14)
ax.tick_params(labelsize = 14)
plt.show()

#Plot of the planet orbits from Mercury to Mars
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x_Me,y_Me,color = 'b')
ax.plot(x_Ve,y_Ve,color = 'b')
ax.plot(x_E,y_E,color = 'r')
ax.plot(x_Ma,y_Ma,color = 'b')

#Increase margins on axes
#ax.set_xmargin(0.1)
#ax.axis('equal')
plt.axis('equal')
plt.axis([min(x_Ma)-1, max(x_Ma)+1, min(y_Ma)-1, max(y_Ma)+1])
ax.set_xlabel('x(t) [AU]', fontsize = 15)
ax.set_ylabel('y(t) [AU]', fontsize = 15)
ax.set_title('The orbits of the planets from Mercury to Mars', fontsize = 16) 
#ax.legend(loc='center', fontsize = 14)
ax.tick_params(labelsize = 14)
plt.show()



#Printing results
print(""" """)
print("""Number of time steps in analysis: n = %d."""%(n))
print("""Number of years in analysis: T = %.1f."""%(T))

print(""" """)
print("""Maximum relative difference in total energy for the system with""")
print("""the Sun included:%.4E."""%(dE_TOT1)) 

print(""" """)
print("""Maximum relative difference in total energy for the system with""")
print("""the Sun excluded: %.4E""" %(dE_TOT2)) 
print(""" """)
print("""Maximum relative difference in angular momentum for the system with""")
print("""the Sun included: %.4E""" %(dE_ANGM)) 

#------------------------------------------------------------------------------

#Results from analysis
#(base) M:\FYS4150-H19\Prosjekt5\Programmer>python Project5f_2.py
#
#Number of time steps in analysis: n = 150000.
#Number of years in analysis: T = 165.0.
#
#Maximum relative difference in total energy for the system with
#the Sun included:4.7775E+01.
#
#Maximum relative difference in total energy for the system with
#the Sun excluded: 4.5228E-03
#
#Maximum relative difference in angular momentum for the system with
#the Sun included: 3.3055E-14
#
#(base) M:\FYS4150-H19\Prosjekt5\Programmer>






#Testprints
#print Solarsystem[1].M
#print Earth.y0
#print Earth.vel[0,0]
#print Alle.Planets[0].M
#print Alle.leng
#print Alle.Planets[1].pos



