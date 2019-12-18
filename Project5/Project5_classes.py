

import numpy as np

class Planet:
    def __init__(self,M,x0,vx0,y0,vy0,n):
        self.M = M         #Mass
        self.x0 = x0       #Position in x direction at t = t0
        self.y0 = y0       #Position in y direction at t = t0
        self.vx0 = vx0     #Velocity in x direction at t = t0
        self.vy0 = vy0     #Velocity in y direction at t = t0
        self.n = n         #number of time steps

        #Matrix pos for x and y component of position
        #Matrix vel for x and y component of speed
        #Matrix acc for x and y component of acceleration
        self.pos = np.zeros([n+1,2])
        self.vel = np.zeros([n+1,2])
        self.acc = np.zeros([n+1,2])

        #Initial values of position and speed in x and y direction 
        #when t = t0
        self.pos[0,:] = np.array([float(x0),float(y0)])
        self.vel[0,:] = np.array([float(vx0),float(vy0)])

        

class System:
    def __init__(self,Planets,T,n):
        self.n = n                #Number of time steps
        self.T = T                #Last time step  T are in number of years
        self.Planets = Planets    #List of the Sun and the planets in the analysis
        n, T, Planets = self.n, self.T, self.Planets

        #Time step between time-grid points
        self.dt = T/float(n)
        dt = self.dt

        #Time vector t[0],t[1],t[2],......,t[n]
        self.t = np.zeros(n+1)
        t = self.t

        self.leng = len(Planets)
        leng = self.leng

        #Function which calculates the acceleration of the sun and all the planets. The
        #acceleration for one object is due to the gravitational pull from all the other objects.
        #Acceleration unit: AU/(year**2)
        def func(i,j):
            leng = self.leng
            ff = np.array([0,0])
            for k in range(leng):
                #An object does not experience a gravitational pull from itself
                if k == j:
                    continue
                    #ff = ff + np.array([0,0]) 
                #The acceleration of each object due to the gravitational pull from all the other
                #objects  
                else: 
                    rx = Planets[j].pos[i,0]-Planets[k].pos[i,0]
                    ry = Planets[j].pos[i,1]-Planets[k].pos[i,1]
                    r = np.sqrt((rx)**2 + (ry)**2)
                    f1 = -4*(np.pi**2)*(Planets[k].M/float(Planets[0].M))/float(r**3)  
                    f2 = np.array([rx,ry])
                    f = f1*f2
                    ff = ff + f
            return ff

        
     
        #The Velocity verlet method
        for i in range(n):
            t[i+1] = t[i] + dt
            #All the planet positions have to be updated before their velocities can be updated
            #(Planet velocities are dependent on updated planet accelerations from the function func,
            #which in turn is dependent of updated planet positions.)

            #Also: If the initial position of the Sun is in origo (x0=y0=0) we let it stay there
            #In a real system, the systems Center of Mass should be positioned
            #in origo and the Sun a little bit outside it. The same applies for the initial velocity
            #of the sun. If it is zero, it stays zero. In a real system the initial velocity of the
            #Sun makes the total momentum of the system exactly zero.
            for j in range(leng):
                if j == 0 and Planets[0].pos[i,0] == 0 and Planets[0].pos[i,1] == 0:
                    Planets[j].pos[i+1] = [0,0]
                else:
                    Planets[j].pos[i+1] = Planets[j].pos[i] + dt*Planets[j].vel[i] + 0.5*(dt**2)*func(i,j)
            for j in range(leng):
                if j == 0 and Planets[0].vel[i,0] == 0 and Planets[0].vel[i,1] == 0:
                    Planets[j].vel[i+1] = [0,0]
                else:
                    Planets[j].vel[i+1] = Planets[j].vel[i] + 0.5*dt*(func(i+1,j) + func(i,j))


            #for j in range(leng):
            #    Planets[j].pos[i+1] = Planets[j].pos[i] + dt*Planets[j].vel[i] + 0.5*dt**2*func(i,j)
            #for j in range(leng):
            #    Planets[j].vel[i+1] = Planets[j].vel[i] + 0.5*dt*(func(i+1,j) + func(i,j))







                


       






            


















    




