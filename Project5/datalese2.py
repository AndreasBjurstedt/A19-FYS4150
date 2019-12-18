
#Reading result files made by Project5g.py and making
#plot of Mercury perihelion percession.
import matplotlib.pyplot as plt
import numpy as np
from math import *
import array as arr



def fileFormat(file):
    f1=open(file, "r")
    lines =f1.readlines()
    f1.close()
    kolonner=0
    rader=0

    for line in lines:

        rader=rader+1
        if(rader==1):
            words=line.split(" ")
            for word in words:
                kolonner =kolonner+1



    return rader,kolonner

def readFile(file):
    rader,kolonner=fileFormat(file)
    A=np.zeros((rader,kolonner))
    f=open(file, "r")
    lines =f.readlines()
    f.close()
    for i in range(rader):
        words=lines[i].split(" ")
        for j in range(kolonner):
            A[i][j]=float(words[j])

    return A,rader,kolonner

data,rader,kolonner=readFile('data5g_E-100y-10xx6.txt')
data_N,rader_N,kolonner_N=readFile('data5g_N-100y-10xx6.txt')
theta=np.zeros(rader)
theta_N=np.zeros(rader_N)
for i in range(rader):
    theta[i]=3600*(180/np.pi)*np.arctan(data[i][2]/data[i][1])
    theta_N[i]=3600*(180/np.pi)*np.arctan(data_N[i][2]/data_N[i][1])
dataT=data.T
dataT_N = data_N.T

t_eksakt = np.array([-5,105])
theta_eksakt= np.array([43,43])

plt.figure()
plt.plot(dataT[0],theta,'-r')
plt.plot(dataT_N[0],theta_N,'-b')
plt.plot(t_eksakt,theta_eksakt,'-g')

plt.xlabel('Time [year]',fontsize = 14)
plt.ylabel('Perihelion angle [arc seconds] ',fontsize = 14)
plt.legend(['Einstein', 'Newton', 'observed (1 century)'])
plt.title('Perihelion precession of Mercury', fontsize = 14)
plt.tight_layout()
#plt.ylim((-798.60,-799))
plt.savefig('merk_theta.png')
