#!/usr/bin/python

#este programa calculara la trayectoria de 2 cuerpos

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os


def twob(y, t, params):
    '''
    funcion para calcular la posicion y velocidad de las particulas
    '''
    mu = params['mu']
    r = y[0:3]
    v = y[3:6]
    drdt = v
    r2 = y[0]**2+y[1]**2
    dvdt = -mu/r2**1.5*r
    return [drdt[0], drdt[1], drdt[2], dvdt[0], dvdt[1], dvdt[2]]

#Propiedades del sistema
m1 = 1
m2 = 0.5
M = m1+m2
MU = 1.0

#condiciones iniciales
y=[5.0,0.0,0.0,-0.1,0.42,0.0]
t=np.linspace(0,100.0,1000)
params=dict(mu=MU)
solution=odeint(twob,y,t,args=(params,))

#BODIES
r=solution[:,:3]
r1=-m2/M*r
r2=m1/M*r

#plt.show()


cont = np.linspace(0,10,len(r1[:,0]))

for i in range(0,len(r1[:,0]),10):
#while i < len(r1[:,0]):

    plt.plot(r2[:,0],r2[:,1],'--b')
    plt.plot(r1[:,0],r1[:,1],'--r')

    plt.plot(r1[i,0],r1[i,1],'.r',markersize=5,label='Body 1')
    plt.plot(r2[i,0],r2[i,1],'.b',markersize=10,label='Body 2')
    plt.grid()
    plt.legend()
    plt.xlim(xmin=np.min(r2[:,0]-1), xmax=np.max(r2[:,0])+1)
    plt.ylim(ymin=np.min(r2[:,1]-1), ymax=np.max(r2[:,1])+1)
    plt.savefig("xvsy%.4f.png"%(cont[i]))
    plt.close()
    print("van %d de %.1f"%(i,len(r1[:,0])))
    #i = i + 10

os.system("convert *.png two_body.gif")
os.system("rm *.png")
