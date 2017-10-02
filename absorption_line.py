#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


def func(x, a, mu, sigma):
    return a*np.exp( (-(x-mu)**2) / (2*sigma**2) )

NOISE = 1/20.

# Generating clean data
x = np.linspace(220, 339 ,266 )
y1 = func(x, 0.225, 257.69, 1.93)

# Adding noise to the data
#yn = y1 + np.random.normal(scale=NOISE,size=len(x))

#param, pcov = opt.curve_fit(func,x,yn)

#print param
#print "std dev amplitud=%f"%np.sqrt(pcov[0][0])
#print "std dev mu=%f"%pcov[1][1]
#print "std dev sigma=%f"%pcov[2][2]


plt.plot(x,y1,'b-')#,label='gaussiana+error')
#plt.plot(x,func(x,param[0],param[1],param[2]),'b-',label='Ajuste')
#plt.legend(loc='best')
plt.show()
#plt.savefig("absorption_line_fit.png")



"""
N = 100

matriz = []
val_a = []
val_mu = []
val_sig = []
sigma_a  = []
sigma_mu  = []
sigma_sig  = []

while len(val_a) < N:
    yn = (y1 + np.random.normal(scale=NOISE,size=len(x)))

    try:
        param, pcov = opt.curve_fit(func,x,yn)
        val_a.append(param[0])
        val_mu.append(param[1])
        val_sig.append(param[2])
        matriz.append(pcov)
        
    except Exception:
        print "Buscando convergencia"
   
for i in range(N):
    sigma_a.append(matriz[i][0][0])
    sigma_mu.append(matriz[i][1][1])
    sigma_sig.append(matriz[i][2][2])
"""    



