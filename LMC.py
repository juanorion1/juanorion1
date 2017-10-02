#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


#Datos iniciales
ahh,amm,ass,dhh,dmm,dss,lyrBest,lyrerror,lMBest  = np.loadtxt("LMC.dat",usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
Name = np.genfromtxt("LMC.dat",usecols=(9),dtype=str)

#Creando una nueva tabla
#dataout = open("LMCF.dat","w")       #Tabla con los datos finales
dataout2 = open("LMC_less100.dat","w")
#dataout.write("#l(y)\tl(m)\tAH\tAM\tAS\tDH\tDM\tDS\tName\n")
dataout2.write("#l(y)\tl(m)\tAH\tAM\tAS\tDH\tDM\tDS\tName\n")

#Arreglos
N   = []                      #Nombre
AG  = []                      #Log de la edad.

"""
for i in range(len(Name)):
        if lyrBest[i] >= 8.3 and lyrBest[i] <= 9.0 and lMBest[i] >= 4:
            
            dataout.write("%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%s\n"%(lyrBest[i],lMBest[i],ahh[i],amm[i],ass[i],dhh[i],dmm[i],dss[i],Name[i]))
"""

for i in range(len(Name)):
        if lyrBest[i] < 8.3 and lMBest[i] >= 4:
            dataout2.write("%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%s\n"%(lyrBest[i],lMBest[i],ahh[i],amm[i],ass[i],dhh[i],dmm[i],dss[i],Name[i]))

#dataout.close()
dataout2.close()
