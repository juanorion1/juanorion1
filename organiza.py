#!/usr/bin/python

#Librerias
import numpy as np


#Tablas entrada
year,mass,AH,AM,AS,DH,DM,DS = np.loadtxt("LMC_less100.dat",usecols=(0,1,2,3,4,5,6,7),unpack = True)
name = np.genfromtxt("LMC_less100.dat",usecols=(8),unpack = True,dtype=str)
name2 = np.genfromtxt("magnitude.dat",usecols=(1),unpack=True,dtype=str)
M_V = np.genfromtxt("magnitude.dat",usecols=(0),unpack=True)
#data3 = np.loadtxt("magnitude.dat",usecols=(0),dtype=float)



#Tablas salida
dataout = open("DataCluster.dat","w")
dataout.write("#l(y)\tl(m)\tAH\tAM\tAS\tDH\tDM\tDS\tM_V\tName\n")

#Ciclo
for i in range(len(year)):
    for j in range(len(M_V)):
        if name[i] == name2[j]:
            
            dataout.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t'%s'\n"%(year[i],mass[i],AH[i],AM[i],AS[i],DH[i],DM[i],DS[i],M_V[j],name[i]))


dataout.close()
