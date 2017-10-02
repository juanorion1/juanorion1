#!/usr/bin/python

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord,Distance

#Arreglo de distancias en arcsec y distancia a LMC
deltaR = np.array([0,5,10,20,30,60])
distLMC  = Distance(distmod=18.50,unit=u.pc) #el modulo de la distancia es de Baumgardt

#Nombre de las tablas
tablas = []
for i in range(len(deltaR)-1):
    tablas.append("distancia_%i-%i_arcsec.dat"%(deltaR[i],deltaR[i+1]))

#Datos paper Bitsakis, T.
ra1,dec1, age,lowage,upage = np.loadtxt("table2.txt",usecols=(1,2,4,5,6),unpack=True)
Name1 = np.genfromtxt("table2.txt",usecols=(0),unpack=True,dtype=str)

#Datos paper Baumgardt
arh,arm,ars,deg,dem,des,age2 = np.loadtxt("DataCluster.dat",usecols=(2,3,4,5,6,7,0),unpack=True)
Name2 = np.genfromtxt("DataCluster.dat",usecols=(8),unpack=True,dtype=str)

#Conversion coordenadas paper Baumgardt
coord = SkyCoord(ra=(arh,arm,ars),dec=(deg,dem,des),unit='deg')
ra2 = coord.ra.value*15
dec2 = coord.dec.value

#Ciclos para determinar la cantidad de cumulos 
for k in range(len(deltaR)-1):
    dataout = open(tablas[k],"w")
    dataout.write("#Las columnas 1-3 son para el catalogo de Baumgardt\n")
    dataout.write("#Las columnas 4-6 para el catalogo de Bitsakis, T.\n\n")
    dataout.write("%10s %10s %10s %10s %10s %10s %18s %18s\n"%("#RA","DEC","Nombre","RA","DEC","Nombre","Distancia[pc]","Distancia[arcsec]"))

    for i in range(len(ra2)): #ciclo para el paper de Baumgardt
        for j in range(len(ra1)): #ciclo para el paper Bitsakis, T.
            c1 = []
            c2 = []
            if (deltaR[k]/3600.)< np.sqrt( (ra2[i] - ra1[j])**2 + (dec2[i] - dec1[j])**2) <= deltaR[k+1]/3600.:
                c1 = (SkyCoord(ra2[i]*u.deg,dec2[i]*u.deg,distance=distLMC))
                c2 = (SkyCoord(ra1[j]*u.deg,dec1[j]*u.deg,distance=distLMC))
                dataout.write("%10f %10f %10s %10f %10f %10s %18f %18f\n"%(ra2[i],dec2[i],Name2[i],ra1[j],dec1[j],Name1[j],c1.separation_3d(c2).value,np.sqrt( (ra2[i] - ra1[j])**2 + (dec2[i] - dec1[j])**2)*3600.))
                

    dataout.close()


#Tomando una distancia maxima de 60 arcsec para cada cumulo de Baumgardt
#se hara un conteo de cuantos cumulos del catalogo de Bitsakis, T hay al rededor
dataout2 = open("Conteo_Cumulos.dat","w")

#Calculando la distancia de cada cumulo de Baumgardt con todos los de Bitsakis, T.
cont = []
for i in range(len(ra2)):
    for j in range(len(ra1)):
        cont.append(np.sqrt( (ra2[i] - ra1[j])**2 + (dec2[i] - dec1[j])**2))

#Teniendo las distancias en arcsec para cada cumulo de Baumgardt
#se organiza una variable de manera tal que separe para cada cumulo de Baumgardt
#la distancia a todos los cumulos de Bitsakis, T.
cont = pd.DataFrame.from_dict(cont)
dist = []
dist = np.append(dist,cont)
dist.shape=(len(ra2),dist.shape[0]/len(ra2))
q = np.where(dist < 60/3600.,dist,0)    #Aqui le pongo cero a todas las distancias que superen los 60 arcsec


#Reemplazo por 1 todas las distancias que son diferentes de cero
dummy = np.zeros((q.shape[0],q.shape[1]))
I = []
dataout2.write("%10s %10s %10s %10s\n"%("#RA","DEC","Nombre","Conteo_cumulos"))
for i in range(q.shape[0]):
    for j in range(q.shape[1]):
        if q[i][j] != 0.0:
            dummy[i][j] = 1.0
            I.append(i)

#Se suma para cada cumulo de Baumgardt todas las entradas de Bitsakis, T.
a = np.sum(dummy,axis=1)
for i in range(len(I)):
    dataout2.write("%10f %10f %10s %10i\n"%(ra2[I[i]],dec2[I[i]],Name2[I[i]],a[I[i]]))

dataout2.close()
