#!/usr/bin/python

#Este programa buscara que cumulos estan en 2 catalogos diferentes
#tambien mirara que tan diferentes tienen los valores de masa y edad

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
import scipy.optimize as opt
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord,Distance
import os

os.system("rm -rf distancia_*.dat")
os.system("rm -rf Conteo_Cumulos.dat")

#Arreglo de distancias en arcsec y distancia a LMC                                                                                                  
 
deltaR = 30                                  #en segundos de arco
distLMC  = Distance(distmod=18.50,unit=u.pc) #el modulo de la distancia es de Baumgar

#Nombre de las tablas                                                                                                                               
 
tablas = []
dataout  = open("distancia_%d_arcsec.dat"%(deltaR),"w")
dataout2 = open("Cumulos_Nuevos.dat","w")

#Steps
Step = np.array(['busqueda','conteo'])

step = Step[0]

if step == 'busqueda':

    #Datos paper Popescu
    
    arh,arm,ars,deg,dem,des,L_age,age,U_age,L_mass,mass,U_mass= np.loadtxt("de_Grijs_LMC.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10,11),unpack=True)

    Name = np.genfromtxt("de_Grijs_LMC.dat",usecols=(12),unpack=True,dtype=str)

    #Datos paper Baumgardt                                                                                                                       
 
    arh2,arm2,ars2,deg2,dem2,des2,age2,e_age2,mass2 = np.loadtxt("LMC.dat",usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    Name2 = np.genfromtxt("LMC.dat",usecols=(9),unpack=True,dtype=str)

    #Conversion coordenadas paper Baumgardt                                                                                                      

    coord = SkyCoord(ra=(arh2,arm2,ars2),dec=(deg2,dem2,des2),unit='deg')
    ra2 = coord.ra.value*15
    dec2 = coord.dec.value

    #Conversion coordenadas paper de Grijs                                                                                                      

    coord1 = SkyCoord(ra=(arh,arm,ars),dec=(deg,dem,des),unit='deg')
    ra1 = coord1.ra.value*15
    dec1 = coord1.dec.value

    #Ciclos para determinar la cantidad de cumulos que hay entre los catalogo
    
    dataout.write("#Las columnas 1-3 son para el catalogo de Baumgardt\n")
    dataout.write("#Las columnas 4-6 para el catalogo de Popescu\n\n")
    dataout.write("%10s %10s %10s %10s %10s %22s %18s %18s\n"%("#RA","DEC","Nombre","RA","DEC","Nombre","Distancia[pc]","Distancia[arcsec]"))

    #Asegurando que solo queden los cumulos con masas mayores a 10 mil masas solares
    B = []
    G = []
    for i in range(len(ra2)):
        if mass2[i] >= 4.0:
            B.append(i)
    for i in range(len(ra1)):
        if mass[i] >= 4.0:
            G.append(i)

    I = []
    J = []
    for i in range(len(B)): #ciclo para el paper de Baumgardt                                                                                   
             for j in range(len(G)): #ciclo para el paper de Grijs.                                                                              
                    c1 = []
                    c2 = []
                    if np.sqrt( (ra2[B[i]] - ra1[G[j]])**2 + (dec2[B[i]] - dec1[G[j]])**2) < deltaR/3600.:
                        I.append(B[i])
                        J.append(G[j])
                        c1 = (SkyCoord(ra2[B[i]]*u.deg,dec2[B[i]]*u.deg,distance=distLMC))
                        c2 = (SkyCoord(ra1[G[j]]*u.deg,dec1[G[j]]*u.deg,distance=distLMC))
                    
                        dataout.write("%10f %10f %10s %10f %10f %22s %18f %18f\n"%(ra2[B[i]],dec2[B[i]],Name2[B[i]],ra1[G[j]],dec1[G[j]],Name[G[j]],c1.separation_3d(c2).value,np.sqrt( (ra2[B[i]] - ra1[G[j]])**2 + (dec2[B[i]] - dec1[G[j]])**2)*3600.))

    #Obtengo los cumulos que no cumplen con la condicion                                                                                         
    #de estar en ambas listas pero que si cumplen la condicion de                                                                                 
    #tener masas mayores a 10 mil masas solares
    dataout2.write("%10s %10s %10s %10s %10s %10s %10s %10s %20s\n"%("#rah","ram","ras","deg","dem","des","Edad","Masa","Nombre"))
    id_cluster = list(set(G) - set(J))
    for i in range(len(id_cluster)):
        if 8.3 <= age[id_cluster[i]] <= 9.0:
            dataout2.write("%10d %10d %10f %10d %10d %10f %10f %10f %20s\n"%(arh[id_cluster[i]],arm[id_cluster[i]],ars[id_cluster[i]],deg[id_cluster[i]],dem[id_cluster[i]],des[id_cluster[i]],age[id_cluster[i]],mass[id_cluster[i]],Name[id_cluster[i]]))
                                                                                                      




    #Realizando el ajuste de los cumulos que son iguales en ambos catalogos                                                     
    #Tambien se grafica este ajuste en el la grafica de masas y se compara                                                                       
    #con la linea de masas iguales, tambien se obtiene el coeficiente de correlacion de Pearson
    x = mass2[I]
    y = mass[J]
    X = np.linspace(4,np.max(mass2[I]),len(x))
    Y = np.linspace(4,np.max(mass[J]),len(x))
    def funcion(x,*p): return p[0]*x + p[1]
    param, pcor = opt.curve_fit(funcion,x,y,p0=(1,0))
    
    a,b = param[0],param[1]
    y0 = a*x + b
    diag = np.diag(pcor)
    r = ss.pearsonr(x,y)
    plt.plot(x,y,'.k')
    plt.plot([],[],'.w', label = "r = %.3f"%(r[0]))
    plt.plot(x,y0,'--r',label='m = %.3f+/-%.3f'%(a,np.sqrt(diag[0])))
    plt.plot(X,Y,'-b',label='m = 1')
    plt.legend(loc='upper left',frameon=False)
    plt.xlabel("Log(mass) Baumgardt")
    plt.ylabel("Log(mass) de Grijs")
    plt.savefig("Mass_Comparation_Baum_Grijs.png")
    plt.close()
    print("Finaliza ciclo para determinar la cantidad de cumulos que hay entre catalogos\n")

    dataout.close()
    dataout2.close()
#step = Step[1]
if step == 'conteo':
    #Tomando una distancia maxima de 60 arcsec para cada cumulo de Baumgardt                                                                              
    #se hara un conteo de cuantos cumulos del catalogo de Bitsakis, T hay al rededor                                                                      
    dataout2 = open("Conteo_Cumulos.dat","w")

    #Calculando la distancia de cada cumulo de Baumgardt con todos los de Bitsakis, T.                                                                    
    cont = []
    RA2,DEC2,NAME2 = [],[],[]

    for i in range(len(B)):
        for j in range(len(G)):

            cont.append(np.sqrt( (ra2[B[i]] - ra1[G[j]])**2 + (dec2[B[i]] - dec1[G[j]])**2))
            RA2.append(ra2[B[i]])
            DEC2.append(dec2[B[i]])
            NAME2.append(Name2[B[i]])

    #Teniendo las distancias en arcsec para cada cumulo de Baumgardt
    #se organiza una variable de manera tal que separe para cada cumulo de Baumgardt 
    #la distancia a todos los cumulos de Popescu
            
    cont = pd.DataFrame.from_dict(cont)
    dist = []
    dist = np.append(dist,cont)
    dist.shape=(len(RA2),dist.shape[0]/len(RA2))
    q = np.where(dist < deltaR/3600.,dist,0)    #Aqui le pongo cero a todas las distancias que superen los 20 arcsec

    #Reemplazo por 1 todas las distancias que son diferentes de cero

    dummy = np.zeros((q.shape[0],q.shape[1]))
    d = []
    dataout2.write("%10s %10s %22s %10s\n"%("#RA","DEC","Nombre","Conteo_cumulos"))
    for i in range(q.shape[0]):
        for j in range(q.shape[1]):
            if q[i][j] != 0.0:
                dummy[i][j] = 1.0
                d.append(i)

    #Se suma para cada cumulo de Baumgardt todas las entradas de Bitsakis, T.                                                                             
    a = np.sum(dummy,axis=1)
    
    for i in range(len(d)):
        dataout2.write("%10f %10f %22s %10i\n"%(RA2[d[i]],DEC2[d[i]],NAME2[d[i]],a[d[i]]))

    dataout2.close()
