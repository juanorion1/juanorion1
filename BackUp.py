#!/usr/bin/python



###################################################################################################################################################
################## Este programa buscara el valor de los pixeles de una imagen y determinara la coordenada de cada pixel. #########################
###################################################################################################################################################

###################################################################################################################################################
###################################################################### Librerias ##################################################################
###################################################################################################################################################

#import pyfits as pf
import time
from astropy.wcs import WCS
from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

start_time = time.time()


###################################################################################################################################################
###################################################################### Constantes #################################################################
###################################################################################################################################################

size_pix  = 8.33e-3                      #Grados
miu       = 18.50                        #Modulo de la distancia considerado en Baumgardt et al 2013
distLMC   = 10**(1+(miu/5.0))            #Unidades en pc
alpha_co  = 6.6                          #Unidades de Msun/(K*km/s*pc**2)
omega     = np.deg2rad(size_pix)**2  #angulo solido del pixel(es un area)
A         = omega*(distLMC*distLMC)      #Unidades de sr**2*pc**2 (A=Superficie de proyeccion sobre la esfera) 


###################################################################################################################################################
############################################################# Cargando columnas e imagenes ########################################################
###################################################################################################################################################


im  = fits.open('LMC_D2.fits')
rms = fits.open('RMS_LMC_D2.fits')
age,mass,arh, arm, ars, deg, dem, des = np.loadtxt("DataCluster.dat",usecols=(0,1,2,3,4,5,6,7),unpack=True)   #AR y DE de los cumulos
name = np.genfromtxt("DataCluster.dat",usecols=(8),unpack=True,dtype=str)

dataout = open("Porcentaje_Masas.dat","w")
dataout.write("#Cx\tCy\tMClus\tAge\tMRegi\t%M\tAR\t\tDE\t\tName\n")


print "listo tablas"

###################################################################################################################################################
##################################################### Cargando los datos de la imagen y el header #################################################
###################################################################################################################################################

scidata    = im[0].data
scirms     = rms[0].data
intensidad = scidata[0].T
intrms     = scirms[0].T
headerrms  = wcs.WCS(rms[0].header)
header     = wcs.WCS(im[0].header)

print "listo header"
###################################################################################################################################################
##################################################################  Arreglos  #####################################################################
###################################################################################################################################################

AR,DE,SN,Name,Lco,Mmol,Age,Mass,ArCum,DeCum,Cx,Cy,Int,IntRMS,y  = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

print "listo arreglos"

###################################################################################################################################################
############################################ Haciendo conversion de coordenadas a pixeles de los cumulos ##########################################
###################################################################################################################################################
for i in range(len(arh)):
    ArCum.append((arh[i]+(arm[i]/60.0)+(ars[i]/3600.0))*15)
    DeCum.append(deg[i]-(dem[i]/60.0)-(des[i]/3600.0))          #Revisar si ambos llevan el menos

for i in range(len(ArCum)):
    y.append(header.all_world2pix(ArCum[i],DeCum[i],0,1))
    q = np.array(y)
    if 0.0<q[i,0]-0.5 <801.0 and 0<q[i,1]-0.5 <736.0:
        Cx.append(q[i,0]-0.5)
        Cy.append(q[i,1]-0.5)
        Age.append(age[i])
        Mass.append(mass[i])
        Name.append(name[i])
        AR.append(ArCum[i])
        DE.append(DeCum[i])

print "listo conversion"



###################################################################################################################################################
############################################ Guardando intensidades de la imagen original y la RMS ################################################
###################################################################################################################################################

for k in range(len(Cx)):
    for i in range(-1,2):	
        for j in range(-1,2):
            if 0 < Cx[k] < 801.0 and 0 < Cy[k] < 736.0:
                Int.append(intensidad[int(i+Cx[k]),int(j+Cy[k])])
                IntRMS.append(intrms[int(i+Cx[k]),int(j+Cy[k])])



###################################################################################################################################################
##############################Cambiando valores de nan por cero en la imagen original, en el RMS cambiandolos por uno #############################
###################################################################################################################################################

dfint = pd.DataFrame.from_dict(Int)                                    #DataFrame intensidad
dfintrms = pd.DataFrame.from_dict(IntRMS)                              #DataFrame intensidad RMS
A_Int=np.array([])
A_IRMS = np.array([])
A_Int = np.append(A_Int,dfint.replace('NaN',0.0))
A_IRMS = np.append(A_IRMS,dfintrms.replace('NaN',1.0))
A_Int.shape=(A_Int.shape[0]/9.,9.)
A_IRMS.shape=(A_IRMS.shape[0]/9.,9.)



###################################################################################################################################################
######################################## Guardando el SN, calculando el valor de la luminosidad y la masa #########################################
###################################################################################################################################################


for i in range(len(A_Int)):
    for j in range(len(A_Int[0])):
        SN.append(A_Int[i][j]/A_IRMS[i][j])


sn = np.array(SN)
sn.shape=(sn.shape[0]/9.,9.)
sn = np.where(sn > 3,sn,0)
INT = sn*A_IRMS
Lco = INT*A
Mmol= Lco*alpha_co
L = np.sum(Lco,axis=1)
M = np.sum(Mmol,axis=1)



###################################################################################################################################################
################################################ Finalizando, se guarda la tabla con el calculo completo ##########################################
###################################################################################################################################################

m=[]
for i in range(len(Mass)):
    m.append(10**Mass[i])

for i in range(len(M)):
    if M[i] > 0:
        dataout.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%s\n"%(int(Cx[i]),int(Cy[i]),Mass[i],Age[i],np.log10(M[i]),(M[i]/10**Mass[i])*100,AR[i],DE[i],Name[i]))


print("---%s seconds---"%((time.time()-start_time)))



dataout.close()
