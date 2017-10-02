#!/usr/bin/python

###################################################################################################################################################
################## Este programa buscara el valor de los pixeles de una imagen y determinara la coordenada de cada pixel. #########################
###################################################################################################################################################

###################################################################################################################################################
###################################################################### Librerias ##################################################################
###################################################################################################################################################

import time
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u,wcs
from astropy.coordinates import Distance,SkyCoord
import numpy as np
import pandas as pd

start_time = time.time()

###################################################################################################################################################
############################################################# Cargando columnas e imagenes ########################################################
###################################################################################################################################################

im  = fits.open('LMC_D2.fits')
rms = fits.open('RMS_LMC_D2.fits')
age,mass,arh, arm, ars, deg, dem, des = np.loadtxt("DataCluster.dat",usecols=(0,1,2,3,4,5,6,7),unpack=True)   #AR y DE de los cumulos
name = np.genfromtxt("DataCluster.dat",usecols=(9),unpack=True,dtype=str)
coord = SkyCoord(ra=(arh,arm,ars),dec=(deg,dem,des),unit='deg')

dataout = open("Porcentaje_Masas.dat","w")
dataout.write("#Cx\tCy\tMClus\tAge\tMRegi\t%M\tAR\t\tDE\t\tName\n")


print "listo tablas"

###################################################################################################################################################
##################################################### Cargando los datos de la imagen y el header #################################################
###################################################################################################################################################

scidata    = im[0].data
scirms     = rms[0].data
intensidad = (scidata[0].T)*u.K*u.km*(u.s*u.sr)**-1
intrms     = (scirms[0].T)*u.K*u.km*(u.s*u.sr)**-1
headerrms  = wcs.WCS(rms[0].header)
header     = wcs.WCS(im[0].header)

print "listo header"
###################################################################################################################################################
###################################################################### Constantes #################################################################
###################################################################################################################################################

val_head  = header.to_header()                         #Grados
size_pix  = val_head[5]*u.deg                          #Tamanio del pixel, sacado directamente del header
miu       = 18.50                                      #Modulo de la distancia considerado en Baumgardt et al 2013
distLMC   = Distance(distmod=miu,unit=u.pc)            #Unidades en pc
alpha_co  = 6.6*u.Msun*u.s*(u.K*u.km*u.pc*u.pc)**-1    #Unidades de Msun/(K*km/s*pc**2)
omega     = (np.deg2rad(size_pix)**2).to(u.sr)         #angulo solido del pixel(es un area)
A         = omega*(distLMC**2)                         #Unidades de sr*pc**2 (A=Superficie de proyeccion sobre la esfera) 

###################################################################################################################################################
##################################################################  Arreglos  #####################################################################
###################################################################################################################################################

AR,DE,SN,Name,Lco,Mmol,Age,Mass,ArCum,DeCum,Cx,Cy,Int,IntRMS = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
print "listo arreglos"

###################################################################################################################################################
############################################ Haciendo conversion de coordenadas a pixeles de los cumulos ##########################################
###################################################################################################################################################
q = np.array(header.all_world2pix(coord.ra*15,coord.dec,0,1))
for i in range(len(coord.dec)):
    if 0.0<q[0,i]-0.5 <im[0].header[3] and 0<q[1,i]-0.5 <im[0].header[4]:
        Cx.append(q[0,i]-0.5)
        Cy.append(q[1,i]-0.5)
        Age.append(age[i])
        Mass.append(mass[i])
        Name.append(name[i])
        AR.append(coord.ra[i].value*15)
        DE.append(coord.dec.value[i])

print "listo conversion"
###################################################################################################################################################
############################################ Guardando intensidades de la imagen original y la RMS ################################################
###################################################################################################################################################

for k in range(len(Cx)):
    for i in range(-1,2):	
        for j in range(-1,2):
            Int.append(intensidad.value[int(i+Cx[k]),int(j+Cy[k])])
            IntRMS.append(intrms.value[int(i+Cx[k]),int(j+Cy[k])])
          
###################################################################################################################################################
##############################Cambiando valores de nan por cero en la imagen original, en el RMS cambiandolos por uno #############################
###################################################################################################################################################

dfint = pd.DataFrame.from_dict(Int)                                    #DataFrame intensidad
dfintrms = pd.DataFrame.from_dict(IntRMS)                              #DataFrame intensidad RMS
A_Int  = []
A_IRMS = []
A_Int = np.append(A_Int,dfint.replace('nan',0.0))
A_IRMS = np.append(A_IRMS,dfintrms.replace('nan',1.0))
A_Int.shape=(A_Int.shape[0]/9.,9.)
A_IRMS.shape=(A_IRMS.shape[0]/9.,9.)

###################################################################################################################################################
######################################## Guardando el SN, calculando el valor de la luminosidad y la masa #########################################
###################################################################################################################################################

SN   = A_Int/A_IRMS
SN   = np.where(SN > 3,SN,0)              #Dejando solamente los valores por encima de 3
INT  = SN*A_IRMS
Lco  = A*INT*intensidad.unit
Mmol = Lco*alpha_co
L    = np.sum(Lco,axis=1)
M    = np.sum(Mmol,axis=1)

###################################################################################################################################################
################################################ Finalizando, se guarda la tabla con el calculo completo ##########################################
###################################################################################################################################################

for i in range(len(M)):
    if M[i] > 0:
        dataout.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%s\n"%(int(Cx[i]),int(Cy[i]),Mass[i],Age[i],np.log10(M.value[i]),(M.value[i]/10**Mass[i])*100,AR[i],DE[i],Name[i]))

print("---%s seconds---"%((time.time()-start_time)))

dataout.close()
