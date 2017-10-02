#!/usr/bin/python

#Este programa determinara el limite inferior de masa molecular que se puede detectar
#dado el ruido, el tamanio del beam, la frecuencia, etc

import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt

#bmi = np.linspace(1,20,100)
#bma = 1.5*bmi

#noise  = np.linspace(1e-3,1e-1,100)

#Constantes
noise         = 1.920703e-02                                             #Jy/beam km/s
alpha_CO      = 4.35*u.Msun/(u.K*u.km*(u.s**-1)*u.pc*u.pc)               #en Msun*(K kms--1 pc**2)**-1
beam_major    = 50*u.arcsec
beam_minor    = 4*u.arcsec
freq          = 115.2081*u.GHz
fwhm_to_sigma = 1./(8*np.log(2))**0.5

#Calculando el area del beam en arsec**2 y la temperatura
beam_area = 2.*np.pi*(beam_major*beam_minor*fwhm_to_sigma**2)            #en arcsec**2 
print beam_area
equiv     = u.brightness_temperature(beam_area, freq)
T         = u.Jy.to(u.K, equivalencies=equiv)                            #Temperatura
print T

#Conversion del tamanio del beam en pc y calculando el area
beam_maj_pc  = np.deg2rad(beam_major.value/3600.)*1.86e6*u.parsec        #beam en rad * distNGC300
print beam_maj_pc
beam_min_pc  = np.deg2rad(beam_minor.value/3600.)*1.86e6*u.parsec        #beam en rad * distNGC300
print beam_min_pc
area_beam_pc = 2.*np.pi*(beam_maj_pc*beam_min_pc*fwhm_to_sigma**2)       #en pc**2   
print area_beam_pc

#Calculando el limite de deteccion de masa, en Msun
Noise        = (noise*T)*u.K*u.km*(u.s**-1)                              #K km/s
print Noise
sigma        = Noise*area_beam_pc                                        #K km/s pc**2
print sigma
Det_lim_Msun = 5*sigma*alpha_CO
print Det_lim_Msun

#print "Detection limit=%.2e Msun"%(Det_lim_Msun.value)
"""
plt.plot(bmi,Det_lim_Msun.value,'-b',label="beam_min")
plt.plot(bma,Det_lim_Msun.value,'-k',label="beam_max")
plt.xlabel("beam size[arcsec]")
plt.ylabel("Detection_limit_Msun [Msun]")
plt.legend()
plt.title("Noise=%.3f"%(noise))
plt.savefig("beam_min_fix_VS_Det_mass.png")
plt.close()
"""
