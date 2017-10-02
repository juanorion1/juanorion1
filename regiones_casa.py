#!/usr/bin/python

#Este programa organizara las coordenadas de AR y DEC para que CASA las organice como regiones.

import numpy as np


Ah, Am, As, Dh, Dm, Ds = np.loadtxt("DataCluster.dat",usecols=(2,3,4,5,6,7), unpack=True)

dataout = open("regiones_less100.crtf","w")
dataout.write("#CRTFv0 CASA Region Text Format version 0\n")

for i in range(len(Ah)):
    dataout.write('ellipse [[%d:%d:%.3f, -0%d.%d.%d.0000], [68.6922arcsec, 168.3284arcsec], 0.00000000deg] coord=J2000, corr=[I], linewidth=4, linestyle=-, symsize=1, symthick=1, color=magenta, font="DejaVu Sans", fontsize=11, fontstyle=normal, usetex=false\n'%(Ah[i],Am[i],As[i],-1*Dh[i],Dm[i],Ds[i]))

dataout.close()
