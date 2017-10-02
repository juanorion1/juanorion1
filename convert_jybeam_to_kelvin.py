#!/usr/bin/python

#Este programa convertira Jy/beam  a kelvin
#Es para cuando tengo un interferometro


import numpy as np
from astropy import units as u

def convert_jy_kelvin(bmajor,bminor,frequ):

    beam_major = bmajor*u.arcsec
    beam_minor = bminor*u.arcsec
    freq       = frequ*u.Hz
    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    beam_area = 2.*np.pi*(beam_major*beam_minor*fwhm_to_sigma**2) 
    equiv = u.brightness_temperature(beam_area, freq)
    T = (u.Jy.to(u.K, equivalencies=equiv))
    return T
