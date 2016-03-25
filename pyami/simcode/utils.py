#!/usr/bin/env python
# original code by Thatte, Anand, Sahlmann, Greenbaum
# utility routines and constants for AMI image simulations reoganized by Sahlmann, Anand 3/2016
# anand@stsci.edu 18 Mar 2016

"""
"""

import numpy as np

cdsreadnoise = 21.0                      # CDS read noise (e-)
readnoise = cdsreadnoise/np.sqrt(2)      # read noise for one frame
darkcurrent = 0.012                      # 0.012 e-/sec 
background = 0.125                       # 0.125 e-/sec 
ips_size = 256                           # holdover from before AMISUB became 80x80
flat_sigma = 0.001                       # flat field error
pixscl = 0.065                           # arcsec/pixel WebbPSF 0.064 - DL 0.065
tframe = 0.0745                          # frame time for NISRAPID on AMI SUB80
amisubfov = 80
dither_stddev_as = 0.015                 # 15 mas placement error one-axis
jitter_stddev_as = 0.007                 # 7 mas level 2 reqt on JWST, arcsec

# Anand's email 2016-02-10
F277W, F380M, F430M, F480M = ("F277W", "F380M", "F430M", "F480M")
ZP = {F277W: 26.14,  
      F380M: 23.75,
      F430M: 23.32,
      F480M: 23.19}




# fast rebin Klaus Pontooppidan found on the web
def krebin(a, shape):
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).sum(-1).sum(1)
# legacy slow rebin rewritten to use fast rebin
def rebin(a = None, rc=(2,2), verbose=None):
        r, c = rc
        R, C = a.shape
        sh = (int(R//r), int(C//c))
        return krebin(a, sh)


def jitter(no_of_jitters, osample):
    """ returns in oversampled pixel units.  
        no_of_jitters is known as nint in STScI terminology
    """ 
    mean_j, sigma_j = 0, jitter_stddev_as * osample / pixscl

    xjit = np.random.normal(mean_j,sigma_j,no_of_jitters)
    xjit_r = [int(round(n, 0)) for n in xjit]

    yjit = np.random.normal(mean_j,sigma_j,no_of_jitters)
    yjit_r = [int(round(n, 0)) for n in yjit]
    return xjit_r, yjit_r


def create_ramp(array, _fov, ngroups, frametime):
    nreadouts = ngroups + 1
    readnoise_cube = np.zeros((nreadouts,int(_fov),int(_fov)), np.float64)
    poisson_noise_cube = np.zeros((nreadouts,int(_fov),int(_fov)), np.float64)
    cumulative_poisson_noise_cube = np.zeros((nreadouts,int(_fov),int(_fov)), np.float64)
    ramp=np.zeros((nreadouts,int(_fov),int(_fov)), np.float64)

    for i in range(nreadouts):
        #calculate poisson noise for single reads, then calculate poisson noise for reads up-the-ramp
        if i == 0:
            poisson_noise_cube[i,:,:] = 0.0
        else:
            poisson_noise_cube[i,:,:] = np.random.poisson(array)
        cumulative_poisson_noise_cube = np.cumsum(poisson_noise_cube, axis=0)
        readnoise_array = np.random.normal(0, readnoise, (int(_fov),int(_fov)))
        readnoise_cube[i,:,:] = readnoise_array
        ramp = cumulative_poisson_noise_cube + readnoise_cube + darkcurrent*frametime*i + background*frametime*i
    return ramp


def create_exposure(utr, ngroups, fov):
    xval = np.zeros((ngroups+1,int(fov),int(fov)), np.float64)
    slope = np.zeros((int(fov),int(fov)), np.float64)
    for i in range(ngroups+1):
        xval[i]=i
    xm=float(ngroups)/2.0
    slope = (np.sum(xval*utr,axis=0)-xm*np.sum(utr,axis=0))/(np.sum(xval**2,axis=0)-ngroups*xm**2)
    return slope


#origin is at bottom left of the image. (ds9?)
def apply_padding_image(a,e_x, e_y, fov, osample):

    err_x = int(e_x)
    err_y = int(e_y)

    if err_x <=  0 and err_y <=  0:
       b = np.pad(a, [(0,abs(err_y)),(0,abs(err_x))],mode='constant')
       c = b[abs(err_y):,abs(err_x):]

    elif err_x >=  0 and err_y <=  0:
       b = np.pad(a, [(0,abs(err_y)),(abs(err_x),0)],mode='constant')
       c = b[abs(err_y):,:(fov*osample)]

    elif err_x <=  0 and err_y >=  0:
       b = np.pad(a, [(abs(err_y),0),(0, abs(err_x))],mode='constant')
       c = b[:(fov*osample),abs(err_x):] 

    elif err_x >=  0 and err_y >=  0:
       b = np.pad(a, [(abs(err_y),0),(abs(err_x), 0)],mode='constant')
       c = b[:(fov*osample),:(fov*osample)]

    return c


#padding of 1 for IPS to avoid division by 0 when divided by IPS flat.
def apply_padding_ips(a,e_x,e_y, fov, osample):

    err_x, err_y = (int(e_x), int(e_y))

    if err_x <=  0 and err_y <=  0:
       b = np.pad(a, [(0,abs(err_y)),(0,abs(err_x))],mode='constant',constant_values=(1,1))
       c = b[abs(err_y):,abs(err_x):]
    elif err_x >=  0 and err_y <=  0:
       b = np.pad(a, [(0,abs(err_y)),(abs(err_x),0)],mode='constant',constant_values=(1,1))
       c = b[abs(err_y):,:(fov*osample)]
    elif err_x <=  0 and err_y >=  0:
       b = np.pad(a, [(abs(err_y),0),(0, abs(err_x))],mode='constant',constant_values=(1,1))
       c = b[:(fov*osample),abs(err_x):] 
    elif err_x >=  0 and err_y >=  0:
       b = np.pad(a, [(abs(err_y),0),(abs(err_x), 0)],mode='constant',constant_values=(1,1))
       c = b[:(fov*osample),:(fov*osample)]
    return c
