#!/usr/bin/env python
#Code to simulate PSF of binary star pair with NIRISS NRM.  
#Jan 19 2016 Thatte mac os x
#Jan 21 linux  fixed compatibility issues re np.pad integer args required  Anand
#Feb 10 2016 version emailed from Anand to Johannes


"""
This code is an updated version of make_binarypair_sim_jan14.py(January 14, 2014 version)
Updates - 4 dithers instead of 9 with updated pointings 
          TFRAME = 0.0745 sec instead of 2.688 sec, fov = 80 instead of 128
          Using 'nint' instead of exposures, using 'ngroups' instead of non_destructive_reads
          x and y are still arbitray similar to the earlier version
          Uses astropy.io.fits instead of pyfits
		  JS+AZG+AS: magA fluxratio sep pa filt ngroups nint Feb 10
"""

import sys
import os
import numpy as np
from astropy.io import fits 
from astropy.io.fits import getheader
import webbpsf as wp
import pysynphot
from poppy import specFromSpectralType
import time

import pyami.simcode.utils as U   
 
osample = 11                             # hardcoded 11 in this file replaced by this global variable


# CREATE NOISELESS POINT SOURCE PSF USING WEBBPSF
def generate_starPSF(FILTER=None, fov=None, osample=None, spectraltype="A0V"):         
    niriss = wp.NIRISS()
    niriss.filter = FILTER
    niriss.pupil_mask = 'MASK_NRM'
    niriss.pupilopd = ("OPD_RevV_niriss_162.fits", 3)
    niriss.pixelscale=0.065
    src = specFromSpectralType(spectraltype)

    
    #Create an oversized array for star PSF. 
    psf_fits = niriss.calcPSF(fov_pixels=fov + 4,oversample=osample,source=src,rebin=False,clobber=True)
    psf_array = psf_fits[0].data
    psf_header = psf_fits[0].header
    print psf_array.sum(), 'sum of star psf'  
    return psf_array, psf_header


def simulate_skydata(_trials, _binarystar_array, _cubename, _dithers,  _x_dith, _y_dith, ngroups, nint, frametime, filt='F430M', outDir = '', tmpDir = '', **kwargs):
    
    offset_x = kwargs['offset_x']
    offset_y = kwargs['offset_y']
    fluxratio = kwargs['fluxratio']
    fov = kwargs['fov']
    starmag = kwargs['starmag']
    magzeropoint = kwargs['magzeropoint']
    sptype = kwargs['sptype']
    osample = kwargs['osample']
    dim = kwargs['dim']
    
    flux = 10**(-(starmag-magzeropoint)/2.5)

    cube = np.zeros((nint,int(fov),int(fov)), np.float64)
    ipsarray_big = np.ones((U.ips_size*osample,U.ips_size*osample))

    #Create target data
    for p in range(_trials):
         print 'Starting trial', p
         #CALCULATE LOCATIONS OF 4 DITHERS WITH 15 MAS ERROR ON 256 X 11 ARRAY
         mean_d, sigma_d = 0, 15.0 * osample/65.0
         x_dith_error = np.random.normal(mean_d,sigma_d, _dithers)
         x_dith_error_r = [int(round(n, 0)) for n in x_dith_error]
         #Accumulate dither errors
         x_dith_error_accum = np.cumsum(x_dith_error_r)
         dither_xcenter = [a + b for a, b in zip(_x_dith, x_dith_error_accum)] 
    
         y_dith_error = np.random.normal(mean_d,sigma_d, _dithers)
         y_dith_error_r = [int(round(n, 0)) for n in y_dith_error]
         #Accumulate dither errors
         y_dith_error_accum = np.cumsum(y_dith_error_r)
         dither_ycenter = [a + b for a, b in zip(_y_dith, y_dith_error_accum)]
    
         print "\tPrinting commanded dither, accumulated error, final dither location for verification"
         print "  "
         print "\tcommanded X dither", _x_dith
         print "\tAccumulated X dither error", x_dith_error_accum
         print "\tdither_xcenter", dither_xcenter      
         print "  "
         print "\tcommanded Y dither", _y_dith
         print "\tAccumulated Y dither error", y_dith_error_accum
         print "\tdither_ycenter", dither_ycenter
         print "  "
    
         #POSITIONAL ERROR = DITHER + JITTER
    
         xjitter = range( _dithers)   #each of the 4 elements is an array of NINT jitters
         yjitter = range( _dithers)   #one set per dither location
         for i in range( _dithers):
             xjitter[i], yjitter[i] = U.jitter(nint, osample)
             print '\t\tx jitter', xjitter[i]
             print '\t\ty jitter', yjitter[i]
         xjitter_array = np.array(xjitter)
         x = range( _dithers)
         yjitter_array = np.array(yjitter)      
         y = range( _dithers)
         total_pos_error_x = range( _dithers)
         total_pos_error_y = range( _dithers)
         for i in range( _dithers):
             x[i] = dither_xcenter[i] + xjitter_array[i]
             y[i] = dither_ycenter[i] + yjitter_array[i] 
             total_pos_error_x[i] = x[i] - _x_dith[i]
             total_pos_error_y[i] = y[i] - _y_dith[i]
             print " "
             print '\ttotal positional error in X', total_pos_error_x[i]
             print '\treal X pointing with dither and jitter', x[i]
             print " "
             print '\ttotal positional error in Y', total_pos_error_y[i]
             print '\treal Y pointing with dither and jitter', y[i] 
             print " "
         
             for k,(ii,jj) in enumerate(zip(x[i],y[i])):
                 print '\ttrial',p ,'dither', i, 'dither location', ii, jj, 'exposure', k
    
                 #ips section for even FOV
                 print "ipsarray_big.shape", ipsarray_big.shape
                 ips_section = ipsarray_big[ii-dim*osample:ii+dim*osample,jj-dim*osample:jj+dim*osample]
    
                 binarystar_ips_array = ips_section * _binarystar_array
                 print '\tjitter error for the exposure - x,y', ii-dither_xcenter[i],jj-dither_ycenter[i]
                 binarystar_ips_array_sh = U.apply_padding_image(binarystar_ips_array,jj-dither_ycenter[i],ii-dither_xcenter[i], fov, osample)
                 fits.writeto(tmpDir+'binarystar_ips_array_sh'+str(k)+'.fits',binarystar_ips_array_sh, clobber = True)
                 print "\tinfo",dither_xcenter[i]-(dither_xcenter[i]//osample)*float(osample),osample-(dither_xcenter[i]-(dither_xcenter[i]//osample)*osample)
                 print "\tinfo",dither_ycenter[i]-(dither_ycenter[i]//osample)*float(osample),osample-(dither_ycenter[i]-(dither_ycenter[i]//osample)*osample) 
                 im = np.pad(binarystar_ips_array_sh, [(int(dither_xcenter[i]-(dither_xcenter[i]//osample)*float(osample)), int(osample-(dither_xcenter[i]-(dither_xcenter[i]//osample)*osample))),
                                                       (int(dither_ycenter[i]-(dither_ycenter[i]//osample)*float(osample)), int(osample-(dither_ycenter[i]-(dither_ycenter[i]//osample)*osample)))],mode='constant')
                 print '\tim.shape is', im.shape
                 rebinned_array_81 = U.rebin(im, (osample,osample))
                 rebinned_array = rebinned_array_81[0:80,0:80]
                 
                 ips_section_sh = U.apply_padding_ips(ips_section,jj-dither_ycenter[i],ii-dither_xcenter[i], fov, osample)
                 fits.writeto(tmpDir+'ips_section_sh'+str(k)+'.fits',ips_section_sh, clobber = True)
                
                 im_ips=np.pad(ips_section_sh, [(int(dither_xcenter[i]-(dither_xcenter[i]//osample)*float(osample)), int(osample-(dither_xcenter[i]-(dither_xcenter[i]//osample)*osample))),
                                                (int(dither_ycenter[i]-(dither_ycenter[i]//osample)*float(osample)), int(osample-(dither_ycenter[i]-(dither_ycenter[i]//osample)*osample)))],mode='constant',constant_values=(1,1))             
                 fits.writeto(tmpDir+'im_ips'+str(k)+'.fits',im_ips, clobber = True)
                 rebinned_ips_flat_81 = U.rebin(im_ips, (osample,osample))/osample**2
                 rebinned_ips_flat = rebinned_ips_flat_81[0:80,0:80]
                 print '\trebinned_ips_flat sum',rebinned_ips_flat.sum()
                 rebinned_array/=rebinned_ips_flat
                 counts_array = flux * rebinned_array * frametime  
                 ramp = U.create_ramp(counts_array, fov, ngroups, frametime)

                 fits.writeto(tmpDir+'ramp.fits',ramp, clobber = True)
                 pflaterror=np.random.normal(0.0, U.flat_sigma, size=(fov,fov))
                 pflat= 1.0 + pflaterror
                 exposure = U.create_exposure(ramp, ngroups, fov)
                 exposure1 = (exposure-U.darkcurrent-U.background)*pflat
                 cube[k,:,:] = exposure1
                 print '\tmax pixel counts', cube[k,:,:].max()
                 print " "                   
    
             outfile = _cubename+str(p)+str(i)+".fits"
             print '\tcreated', _cubename+str(p)+str(i)+'.fits'
             (year, month, day, hour, minute, second, weekday, DOY, DST) =  time.gmtime()
             fitsobj = fits.HDUList()
             hdu = fits.PrimaryHDU(  )
             printhdr = hdu.header
      
             printhdr['INSTRUME'] =  'NIRISS'
             printhdr['pixscl'] = U.pixscl, 'Pixel scale (arcsec/pixel)'
             printhdr['NRMNAME'] =  'G7S6', 'Tuthill Anand Beaulieu Lightsey'
             printhdr['starmag'] = starmag,'Star magnitude'
             printhdr['sptype'] = sptype
             printhdr['NRM_X_A1'] =  0.00000, 'X coordinate (m) of NRM sub-aperture 0'          
             printhdr['NRM_Y_A1'] = -2.64000, 'Y coordinate (m) of NRM sub-aperture 0'         
             printhdr['NRM_X_A2'] = -2.28631, 'X coordinate (m) of NRM sub-aperture 1'          
             printhdr['NRM_Y_A2'] =  0.00000, 'Y coordinate (m) of NRM sub-aperture 1'          
             printhdr['NRM_X_A3'] =  2.28631, 'X coordinate (m) of NRM sub-aperture 2'          
             printhdr['NRM_Y_A3'] = -1.32000, 'Y coordinate (m) of NRM sub-aperture 2'          
             printhdr['NRM_X_A4'] = -2.28631, 'X coordinate (m) of NRM sub-aperture 3'          
             printhdr['NRM_Y_A4'] =  1.32000, 'Y coordinate (m) of NRM sub-aperture 3'          
             printhdr['NRM_X_A5'] = -1.14315, 'X coordinate (m) of NRM sub-aperture 4'          
             printhdr['NRM_Y_A5'] =  1.98000, 'Y coordinate (m) of NRM sub-aperture 4'          
             printhdr['NRM_X_A6'] =  2.28631, 'X coordinate (m) of NRM sub-aperture 5'          
             printhdr['NRM_Y_A6'] =  1.32000, 'Y coordinate (m) of NRM sub-aperture 5'          
             printhdr['NRM_X_A7'] =  1.14315, 'X coordinate (m) of NRM sub-aperture 6'          
             printhdr['NRM_Y_A7'] =  1.98000, 'Y coordinate (m) of NRM sub-aperture 6'   
             printhdr['nframe'] = 1,'Readout number of frames'  
             printhdr['ngroup'] = ngroups,'Readout number of groups'  
             printhdr['framtime'] = frametime, 'one(utr=1)/first-to-last(utr=0) (s)'
             printhdr['units'] = 'photoelectrons'
             printhdr['COMP_DX'] = offset_x*U.pixscl/float(osample), 'Companion separation in X, arcsec'
             printhdr['COMP_DY'] = offset_y*U.pixscl/float(osample), 'Companion separation in Y, arcsec'      
             printhdr['COMP_FR'] = fluxratio, 'Companion flux ratio' 
             printhdr['ffe_err'] =  U.flat_sigma*100, '% Flat field error stddev'
             printhdr['jitter'] =   U.jitter_stddev_as*1000, '1-axis jitter stddev mas'
             printhdr['dith_err'] = U.dither_stddev_as*1000, '1-axis dither placement stddev mas'
             printhdr['dithx%d'%i] = _x_dith[i]/osample, 'Commanded X dither (detpix in ipsarray)'
             printhdr['dithy%d'%i] = _y_dith[i]/osample, 'Commanded Y dither (detpix in ipsarray)'
             printhdr['dithx_r%d'%i] = dither_xcenter[i]/float(osample), 'Real X dither (detpix in ipsarray)'
             printhdr['dithy_r%d'%i] = dither_ycenter[i]/float(osample), 'Real Y dither (detpix in ipsarray)'
             printhdr['CODESRC'] = 'make_binary.py', '[thatte anand jsahlmann]@stsci.edu'

             # Append the header from psf_star.fits, likely created by WebbPSF
             wp_header = getheader(outDir +'star_array_fov80_%s.fits'%filt.lower())
             printhdr.extend(wp_header, update=True)

             #Delete and over-write keyword values written by WebbPSF
             del printhdr['PLANE1']
             del printhdr['DET_SAMP']
             printhdr['oversamp']= osample, 'Oversampling factor for MFT'
             printhdr['AUTHOR'] = '%s@%s' % (os.getenv('USER'), os.getenv('HOST')), 'username@host for calculation'
             printhdr['DATE'] = '%4d-%02d-%02dT%02d:%02d:%02d' %  (year, month, day, hour, minute, second), 'Date of calculation'
             hdu.data = cube
             fitsobj.append( hdu )
             fitsobj.writeto(outDir+outfile, clobber = True)
             fitsobj.close()

             print "\nPeak pixel and total e- in each slice:"
             for i in range(cube.shape[0]):
                 print i, " %.1e"%cube[i,:,:].max(), " %.1e"%cube[i,:,:].sum()
             print ""
