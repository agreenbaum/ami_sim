#!/usr/bin/env python
# Code to calculate NRM exposure parameters
# Deepashri Thatte and Anand Sivaramakrishnan
# Double-checked & corrected 24 Nov 2015 for apt values of NINT, NGROUP calculation
# Generation of WebbPSF files was removed.
# Code expects fov 255 and fov 31 detector pixel files to support calculation.  Oversample 3 used initially.

"""
"ETC" code and the required PSFs (they'll unpack on a mac to the right directory name, in the same directory as  the code sits.  You can repoint the code by editing DATADIR's value).  It just tells how long it takes to collect a required number of photons.  It uses the 80x80 subarray, and 35000e saturation for the brightest pixel (it's all in there somewhere).

The "user manual" consists of:

Use the code as follows (default values kick in if you lack enough arguments):

./NIRISSami_apt_calc_2016jan19.py 6.33 1e10
Opening psf file  ./NIRISSami_apt_calcPSF/F277W_31_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F277W_255_A0V_det.fits
NIRISSami_apt_calc_2016jan15.py:207: RuntimeWarning: divide by zero encountered in double_scalars
  nint_apt = nint_ + (ngroups_ - int(ngroups_)) * nint_ / int(ngroups_)
Opening psf file  ./NIRISSami_apt_calcPSF/F380M_31_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F380M_255_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F430M_31_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F430M_255_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F480M_31_A0V_det.fits
Opening psf file  ./NIRISSami_apt_calcPSF/F480M_255_A0V_det.fits

	Using Kevin Volk's zero points, based on V band photometry.  Replace with K=9??
	F277W 26.14
	F380M 23.75
	F430M 23.32
	F480M 23.19
	Target magnitude is 6.33, total  number of detected photons is 1e+10

	F277W: 
		87.5% of nrmflux in 31x31,  0.57% peak frac wrt fullpsf.sum,  3.80% peak frac wrt nrmpsf.sum,  15.04% nrmtot in 255x255  
		CR 8.4e+07  cpfCR 4.801e+05  cptot 4.3e+08  cp_e/frame 3.6e+04  Tsat 7.3e-02, 
		NINT 12410.981825 NINT(updated to account for lost partial group per integration) inf,
		NGROUPS 0.978637 totnframes 12145,
		=====FOR APT INPUT: NGROUPS = 0 NINT = inf=====
		datacollect 15.08110 minutes  
		BRIGHTNESS LIMIT EXCEEDED: some pixels may saturate, use NGROUPS = NINT = 1 in APT if desired.
	F380M: 
		84.8% of nrmflux in 31x31,  0.35% peak frac wrt fullpsf.sum,  2.32% peak frac wrt nrmpsf.sum,  14.95% nrmtot in 255x255  
		CR 9.3e+06  cpfCR 3.220e+04  cptot 2.7e+08  cp_e/frame 2.4e+03  Tsat 1.1e+00, 
		NINT 7807.319092 NINT(updated to account for lost partial group per integration) 8135.786935,
		NGROUPS 14.589005 totnframes 113901,
		=====FOR APT INPUT: NGROUPS = 14 NINT = 8136=====
		datacollect 141.42710 minutes  
	F430M: 
		84.2% of nrmflux in 31x31,  0.28% peak frac wrt fullpsf.sum,  1.88% peak frac wrt nrmpsf.sum,  14.91% nrmtot in 255x255  
		CR 6.3e+06  cpfCR 1.752e+04  cptot 2.2e+08  cp_e/frame 1.3e+03  Tsat 2.0e+00, 
		NINT 6376.864995 NINT(updated to account for lost partial group per integration) 6578.147779,
		NGROUPS 26.820678 totnframes 171031,
		=====FOR APT INPUT: NGROUPS = 26 NINT = 6579=====
		datacollect 212.36454 minutes  
	F480M: 
		82.8% of nrmflux in 31x31,  0.23% peak frac wrt fullpsf.sum,  1.55% peak frac wrt nrmpsf.sum,  14.87% nrmtot in 255x255  
		CR 5.5e+06  cpfCR 1.282e+04  cptot 1.9e+08  cp_e/frame 9.6e+02  Tsat 2.7e+00, 
		NINT 5362.319208 NINT(updated to account for lost partial group per integration) 5458.589561,
		NGROUPS 36.646312 totnframes 196509,
		=====FOR APT INPUT: NGROUPS = 36 NINT = 5459=====
		datacollect 243.99895 minutes  
	Target magnitude is 6.33, total  number of detected photons is 1e+10



	pure synphot, without Kevin's niriss throughput effects:
	anand@ati.st-visitor.org:62  ./NIRISSami_calc_zp.py
	F430M zero point: 23.7854821156
	F480M zero point: 24.1274958381
	F380M zero point: 24.1072792513
	F277W zero point: 26.532110021

	Note Kevin's numbers are:
	F430M: 23.32, # this is ~0.5 mag brighter than the pure synphot calculation

mag lim s centered PSF exact, A0V, 35000e-,   corner psf just 1.5 mag or /= 4 less, using Kevin+Deepashri ZPs
	F277W 6.7 (~5.2 if exactly pixel cornered PSF)
	F380M 3.5 (~2.0 if exactly pixel cornered PSF)
	F430M 3.1 (~1.6 if exactly pixel cornered PSF)
	F480M 2.7 (~1.2 if exactly pixel cornered PSF)


April 1 2015 - 

Hi Anand,
Did you mean to say the following four zeropoints?

In [263]: zp_277w=9+(2.5*np.log10(7161507.0182116))

In [264]: zp_277w
Out[264]: 26.137511054661204

In [265]: zp_380m=9+(2.5*np.log10(794466.0863508))

In [266]: zp_380m
Out[266]: 23.750188407649986

In [267]: zp_480m=9+(2.5*np.log10(473350.0643434))

In [268]: zp_480m
Out[268]: 23.187956101219747

In [269]: zp_430m=9+(2.5*np.log10(535999.4518994))

In [270]: zp_430m
Out[270]: 23.322910863983846


These are using the count rates in the table (attached) that Kevin gave last year.

Thanks
Deepashri

rcp -p witserv.stsci.edu:/witserv/data18/martel/CV1RR/OTP74/

:::::::

anand@imac24:16  
"""
"""
nint_        Exact value of NINT corresponding to exact (float) value of NGROUPS
nint         Updated value of NINT corresponding to integer value of NGROUPS to account for loss of last partial group per integration 
nint_ceil    Integer value of NINT that completes last or only partial integration ----------> THIS IS NINT FOR APT USE
ngroups_     exact (float) value of NGROUPS
ngroups      integer value of NGROUPS (lost partial group gets collected in additional integrations represented by nint) ---------> THIS IS NGROUPS FOR APT USE
"""

import os, sys
import numpy as np
import astropy.io.fits as pyfits

# SRC = "A0V"
# f277w MAG_T = 6.35  for <35ke/TFRAME in cp pixctr   DTfull 11.74  Approx and old.  Numbers in 'live' code are more recent
# f380m MAG_T = 3.42  for <35ke/TFRAME in cp pixctr   DTfull  8.78
# f430m MAG_T = 2.76  for <35ke/TFRAME in cp pixctr   DTfull  8.15
# f480m MAG_T = 2.42  for <35ke/TFRAME in cp pixctr   DTfull  7.80
# MAG_T = 9.0
# TOT_E = 1.0e10 # Greeenbaum et al. ApJ, Ireland MNRAS  contrast (1e4) ~ sqrt(totalphotons)/10
TFRAME = 0.0745 # seconds Alex Fullerton
# SAT_E = 35000.0  # Kevin Volk
infostr = """   old...
	pure synphot, without Kevin's niriss throughput effects:
	anand@ati.st-visitor.org:62  ./NIRISSami_calc_zp.py
	F277W zero point: 26.53
	F380M zero point: 24.11
	F430M zero point: 23.79
	F480M zero point: 24.13

	Kevin Volk's numbers - 
	F277W zero point: 26.18
	F380M zero point: 23.75
	F430M zero point: 23.32
	F480M zero point: 23.19

	Brightness lims pixel-centered PSF exact  (pixel-corner psf approx 1/4 of this,  1.5 mag less)
	F277W 6.7 (~5.2 if exactly pixel cornered PSF)
	F380M 3.5 (~2.0 if exactly pixel cornered PSF)
	F430M 3.1 (~1.6 if exactly pixel cornered PSF)
	F480M 2.7 (~1.2 if exactly pixel cornered PSF)


	2015.10.05  Dean Hines asks for Jy brightness limits:

 Doris Daou's
 http://www.stsci.edu/hst/nicmos/tools/conversion_form.html
 Temperature of the blackbody =    5500.0
 0mag  at L' band center 3.76um  Flux=2.50E+02 Jy
 0mag  at M  band center 4.76um  Flux=1.50E+02 Jy

 so 
 	use L1 mag for 277 380
 	use M  mag for 430 480

	

"""

# F277W, F380M, F430M, F480M = ("F277W", "F380M", "F430M", "F480M")
# """ZP = {F277W: 26.53,  # pure synphot
#       F380M: 24.11,
# 	  F430M: 23.32,
# 	  F480M: 24.13}
# 	  """
# ZP = {F277W: 26.14,  # overwrite with Kevin's
#       F380M: 23.75,
# 	  F430M: 23.32,
# 	  F480M: 23.19}

# OPDFILE = ("OPD_RevV_niriss_162.fits", 3)
# DATADIR = "./NIRISSami_apt_calcPSF/"



def apt_values(ngroups_, nint_):
 	"""
	Corrected May 26 2015 Thatte: 
	report="\n\t\tCR %.1e NINT %.1f NGROUPS_ceil %d NGROUPS_float %f NINT_new %f" %
	(cr,nint_,np.ceil(ngroups_),ngroups_,(nint_+(ngroups_-int(ngroups_))*nint_/int(ngroups_)))

	Compensates for loss of fractional part of float ngroup_

	Sum up how many of the fractional parts of the group are lost over nint_ integrations, 
	then add these back to the original number of integrations (nint_).  Returns integer values for APT use

	"""
	ngroups_apt = int(ngroups_)
	nint_apt = nint_ + (ngroups_ - int(ngroups_)) * nint_ / int(ngroups_)
        nint_apt_ceil = np.ceil(nint_apt)


	return ngroups_apt, nint_apt, nint_apt_ceil

def generatePSF(filt=None, fov=None, osample=5, cr=None, tot_e=None, sat_e=None, SRC=None, return_params=0, DATADIR = "./pyami/etc/NIRISSami_apt_calcPSF/"):

	# 	print os.getcwd()
	if os.access(DATADIR+'%s_%d_%s_det.fits'%(filt,fov,SRC), os.F_OK) == True:

		resfov = pyfits.open(DATADIR+'%s_%d_%s_det.fits'%(filt,fov,SRC))
		resbig = pyfits.open(DATADIR+'%s_%d_%s_det.fits'%(filt,255,SRC))
# 		print "Opening psf file ", DATADIR+'%s_%d_%s_det.fits'%(filt,fov,SRC)
# 		print "Opening psf file ", DATADIR+'%s_%d_%s_det.fits'%(filt,255,SRC)
		nrmfov, hdrfov = (resfov[0].data, resfov[0].header)
		nrmbig, hdrbig = (resbig[0].data, resbig[0].header)
		readfromdisk = True
	else:
		sys.exit("Missing PSF %s directory or PSF files for spectral type/filter %s/%s" % (DATADIR, SRC, filt))

	# add to hdr of small nrm psf
	
	fract_31 = nrmfov.sum()/nrmbig.sum() # fraction of flux in small nrm array compared to "almost infinite" nrm array
	cpf = nrmfov.max() # CPF  compared to CLEAR aperture large psf array total being 1
	nrmtot = nrmbig.sum() # double check NRM throughput is about 15%, clear psf total 1
	cpfnrm = cpf / nrmtot # Central pixel fraction of NRM array

	# total desired electrons in CLEAR given our NRM total electron needs.
	# 1/(nrmtot=0.15)(frac_31=0.85) = Factor of about 8 for F430M
	tot_e_full = tot_e / (nrmtot*fract_31)

	# total number of electrons in central pixel of NRM array given total number in CLEAR psf
	cptot = cpf * tot_e_full

	# NINT = number of ramps
	nint_ = cptot / sat_e

	# given CR, cr*cpf in central pixel: need this many frames to reach saturation
	# NISRAPID NGROUPS = NFRAMES:
	Tsat = sat_e/(cr*cpf)
	ngroups_ =  Tsat/TFRAME # also =  sat_e/(TFRAME * cr * cpf)

	# to reach sat_e at given countrate in CP,  cr in CP is (cpf*cr)
	cp_e_per_frame = TFRAME * (cpf*cr)
	totnframes = cptot / cp_e_per_frame



	hdrfov["fluxfrac"] = (fract_31, "frac of WP power in this array")
	hdrfov["cpf"] = (cpf, "nrm peak wrt CLEAR psf total of 1")
	hdrfov["cpfnrm"] = (cpfnrm, "nrm peak wrt Inf. nrm psf total")
	hdrfov["nrmtot"] = (nrmtot, "WP nrm psf total is normed to mask thruput")
	hdrfov["tote"] = (tot_e, "total number of electrons required")
	hdrfov["totefull"] = (tot_e_full, "total number of electrons if CLEAR used")
	hdrfov["cptot"] = (cptot, "total in central pixel for TOT_E electrons")
	hdrfov["sate"] = (sat_e, "saturation limit in electrons")
	hdrfov["Tsat"] = (Tsat, "time for 1 integration, reaches saturation in CP")
	hdrfov["NINT"] = (nint_, "number of integrations required to reach TOT_E")
	hdrfov["NGROUPS"] = (np.ceil(ngroups_), "number of integrations required to reach TOT_E")
	hdrfov["TFRAME"] = (TFRAME, "single frame read time")
	hdrfov["totnfram"] = (totnframes, "total number of frames of data read")
	hdrfov["GAIN"] = (1.5, "electrons/ADU - not used here")

	ngroups, nint, nint_ceil = apt_values(ngroups_, nint_)

	report = "\t%s: \n\t\t%.1f%% of nrmflux in 31x31,  " % (filt, 100.0*fract_31) + \
	"%.2f%% peak frac wrt fullpsf.sum,  " % (100.0*cpf) + \
	"%.2f%% peak frac wrt nrmpsf.sum,  " % (100.0*cpfnrm) + \
	"%.2f%% nrmtot in 255x255  " % (100.0*nrmbig.sum()) + \
	"\n\t\tCRclearp %.1e  cpfCR %.3e  cptot %.1e  cp_e/frame %.1e  Tsat %5.1e, " %( cr, cpf*cr, cptot, cp_e_per_frame, Tsat) +\
        "\n\t\tNINT %f NINT(updated to account for lost partial group per integration) %f,"  %( nint_, nint) +\
        "\n\t\tNGROUPS %f totnframes %d,"  %( ngroups_,totnframes) +\
        "\n\t\t=====FOR APT INPUT: NGROUPS = %.0f NINT = %.0f=====" % ( ngroups, nint_ceil ) + \
	"\n\t\tdatacollect %.1f minutes = %.2f hours " % (nint_ * ngroups_ * TFRAME / 60.0, nint_ * ngroups_ * TFRAME / 3600.) #+ \
# 	"\n\t\tdatacollect %.5f minutes  "%(nint_ * ngroups_ * TFRAME / 60.0) #+ \
	#" %.5f hours  "                   %(nint_ * ngroups_ * TFRAME / 3600.0) + \
	#" %.5f days \n  "                 %(nint_ * ngroups_ * TFRAME / (24*3600.0))
	# times calculated with exact float values on nint_, ngroup_
	if np.isinf(nint_ceil):
           report = report + "\n\t\tBRIGHTNESS LIMIT EXCEEDED: some pixels may saturate, use NGROUPS = NINT = 1 in APT if desired."


	if return_params:
		params = [np.int(ngroups), nint, np.int(nint_ceil)]
		return report,params
	else:
		return report


def driver(tot_e, sat_e, mag_t):
	reports = []
	reports.append(generatePSF(filt=F277W, fov=31, osample=3, cr=cr_from_mag(mag_t, ZP[F277W]), tot_e=tot_e, sat_e=sat_e))
	reports.append(generatePSF(filt=F380M, fov=31, osample=3, cr=cr_from_mag(mag_t, ZP[F380M]), tot_e=tot_e, sat_e=sat_e))
	reports.append(generatePSF(filt=F430M, fov=31, osample=3, cr=cr_from_mag(mag_t, ZP[F430M]), tot_e=tot_e, sat_e=sat_e))
	reports.append(generatePSF(filt=F480M, fov=31, osample=3, cr=cr_from_mag(mag_t, ZP[F480M]), tot_e=tot_e, sat_e=sat_e))
	return reports
	
def cr_from_mag(M, zp=None):
	"""
	Flux of M magnitude star in counts/sec given a zero point
	zp is the zero point when using CLEARP
	"""
	return pow(10,-(M-zp)/2.5)

#"""
#	2015.10.05  Dean Hines asks for Jy brightness limits:
#
# Doris Daou's
# http://www.stsci.edu/hst/nicmos/tools/conversion_form.html
# Temperature of the blackbody =    5500.0
# 0mag  at L' band center 3.76um  Flux=2.50E+02 Jy
# 0mag  at M  band center 4.76um  Flux=1.50E+02 Jy
#
# so 
# 	use L1 mag for 277 380
# 	use M  mag for 430 480
#"""
#def vega2mjy():
#	jy_0mag = {"L1": 2.5e2, "M": 1.5e2}
#	blim_mag = {F277W: (6.7, "L1"), F380M: (3.8, "L1"), F430M: (3.1, "M"), F480M: (2.7, "M")}
#	print "filt   vega factor   mJy"
#	for bkey in sorted(blim_mag, reverse=True):
#		print bkey, " %.2f"%blim_mag[bkey][0], "%.1e"%pow(10, -blim_mag[bkey][0]),  "%.2e"%(1000.0 * pow(10, -blim_mag[bkey][0]) * jy_0mag[blim_mag[bkey][1]])
"""
>>> vega2mjy()filt   vega factor   mJy
F480M  2.70 2.0e-03 4.99e+02
F430M  3.10 7.9e-04 1.99e+02
F380M  3.80 1.6e-04 2.38e+01
F277W  6.70 2.0e-07 2.99e-02
>>> 
"""





# if __name__ == "__main__":
# 
# 	#vega2mjy()
# 	#sys.exit("Done Dean's calculation")
# 
# 	if len(sys.argv) > 1:
# 		MAG_T  = float(sys.argv[1])
# 	if len(sys.argv) > 2:
# 		TOT_E = float(sys.argv[2])
# 	tellme = driver(TOT_E,  SAT_E, MAG_T)
# 	print
# 	print "\t", "Using Kevin Volk's zero points, based on V band photometry.  Replace with K=9??"
# 	print "\t", F277W, ZP[F277W]
# 	print "\t", F380M, ZP[F380M]
# 	print "\t", F430M, ZP[F430M]
# 	print "\t", F480M, ZP[F480M]
# 
# 	print "\tTarget magnitude is %.2f, total  number of detected photons is %.0e\n" % (MAG_T, TOT_E)
# 	for t in tellme:
# 		print t
# 	print "\tTarget magnitude is %.2f, total  number of detected photons is %.0e\n" % (MAG_T, TOT_E)
# 	#print infostr
