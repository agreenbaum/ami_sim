Get the source files for the AMI data simulation programs, e.g.: 

> git clone https://github.com/agreenbaum/ami_sim

> cd ~/ami_sim

You should be able to run three programs:
> python driver_binary.py --help

> python driver_scene.py --help

> python ami_etc.py --help


# driver_binary  -  creates simulated binary star data [better to use driver_scene for more general simulation]
----------------

> python driver_binary.py F430M 11.0 2e6 -t simulatedData/ -o 1 -fr 0.16 -dx 2.0 -dy 2.0


It will (if needed) create a new directory (~/simulatedData) in your home
directory which contains the output data.

It requires WebbPSF and pysynphot python modules and supporting data files.



# driver_scene  -  creates simulated data of sky.fits using psf.fits, e.g.:
---------------

> python driver_scene.py -t simulatedData/ -o 0 -utr 0 -f F430M \
                         -p [absolute path to psf file]/psf.fits \
                         -s [absolute path to psf file]/sky.fits \
                         -os 11 -I 10 -G 2 [-c 0]  -cr 2.1e6  -v 1\
                         --random_seed 42

You need to create a directory (eg simulatedData in the above example) in your home directory for the saved files. You also need to point to: 

	- Absolute path to the sky scene (eg sky.fits, normalized on-the-fly to cr in photons/s on 25m^2 in filter bandpass)
	- Absolute path to the PSF file (eg psf.fits) (PSF.sum = NRM area / full aperture area)

The output "data cube"" files will be named using your input filename,  i.e. 
	- t_sky__psf.fits  for the target
	- c_sky__psf.fits for the calibrator (-c flag omitted, or -c 1)


** Only oversampling of 11 is tested so far - use the flag "-os 11" **
** Oversampling 5 testing in progress **

NOTES: 
	sky.fits must be square, and an odd number (<80) of 65 mas detector pixels on a side
	Generate your sky.fits using an odd oversampling (eg 11)
	sky.fits in units of detected e-/s using 25m^2 collecting area through the filter
	Create a PSF file (eg using WebbPSF) with the shape and pixel scale of sky.fits
	To generate a calibrator star data place the total of your sky.fits into a single oversampled pixel in the same-sized cal.fits input file.
	The user chooses NINT and NGROUP.  The accompanying ETC can estimate operational values, or choose your own.
		NGROUP+1 is the number of detector readouts up the ramp
		NINT is the number of full ramps

The following two commands simulate a binary companion without a disk (user-created sky scene) and a binary companion with a disk (also user-created sky scene sky scene).  Both simulations generate a calibrator observation also ("c_*.fits").  For the -cr (--countrate) parameter, use the "CRclearp"" number that ami_etc.py gives you for a given point source target magnitude and filter.

	> python driver_scene.py -t simDataLkCa15/bin/ -o 0 -utr 0 -f F430M \
          -p PSF_5x_MASK_NRM_F430M_opd162_3_A0V.fits -s binary_cr0.01_pa0.0.fits \
          -O 5 -I 3 -G 30 -utr 0 -c 1 -cr 3e6

	> python driver_scene.py -t simDataLkCa15/bnd/ -o 0 -utr 0 -f F430M \
          -p PSF_5x_MASK_NRM_F430M_opd162_3_A0V.fits -s binaryanddisk_cr0.01_pa0.0.fits \
          -O 5 -I 3 -G 30 -utr 0 -c 1 -cr 3e6


# ami_etc  -   ETC script
----------

> python ami_etc.py --help

and as an example, for an A0V star

> python ami_etc.py f277w 7.5 1e8 

 or
 
> python ami_etc.py f277w 7.5 1e8  -st M5V
 
Only these two spectral types are supported.

This will report exposure parameters for the Astronomer's Proposal Tool (APT) and other useful quantities using the NIRISS team's private ETC, which is not the official STScI ETC



Maintainers: Deepashri Thatte, Anand Sivaramakrishnan, Johannes Sahlmann, Alexandra Greenbaum

