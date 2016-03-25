Get the source files for the AMI data simulation programs, e.g.: 

> git clone https://github.com/agreenbaum/ami_sim

cd into the ami_sim directory

You should be able to run three programs:

> python driver_binary.py --help
> python driver_scene.py --help
> python ami_etc.py --help


# driver_binary  -  creates simulated binary star data
----------------

> python driver_binary.py F430M 11.0 2e6 -t simulatedData/ -o 1 -fr 0.16 -dx 2.0 -dy 2.0


It will (if needed) create a new directory (~/simulatedData) in your home
directory which contains the output data.

It requires WebbPSF and pysynphot python modules and supporting data files.



# driver_scene  -  creates simulated data of sky.fits using psf.fits, e.g.:
---------------

> python driver_scene.py -t simulatedData/ -o 0 -utr 0 -f \
    F430M -p psf.fits -s sky.fits -O 11 -I 10 -G 2
    
> python driver_scene.py -t simulatedData/ -o 0 -utr 0 -f 
    F430M -p psf.fits -s cal.fits -O 11 -I 10 -G 2

	** Only oversampling of 11 is tested so far - use the flag "-O 11" **
	
These will (if needed) create a new directory (~/simulatedData) in your home
directory which contains the output simulated data file named with your input fits filename roots,  i.e. sky__psf.fits.


NOTES: 
	sky.fits must be square, and an odd number (<80) of 65 mas detector pixels on a side
	Generate your sky.fits using an odd oversampling (eg 11)
	sky.fits in units of detected e-/s using 25m^2 collecting area through the filter
	Create a PSF file (eg using WebbPSF) with the shape and pixel scale of sky.fits
	To generate a calibrator star data place the total of your sky.fits into a single oversampled pixel in the same-sized cal.fits input file.
	The user chooses NINT and NGROUP.  The accompanying ETC can estimate operational values, or choose your own.
		NGROUP+1 is the number of detector readouts up the ramp
		NINT is the number of full ramps


# ami_etc  -   ETC script
----------

> python ami_etc.py --help

and as an example, for an A0V star

> python ami_etc.py f277w 7.5 1e8 

 or
 
> python ami_etc.py f277w 7.5 1e8  -st M5V
 
Only these two spectral types are supported.

This will report exposure parameters for the Astronomer's Proposal Tool (APT) and other useful quantities using the NIRISS team's private ETC, which is not the official STScI ETC



Maintainers (in alphabetic order): Greenbaum, Sahlmann, Sivaramakrishnan, Thatte

