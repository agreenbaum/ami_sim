#!/usr/bin/env python
# Code to run observation simulation for binary:
# 1. run ETC to get desired electron count
# 2. simulate AMI data of a binary with desired separation and flux ratio
# 
# 2016-02-15 Johannes Sahlmann
# based on ETC and binary simulation codes from Deepashri Thatte and Anand Sivaramakrishnan


import sys, os, argparse
import numpy as np
from astropy.io import fits 


def main(argv):

    parser = argparse.ArgumentParser(description="This script provides a basic example of how to use the NIRISS-AMI ETC and binary simulation code.")
    parser.add_argument('-t','--targetDir',  type=str, default='niriss-ami_out/', help='path (relative to home directory) to output directory for file writing')
    parser.add_argument('-o','--overwrite',  type=int, default='0', help='Overwrite yes/no, must be either 0 or 1', choices=[0,1])
    parser.add_argument('-st','--spectralType',  type=str, default='A0V', help='Spectral type of target, only certain types will work')
    parser.add_argument('-fr','--fluxRatio',  type=float, default='0.16', help='Flux ratio between primary and secondary in selected filter')
    parser.add_argument('-dx','--deltaX',  type=float, default='1.2', help='Separation in X in PIXEL between primary and secondary')
    parser.add_argument('-dy','--deltaY',  type=float, default='1.0', help='Separation in Y in PIXEL between primary and secondary')
    parser.add_argument('-utr','--uptheramp',  type=int, default='0', help='Generate up-the-ramp fits file? yes/no, must be either 0 or 1', choices=[0,1])
    
    parser.add_argument('filter', type=str, help='filter to be used, must be one of F277W,F380M,F430M,F480M', choices=["F277W", "F380M", "F430M", "F480M"])
    parser.add_argument('targetMagnitude', type=float, help='Target apparent magnitude in selected filter')
    parser.add_argument('totalElectrons', type=float, help='Requested total number of collected electrons')
    parser.add_argument('-O','--oversample', type=int, default='11', help='detector pixel oversampling (odd)', choices=range(1,12,2))
	
    
    args = parser.parse_args(sys.argv[1:])

    # JSA: only for development purposes
    
    if 0==1:
        import pyami.etc.NIRISSami_apt_calc_v2_JSA as etc        
        import pyami.simcode.make_binary as binsim        
        if ('pyami' in sys.modules):
            import pyami; reload(pyami); 
            reload(pyami.etc.NIRISSami_apt_calc_v2_JSA);
            reload(pyami.simcode.make_binary);
            
    import pyami.etc.NIRISSami_apt_calc_v2_JSA as etc 
    import pyami.simcode.utils as U   
    import pyami.simcode.make_binary as binsim        
    
    pathname = os.path.dirname(sys.argv[0])
    fullPath = os.path.abspath(pathname)
    pyamiDataDir = fullPath + '/pyami/etc/NIRISSami_apt_calcPSF/';

    targetDir = args.targetDir;   #     /astro/projects/JWST/NIRISS/AMI/simulatedData/        
    outDir0 = os.getenv('HOME') + '/' + targetDir;

    overwrite = args.overwrite; #0;
    osample = osample = args.oversample

    # ETC inputs
    MAG_T  = args.targetMagnitude; #11.0
    TOT_E = args.totalElectrons; # 1e6
    SAT_E = 35e3
    
    filt = args.filter
    sptype = args.spectralType

    uptheramp = args.uptheramp

    fluxratio = args.fluxRatio
    dx_pix = args.deltaX; # in PIXEL
    dy_pix = args.deltaY; # in PIXEL

    filters = np.array([filt])
    

    for filter in filters:
        #   run ETC prototype
        mag_t = MAG_T
        tot_e = TOT_E
        sat_e = SAT_E

        report, params = etc.generatePSF(filt=filter, fov=31, osample=3, cr=etc.cr_from_mag(mag_t, U.ZP[filter]), tot_e=tot_e, sat_e=sat_e, SRC = sptype, return_params=1, DATADIR=pyamiDataDir)
        ngroups, nint, nint_ceil = params;
        print report
        print "\tTarget magnitude is %.2f, total  number of detected photons is %.0e\n" % (MAG_T, TOT_E)
        
        #   make binary simulation and generate images  
        if 1==1:
            fov = U.amisubfov
            tframe = 0.0745 # SUB80 NISRAPID
            dim = fov/2.0
            print "oversampling set in top level driver to 11" 
            trials = 1
        
            #             fluxratio = 0.16
            #define offsets of the companion from the star at the center
            offset_x = dx_pix * osample; #-14.0   #companion offset in x in oversampled pixels
            offset_y = dy_pix * osample; #11.0    #companion offset in y in oversampled pixels

            outDir = outDir0 + '%s/' % (filter);
            tmpDir = outDir0 + 'tmp/'
            for dd in [outDir,tmpDir]:
                if not os.path.exists(dd):
                    os.makedirs(dd)

            if not os.path.exists(outDir):
                os.makedirs(outDir)
            starArrayFile = outDir+'star_array_fov80_%s.fits'%filter.lower();
        
            # FEEDER FOR SIMULATION - star_array
            #f ( (not os.path.isfile(starArrayFile)) or (overwrite == 1) ) :
            if not os.path.isfile(starArrayFile):
                print "Driver creating or overwriting PSF file", starArrayFile
                star_array,star_header = binsim.generate_starPSF(FILTER=filter,fov=80, osample=11, spectraltype=sptype)
                fits.PrimaryHDU(data=star_array, header=star_header).writeto(starArrayFile, clobber = True)
            else:
                print "Driver found appropriate PSF file on disk:", starArrayFile
                star_array,star_header = fits.getdata(starArrayFile,header=True)

            calstar_array = star_array[((fov+4)/2-dim)*osample:((fov+4)/2+dim)*osample,((fov+4)/2-dim)*osample:((fov+4)/2+dim)*osample]

            # --
            # DEFINE DITHER POINTING in det pixels
            #
            ipsoffset = U.ips_size//2 - (U.amisubfov)//2
            print "ipsoffset %d = U.ips_size//2 %d - (U.amisubfov)//2 %d" % (ipsoffset, U.ips_size//2, (U.amisubfov)//2)

            x_dith = [17,57,17,53]; # see email Deepashri 2016-02-10

            x_dith[:] = [a + ipsoffset for a in x_dith]
            x_dith[:] = [(a*osample - osample//2-1) for a in x_dith]

            y_dith = [58,58,18,18]
            y_dith[:] = [a + ipsoffset for a in y_dith]
            y_dith[:] = [(a*osample - osample//2-1) for a in y_dith]

            dithers = len(x_dith)

            #======================================
            #[2] CREATE NOISELESS BINARY STAR ARRAY
            #======================================


            #define the edges of star PSF
            a = 2 * osample
            b = (2 + fov) * osample 

            # FEEDER FOR SIM!!!
            try: 
                binarystar_array = star_array[a:b, a:b] + fluxratio *  star_array[a+offset_y:b+offset_y, a+offset_x:b+offset_x]
            except ValueError:
                sys.exit("Error: binary separation too great?")

            fits.writeto(outDir+'psf_binarystar.fits',binarystar_array, clobber = True)

            magzeropoint = U.ZP[filter]
            # additional arguments passed to simulate_skydata:   >>>> NEEDS CLEANUP <<<<
            kwargs = {'offset_x': offset_x, 'offset_y':offset_y, 'fluxratio':fluxratio, 'starmag':MAG_T , 'magzeropoint':magzeropoint, 'sptype':sptype, 'fov':fov, 'osample':osample, 'dim':dim};

            if uptheramp:
                binsim.simulate_skydata(trials, binarystar_array, "tcube", dithers, x_dith, y_dith, ngroups, nint_ceil, U.tframe, filt=filter, outDir = outDir,tmpDir=tmpDir, **kwargs)
                binsim.simulate_skydata(trials, calstar_array   , "ccube", dithers, x_dith, y_dith, ngroups, nint_ceil, U.tframe, filt=filter, outDir = outDir,tmpDir=tmpDir, **kwargs)
            else: 
                frametime = ngroups * U.tframe # What is the 'total photon collection time'? Put that here
                ngroups = 1 # and just perform the reset-readout and the final readout
                binsim.simulate_skydata(trials, binarystar_array, "tcube", dithers, x_dith, y_dith, ngroups, nint_ceil, frametime, filt=filter, outDir = outDir,tmpDir=tmpDir, **kwargs)
                binsim.simulate_skydata(trials, calstar_array   , "ccube", dithers, x_dith, y_dith, ngroups, nint_ceil, frametime, filt=filter, outDir = outDir,tmpDir=tmpDir, **kwargs)


    
if __name__ == "__main__":
    main(sys.argv[1:])    
    
sys.exit(0)


