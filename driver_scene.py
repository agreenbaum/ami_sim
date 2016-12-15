#!/usr/bin/env python
# Code to run observation simulation for binary:
# 1. run ETC to get desired electron count
# 2. simulate AMI data of a binary with desired separation and flux ratio
# 
# 2016-02-15 Johannes Sahlmann
# based on ETC and binary simulation codes from Deepashri Thatte and Anand Sivaramakrishnan
#
# 2016-09-27 anand@stsci.edu 
#     Create calibrator observation with same readout as target,  same total flux
#     Places calibrator at location of first maximum of target (in np.where(tgt==tgt.max()) list)
#     rework driver_scene to create cal and tgt data cubes
#     rework make_scene to match reworked make_binary 
#         (utils create_ramp & create_integration changed,
#         new up-the-ramp handling to match binary simulation)
#     name cal and tgt data cubes c_* and t_*


import sys, os, argparse
import numpy as np
from astropy.io import fits 


def main(argv):

    print argv

    parser = argparse.ArgumentParser(description="This script simulates the observation of a sky scene with a (Webb)PSF supplied by the user")
    parser.add_argument('-t','--targetDir',  type=str, default='niriss-ami_out/', help='output dir path (relative to ~)')
    parser.add_argument('-o','--overwrite',  type=int, default='0', help='overwrite yes/no, default 0 (no)', choices=[0,1])
    parser.add_argument('-utr','--uptheramp',  type=int, default='0', help='generate up-the-ramp fits file? yes/no, default 0 (no)', choices=[0,1])
    parser.add_argument('-f', '--filter', type=str, help='filter name (upper/lower case)', choices=["F277W", "F380M", "F430M", "F480M"])
    parser.add_argument('-p','--psf', type=str, help='oversampled PSF fits file. Spectral type set in this')
    parser.add_argument('-s','--sky', type=str, help='oversampled sky scene fits file, normalized to sum to unity')
    parser.add_argument('-O','--oversample', type=int, help='sky scene oversampling (odd)', choices=range(1,12,2))
    parser.add_argument('-I','--nint', type=int, default=1, help='number of integrations (IR community calls these exposures sometimes)')
    parser.add_argument('-G','--ngroups', type=int, default=1, help='number of up-the-ramp readouts')
    parser.add_argument('-c','--calibrator', type=int, default=1, help='create calibrator observation yes/no default 1 (yes)', choices=[0,1])
    parser.add_argument('-cr','--countrate', type=float, help='Photon count rate on 25m^2 per sec in the bandpass (CRclearp in ami_etc output)',)
    
    args = parser.parse_args(sys.argv[1:])

    import pyami.simcode.make_scene as scenesim        
    import pyami.simcode.utils as U        

    pathname = os.path.dirname(sys.argv[0])
    fullPath = os.path.abspath(pathname)

    targetDir = args.targetDir 
    outDir0 = os.getenv('HOME') + '/' + targetDir;

    overwrite = args.overwrite #0;
    uptheramp = args.uptheramp
    calibrator = args.calibrator

    print "countrate input as %.2e photons/sec on 25m^2 primary in filter bandpass" % args.countrate
    countrate = args.countrate

    filt = args.filter
    psffile = args.psf
    skyfile = args.sky
    osample = args.oversample

    nint = args.nint        # TBD: calculate internally to save the user prep time doing ETC work
    ngroups = args.ngroups  # TBD: calculate internally to save the user prep time doing ETC work
    # rebin sky_conv_psf image to detector scale, use max of detector array to calculate nint, ngroups, data-collect-time

        
    # generate images  
    if 1==1:
        print "oversampling set in top level driver to 11" 
        trials = 1
        
        outDir = outDir0 + '%s/' % (filt);
        tmpDir = outDir0 + 'tmp/'
        # NB outDir must exist to contain input files - clean up organization later?
        for dd in [outDir,tmpDir]:
            if not os.path.exists(dd):
                os.makedirs(dd)


        # FEEDER FOR SIMULATION - read in pre-made psf made by WebbPSF (or any other way)
        # File sizes: 
        psfdata, psfhdr = fits.getdata(outDir0+psffile, header=True)
        skydata, skyhdr = fits.getdata(outDir0+skyfile, header=True)
        skydata = skydata / skydata.sum()  # normalize sky data total to unity!
        skydata = skydata * countrate
        print "psfdata", psfdata.shape, "totals %.2e (NRM throughput / full aperture throughput)"%psfdata.sum()
        print "skydata", skydata.shape, "totals %.2e (photons / s on 25^m in band)"%skydata.sum()

        # Note: to generate a calibration star observation, use a 'delta function' single positive
        # pixel in an otherwise zero-filled array as your sky fits file.  No longer do we m
        # match total CR in skydata, but we match the CR to that of the dominant point source
        # in skydata.  Thus if the skydata point source accounts for 0.99 of the skydata array 
        # then we place a unit delta function into caldata, and request a count rate of 
        # 0.99 * the user-requested count rate when creating calibrator data.
        # WARNING - WE ASSUME A SINGLE DOMINANT POINT SOURCE IN THE SKY DATA.

        caldata = np.zeros(skydata.shape, np.float64)
        maxloc = np.where(skydata==skydata.max())
        ptsrcfraction = skydata[maxloc]/skydata.sum()
        #print maxloc, ptsrcfraction
        caldata[maxloc[0][0], maxloc[1][0]] = ptsrcfraction * countrate
        #print "caldata[maxloc],  skydata[maxloc], ratio: ",  (caldata[maxloc],  skydata[maxloc], 
        print """

          We match the CR to that of the dominant point source in skydata.
          Thus if the skydata point source accounts for 0.99 of of target flux,
          then we place a unit delta function into caldata, and request a count
          rate of 0.99 * the user-requested count rate for creating calibrator
          data.  
          
          WARNING - WE ASSUME A SINGLE DOMINANT POINT SOURCE IN THE SKY DATA.  

          """

        """ fov = 80
        dim = fov/2.0 """

        
        # --
        # DEFINE DITHER POINTING in det pixels
        ipsoffset = U.ips_size//2 - (skydata.shape[0]//osample)//2
        x_dith, y_dith = [(skydata.shape[0]//2)/osample + ipsoffset,], \
                         [(skydata.shape[0]//2)/osample + ipsoffset,]
        dithers = len(x_dith)
        print "x_dith, y_dith", x_dith, y_dith

        # now convert to oversampled pixels for the calculation:
        x_dith[:] = [(x*osample - osample//2+1) for x in x_dith]
        y_dith[:] = [(y*osample - osample//2+1) for y in y_dith]

        cubename = "t_" + skyfile.replace(".fits","__") + psffile
        scenesim.simulate_scenedata(trials, 
                                    skydata, psfdata, psfhdr, cubename, osample,
                                    dithers, x_dith, y_dith,
                                    ngroups, nint, U.tframe, filt,
                                    outDir, tmpDir, uptheramp)

        if calibrator:
            cubename = "c_" + skyfile.replace(".fits","__") + psffile
            scenesim.simulate_scenedata(trials, 
                                        caldata, psfdata, psfhdr, cubename, osample,
                                        dithers, x_dith, y_dith,
                                        ngroups, nint, U.tframe, filt,
                                        outDir, tmpDir, uptheramp)

if __name__ == "__main__":
    print sys.argv
    main(sys.argv[1:])    
    
sys.exit(0)

