from __future__ import print_function

import unittest, os, glob
import numpy as np
from astropy.io import fits

import driver_scene
from pyami.simcode import utils


class DriverSceneTestCase(unittest.TestCase):
    def setUp(self):
    
        # setup parameters for simulation, most are passed on to driver_scene
        verbose = 0
        overwrite = 1

        apply_dither = 0
        apply_jitter = 0
        include_detection_noise = 0
        uniform_flatfield = 1

        OVERSAMPLE = 11
        monochromatic_wavelength_m = 4.3e-6 

        create_calibrator = 0
        random_seed = 124
        flip = False
        mask = 'MASK_NRM'
        filter = 'F430M'
        n_image = 77

        filter_name = 'Monochromatic '+np.str(monochromatic_wavelength_m)

        # directory containing the static test data
        data_dir = os.path.join(os.path.dirname(__file__),'test_data')

        # define output directory for dynamically generated data
        out_dir = os.path.join(os.path.dirname(__file__),'tmp_data')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        name_seed = 'PSF_NIRISS_%s_%s'%(mask,filter)
        point_source_image_name = 'point_source_image.fits'
        point_source_image = os.path.join(out_dir,point_source_image_name) 

        psf_image_name = name_seed + '_reference.fits'
        psf_image = os.path.join(data_dir,psf_image_name)
        psf_image_without_oversampling = os.path.join(data_dir,psf_image_name.replace('.fits','_without_oversampling.fits'))
    
        # generate delta-function like image to feed as point-source into scene_sim    
        n_image2 = n_image * OVERSAMPLE
        data = np.zeros((n_image2,n_image2))
        bright_pixel_index = (np.int(np.floor(n_image2/2.)),np.int(np.floor(n_image2/2.)))
        data[bright_pixel_index] = 1.
        fits.writeto(point_source_image, data, clobber=True)
        
        
        # following code can be used to regenerate reference PSF file, but it depends on nrm_analysis        
        if 0:
            from nrm_analysis.fringefitting.LG_Model import NRM_Model
            jw = NRM_Model(mask='jwst',holeshape="hex",flip=flip)
            jw.simulate(fov=n_image, bandpass=monochromatic_wavelength_m, over=OVERSAMPLE, pixel = mas2rad(PIXELSCL_arcsec*1000))
            fits.writeto(psf_image,jw.psf_over, clobber=True)
            header = fits.getheader(psf_image)
            header['PIXELSCL'] = PIXELSCL_arcsec/OVERSAMPLE
            header['FILTER'] = filter_name
            header['PUPIL'] = mask 
            fits.update(psf_image,jw.psf_over/10000./28., header=header)
    
            # PSF without oversampling
            fits.writeto(psf_image_without_oversampling,jw.psf, clobber=True)
            header = fits.getheader(psf_image_without_oversampling)
            header['PIXELSCL'] = PIXELSCL_arcsec
            header['FILTER'] = filter_name
            header['PUPIL'] = mask 
            fits.update(psf_image_without_oversampling,jw.psf, header=header)
    
    
        NGROUP = 1
        NINT = 1
        COUNTRATE = 5e8

#     run driver_scene to generate target and calibrator images    
        driver_scene.main(['--output_absolute_path','%s'%out_dir,\
                           '--overwrite','%d'%overwrite,'-utr','0', '-f','%s'%filter,'-v','%d'%verbose,'--apply_dither','%d'%apply_dither,'--apply_jitter','%d'%apply_jitter,\
                           '-p',psf_image,'-s','%s'%point_source_image,'-os','%d'%OVERSAMPLE,'-cr','%e'%COUNTRATE,'--include_detection_noise','%d'%include_detection_noise,\
                           '-tag',('%3.2e_NGROUP%d'%(COUNTRATE,NGROUP)).replace('.','p'),'--uniform_flatfield','%d'%uniform_flatfield,'--random_seed','%d'%random_seed,\
                           '--create_calibrator','%d'%create_calibrator,'--nint','%d'%NINT,'--ngroups','%d'%NGROUP])
        
        # directory where simulated data was written
        filter_dir = os.path.join(out_dir,filter)
        self.filter_dir = filter_dir    

        # list of files produced for target
        file_list = glob.glob(os.path.join(filter_dir,'%s*.fits' % 't_' ));
            
        self.simulated_image = file_list[0]
        self.psf_image_without_oversampling = psf_image_without_oversampling

                
    def test_image_ratio(self):
    
        # load simulated image  data
        t_sim,t_sim_header = fits.getdata(self.simulated_image,header=True)
        
        # load PSF data
        psf_without_oversampling = fits.getdata(self.psf_image_without_oversampling)
        
        ratio_data = t_sim/psf_without_oversampling;
        fits.writeto(os.path.join(self.filter_dir,'ratio_image.fits'), ratio_data, clobber=True)
		#         print('Noise of ratio image %3.3f' % np.std(ratio_data))
        
        self.assertTrue(np.std(ratio_data) < 10., 'driver_scene simulation of point source is inaccurate')

if __name__ == '__main__':
    unittest.main()
    
    
    
    
    