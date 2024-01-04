The code is run by typing in terminal: 
SPADE.py config_file Min_Dec Max_Dec

The pipeline requires the following information:
Headers: FILTER, EXPTIME, TIME,RA,DEC,NAXIS1,NAXIS2
Sextractor Params: NUMBER,FLUX_ISOCOR,MAG_ISOCOR,MAGERR_ISOCOR,FLUX_APER(1),FLUXERR_APER(1), X_IMAGE, Y_IMAGE, FLAGS, FWHM_IMAGE, ELONGATION, FLUX_RADIUS, SNR_WIN, VIGNET, MAG_PSF

An example config_file is placed in this pipeline. Build one for your own equipment. It allows you to control the entire pipeline, enabling/disabling different parts and control different parameters that affect them.

Pipeline explanation:
An illustration of the pipeline is available on this folder.



