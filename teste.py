def run_sextractor (image, cat_out, file_config, file_params, pixscale,
                    header, fit_psf=False, return_fwhm_elong=True, fraction=1.0,
                    fwhm=5.0, update_vignet=False, imtype=None, fits_mask=None, 
                    npasses=2, tel=None, set_zogy=None, nthreads=0,
                    Scorr_mode=None, image_analysis=None, std_Scorr=1):

    
    """Function that runs SExtractor on [image], and saves the output
       catalog in [cat_out], using the configuration file
       [file_config] and the parameters defining the output recorded
       in the catalogue [file_params]. If [fit_psf] is True,
       SExtractor will perform PSF fitting photometry using the PSF
       built by PSFex. If [return_fwhm_elong] is True, an estimate of
       the image median FWHM and ELONGATION and their standard
       deviations are returned using SExtractor's seeing estimate of
       the detected sources; if False, it will return zeros. If
       [fraction] is less than the default 1.0, SExtractor will be run
       on a fraction [fraction] of the area of the full
       image. Sextractor will use the input value [fwhm], which is
       important for the star-galaxy classification. If [save-bkg] is
       True, the background image, its standard deviation and the
       -OBJECTS image (background-subtracted image with all objects
       masked with zero values), all produced by SExtractor, will be
       saved. If [set_zogy.bkg_method] is not set to 1 (use
       SExtractor's background), then improve the estimates of the
       background and its standard deviation.

    """

    if get_par(set_zogy.timing,tel): t = time.time()
    log.info('executing run_sextractor ...')

    base = image.replace('.fits','')
    
    # make copy of input image, as input image will get background
    # subtracted
    if get_par(set_zogy.keep_tmp,tel):
        image_orig = '{}_orig.fits'.format(base)
        if not os.path.exists(image_orig):
            shutil.copy2 (image, image_orig)

        
    # if fraction less than one, run SExtractor on specified fraction of
    # the image
    if fraction < 1:

        # read input image and header
        data, header = read_hdulist (image, get_header=True)
        # get input image size from header
        xsize, ysize = header['NAXIS1'], header['NAXIS2']
        
        # determine cutout from [fraction]
        center_x = int(xsize/2+0.5)
        center_y = int(ysize/2+0.5)
        halfsize_x = int((xsize * np.sqrt(fraction))/2.+0.5)
        halfsize_y = int((ysize * np.sqrt(fraction))/2.+0.5)
        data_fraction = data[center_y-halfsize_y:center_y+halfsize_y,
                             center_x-halfsize_x:center_x+halfsize_x]

        # write small image to fits
        image_fraction = '{}_fraction.fits'.format(base)
        header['DATEFILE'] = (Time.now().isot, 'UTC date of writing file')
        fits.writeto(image_fraction, data_fraction.astype('float32'), header, 
                     overwrite=True)

        # make image point to image_fraction
        image = image_fraction
        cat_out = '{}_ldac_fraction.fits'.format(base)


    # the input fwhm determines the SEEING_FWHM (important for
    # star/galaxy separation) and the radii of the apertures used for
    # aperture photometry. If fwhm is not provided as input, it will
    # assume fwhm=5.0 pixels.
    fwhm = float('{:.2f}'.format(fwhm))
    # determine seeing
    seeing = fwhm * pixscale
    # prepare aperture diameter string to provide to SExtractor 
    apphot_diams = np.array(get_par(set_zogy.apphot_radii,tel)) * 2 * fwhm
    apphot_diams_str = ','.join(apphot_diams.astype(str))

    if update_vignet or fits_mask is not None:

        # create updated parameters file starting from original file
        file_params_edited = '{}_params.txt'.format(base)
        shutil.copy2(file_params, file_params_edited)
    
        # update size of VIGNET
        if update_vignet:
            size_vignet = get_par(set_zogy.size_vignet,tel)
            # write vignet_size to header
            header['S-VIGNET'] = (size_vignet, '[pix] size square VIGNET used '
                                  'in SExtractor')
            # append the VIGNET size to the temporary SExtractor
            # parameter file created above
            size_vignet_str = str((size_vignet, size_vignet))
            with open(file_params_edited, 'a') as myfile:
                myfile.write('VIGNET{}\n'.format(size_vignet_str))

        if fits_mask is not None:
            # and add line in parameter file to include IMAFLAG_ISO
            with open(file_params_edited, 'a') as myfile:
                myfile.write('IMAFLAGS_ISO\n')

        # use edited version for parameter file
        file_params = file_params_edited


    # check if background was already subtracted from image
    if 'BKG-SUB' in header:
        bkg_sub = header['BKG-SUB']
    else:
        header['BKG-SUB'] = (False, 'sky background was subtracted?')
        bkg_sub = False


    # determine whether interpolation is allowed across different
    # channels in [mini2back] using function get_Xchan_bool
    chancorr = get_par(set_zogy.MLBG_chancorr,tel)
    interp_Xchan = get_Xchan_bool (tel, chancorr, imtype)
    interp_Xchan_std = get_Xchan_bool (tel, chancorr, imtype, std=True)


    # get gain from header
    gain = read_header(header, ['gain'])

    
    log.info('background already subtracted?: {}'.format(bkg_sub))

    # do not apply weighting
    apply_weight = False

    # initialize cmd_dict dictionary that is filled in below for the
    # different options, which will be converted to a cmd list to be
    # executed from the unix command line; make sure it is an
    # insertion-ordered dictionary (default for dictionaries starting
    # from python3.6)
    cmd_dict = collections.OrderedDict()
    cmd_dict['source-extractor'] = image
    if image_analysis is not None:
        cmd_dict['source-extractor'] = '{},{}'.format(image, image_analysis)
    cmd_dict['-c'] = file_config
    cmd_dict['-BACK_TYPE'] = 'AUTO'
    cmd_dict['-BACK_VALUE'] = '0.0'
    cmd_dict['-BACK_SIZE'] = str(get_par(set_zogy.bkg_boxsize,tel))
    cmd_dict['-BACK_FILTERSIZE'] = str(get_par(set_zogy.bkg_filtersize,tel))
    cmd_dict['-BACKPHOTO_TYPE'] = get_par(set_zogy.bkg_phototype,tel).upper()
    cmd_dict['-VERBOSE_TYPE'] = 'QUIET'
    cmd_dict['-CATALOG_NAME'] = cat_out
    cmd_dict['-PARAMETERS_NAME'] = file_params
    cmd_dict['-PIXEL_SCALE'] = str(pixscale)
    cmd_dict['-SEEING_FWHM'] = str(seeing)
    cmd_dict['-PHOT_APERTURES'] = apphot_diams_str
    cmd_dict['-NTHREADS'] = str(nthreads)
    cmd_dict['-FILTER_NAME'] = get_par(set_zogy.sex_det_filt,tel)
    starnnw_name = '{}default.nnw'.format(get_par(set_zogy.cfg_dir,tel))
    cmd_dict['-STARNNW_NAME'] = starnnw_name


    for npass in range(npasses):
            
        if npass==0:        

            if apply_weight:
                cmd_dict['-WEIGHT_TYPE'] = 'BACKGROUND'
                
        else:
            
            # for FWHM estimate, or if background method is set to 1
            # (SExtractor method) or background had already been
            # subtracted, or Scorr_mode is not None: no need to do
            # multiple passes
            if (return_fwhm_elong or
                get_par(set_zogy.bkg_method,tel)==1 or
                bkg_sub or
                Scorr_mode is not None):
                break

            log.info ('running 2nd pass of SExtractor')

            # save catalog of 1st/initial run to different filename
            if get_par(set_zogy.keep_tmp,tel):
                os.rename (cat_out, cat_out.replace('ldac', 'ldac_init'))


            if apply_weight:
                cmd_dict['-WEIGHT_TYPE'] = 'MAP_RMS'
                cmd_dict['-WEIGHT_IMAGE'] = fits_bkg_std


        # add commands to produce BACKGROUND, BACKGROUND_RMS and
        # background-subtracted image with all pixels where objects
        # were detected set to zero (-OBJECTS). These are used to
        # build an improved background map.
        if (not return_fwhm_elong and
            npass==0 and
            not bkg_sub and
            Scorr_mode is None):

            fits_bkg = '{}_bkg.fits'.format(base)
            fits_bkg_std = '{}_bkg_std.fits'.format(base)
            fits_objmask = '{}_objmask.fits'.format(base)
            cmd_dict['-CHECKIMAGE_TYPE'] = 'BACKGROUND,BACKGROUND_RMS,-OBJECTS'
            image_names = '{},{},{}'.format(fits_bkg, fits_bkg_std, fits_objmask)
            cmd_dict['-CHECKIMAGE_NAME'] = image_names

        elif bkg_sub or npass>0:
            
            # if background was already subtracted or this is the 2nd
            # pass, still let SExtractor produce the object mask image
            fits_objmask = '{}_objmask.fits'.format(base)
            cmd_dict['-CHECKIMAGE_TYPE'] = '-OBJECTS'
            image_names = '{}'.format(fits_objmask)
            cmd_dict['-CHECKIMAGE_NAME'] = image_names


        # in case of fwhm/elongation estimate: only consider higher
        # S/N detections
        if return_fwhm_elong:
            cmd_dict['-DETECT_THRESH'] = str(get_par(
                set_zogy.fwhm_detect_thresh,tel))

        elif Scorr_mode is not None:
            # detection threshold for transient detection of at least
            # 3 pixels above 2/3 * set_zogy.transient_nsigma; scale
            # with std_Scorr as source extractor will use actual Scorr
            # STD
            det_th = 0.66 * get_par(set_zogy.transient_nsigma,tel) / std_Scorr
            ana_th = det_th
            # increase in 'init' mode - exact detection threshold is
            # only important for background-subtracted image
            if Scorr_mode=='init':
                det_th *= 5

            cmd_dict['-DETECT_THRESH'] = str(det_th)
            cmd_dict['-DETECT_MINAREA'] = '3'
            cmd_dict['-DETECT_MAXAREA'] = '2000'
            cmd_dict['-ANALYSIS_THRESH'] = str(ana_th)
            cmd_dict['-CATALOG_TYPE'] = 'FITS_1.0'
            cmd_dict['-FILTER'] = 'N'
            #cmd_dict['-FILTER_THRESH'] = '2.0'
            #cmd_dict['-FILTER_NAME'] = '/Users/pmv/ZOGY/Config/gauss_3.0_5x5.conv'
            cmd_dict['-DEBLEND_NTHRESH'] = '8'
            cmd_dict['-DEBLEND_MINCONT'] = '0.1'
            # save background-subtracted Scorr image only in 'init' mode
            if Scorr_mode=='init':
                fits_bkgsub = '{}_bkgsub.fits'.format(base)
                cmd_dict['-CHECKIMAGE_TYPE'] = '-BACKGROUND'
                cmd_dict['-CHECKIMAGE_NAME'] = fits_bkgsub
                

        # provide PSF file from PSFex
        if fit_psf:
            cmd_dict['-PSF_NAME'] = '{}_psf.fits'.format(base)


        # provide mask image and type if not None
        if fits_mask is not None:
            log.info('mask: {}'.format(fits_mask))
            cmd_dict['-FLAG_IMAGE'] = fits_mask
            cmd_dict['-FLAG_TYPE'] = 'OR'

            
        # convert cmd_dict to list
        cmd_list = list(itertools.chain.from_iterable(list(cmd_dict.items())))

        # log cmd executed
        cmd_str = ' '.join(cmd_list)
        log.info('SExtractor command executed:\n{}'.format(cmd_str))
        
        # run command
        result = subprocess.run(cmd_list, capture_output=True)
        status = result.returncode
        log.info('stdout: {}'.format(result.stdout.decode('UTF-8')))
        log.info('stderr: {}'.format(result.stderr.decode('UTF-8')))
        log.info('status: {}'.format(status))

        if status != 0:
            msg = 'SExtractor failed with exit code {}'.format(status)
            log.exception(msg)
            raise Exception(msg)

        if get_par(set_zogy.timing,tel):
            log_timing_memory (t0=t, label='run_sextractor before get_back')


        # improve background and its standard deviation estimate if
        # [set_zogy.bkg_method] not set to 1 (= use background
        # determined by SExtractor) and background had not already
        # been subtracted
        if (not return_fwhm_elong and
            npass==0 and
            get_par(set_zogy.bkg_method,tel)==2 and
            not bkg_sub and
            Scorr_mode is None):

            # read in input image; do not read header as that would
            # overwrite the header updates already done in this
            # function
            data = read_hdulist (image, dtype='float32')

            if fits_mask is not None:
                data_mask = read_hdulist (fits_mask, dtype='uint8')
            
            # construct background image using [get_back]; for ML/BG:
            # if the global background was already subtracted,
            # [data_bkg_mini] will be zero
            data_bkg_mini, data_bkg_std_mini = get_back (
                data, header, fits_objmask, fits_mask=fits_mask,
                tel=tel, set_zogy=set_zogy, imtype=imtype)


            # write these filtered meshes to fits and update their
            # headers with [set_zogy.bkg_boxsize]
            fits_tmp = '{}_bkg_mini.fits'.format(base)
            fits.writeto(fits_tmp, data_bkg_mini, overwrite=True)
            bkg_size = get_par(set_zogy.bkg_boxsize,tel)
            txt_tmp = '[pix] background boxsize used to create this image'
            fits.setval(fits_tmp, 'BKG-SIZE', value=bkg_size, comment=txt_tmp)


            # bkg STD
            fits_tmp = '{}_bkg_std_mini.fits'.format(base)
            fits.writeto(fits_tmp, data_bkg_std_mini, overwrite=True)
            fits.setval(fits_tmp, 'BKG-SIZE', value=bkg_size, comment=txt_tmp)


            # now use function [mini2back] to turn filtered mesh of median
            # and std of background regions into full background image and
            # its standard deviation
            bkg_size = get_par(set_zogy.bkg_boxsize,tel)
            data_bkg = mini2back (data_bkg_mini, data.shape, order_interp=2,
                                  bkg_boxsize=bkg_size,
                                  interp_Xchan=interp_Xchan,
                                  timing=get_par(set_zogy.timing,tel))
            data_bkg_std = mini2back (data_bkg_std_mini, data.shape,
                                      order_interp=1, bkg_boxsize=bkg_size,
                                      interp_Xchan=interp_Xchan_std,
                                      timing=get_par(set_zogy.timing,tel))


            # subtract the global background
            data -= data_bkg
            header['BKG-SUB'] = (True, 'sky background was subtracted?')
                
            header['S-BKG'] = (gain * np.median(data_bkg_mini), '[e-] median '
                               'background full image')
            header['S-BKGSTD'] = (gain * np.median(data_bkg_std_mini), '[e-] '
                                  'sigma (STD) background full image')


            # best to ensure that edge pixels are set to zero
            if fits_mask is not None:
                value_edge = get_par(set_zogy.mask_value['edge'],tel)
                mask_edge = (data_mask & value_edge == value_edge)
                data[mask_edge] = 0

            # save to fits; this image will be used to feed to
            # SExtractor in the next pass
            fits.writeto(image, data, header, overwrite=True)

            # write the improved background and standard deviation to fits
            # overwriting the fits images produced by SExtractor
            fits.writeto(fits_bkg, data_bkg, header, overwrite=True)
            fits.writeto(fits_bkg_std, data_bkg_std, header, overwrite=True)



    # now that catalog has been created, rename some of the columns
    # using the function rename_catcols. N.B.: except for when
    # run_sextractor is being used for the transient extraction
    # (Scorr_mode not None) this renaming was moved to the function
    # [ldac2fits] after converting the LDAC to normal fits catalog, to
    # avoid opening the potentially very large LDAC table just to
    # change some column names
    if Scorr_mode is not None:
        rename_catcols(cat_out)


    if return_fwhm_elong:
        # get estimate of seeing and elongation from output catalog
        fwhm, fwhm_std, elong, elong_std = get_fwhm(
            cat_out, get_par(set_zogy.fwhm_frac,tel),
            class_sort=get_par(set_zogy.fwhm_class_sort,tel), get_elong=True)

    else:
        fwhm = 0.
        fwhm_std = 0.
        elong = 0.
        elong_std = 0.



    # add number of objects detected (=number of catalog rows) to header
    header_catout = read_hdulist (cat_out, get_data=False, get_header=True)
    nobjects = header_catout['NAXIS2']

    header['S-NOBJ'] = (nobjects, 'number of objects detected by SExtractor')
    log.info('number of objects detected by SExtractor: {}'.format(nobjects))


    # also add header keyword(s) regarding background in case
    # background was determined by source-extractor
    if (get_par(set_zogy.bkg_method,tel)==1 and not return_fwhm_elong):
        data_bkg = read_hdulist (fits_bkg)
        header['S-BKG'] = (gain * np.median(data_bkg), '[e-] median background '
                           'full image')
        data_bkg_std = read_hdulist (fits_bkg_std)
        header['S-BKGSTD'] = (gain * np.median(data_bkg_std), '[e-] sigma (STD) '
                              'background full image')


    # remove file(s) if not keeping intermediate/temporary files
    if not get_par(set_zogy.keep_tmp,tel):

        if 'image_fraction' in locals():
            remove_files ([image_fraction], verbose=True)

        #if 'fits_objmask' in locals():
        #    remove_files ([fits_objmask])

        if fraction < 1:
            remove_files ([cat_out], verbose=True)



    if get_par(set_zogy.timing,tel):
        log_timing_memory (t0=t, label='run_sextractor')

    return fwhm, fwhm_std, elong, elong_std

def run_wcs (image_in, ra, dec, pixscale, width, height, header, imtype):

    if get_par(set_zogy.timing,tel): t = time.time()
    log.info('executing run_wcs ...')
    
    varyfrac = get_par(set_zogy.pixscale_varyfrac,tel)
    scale_low = (1 - varyfrac) * pixscale
    scale_high = (1 + varyfrac) * pixscale

    base = image_in.replace('.fits','')
    sexcat = '{}_cat.fits'.format(base)

    # read SExtractor catalogue (this is also used further down below
    # in this function)
    data_sexcat = read_hdulist (sexcat)
    #nobjects = data_sexcat.shape[0]
    #header['S-NOBJ'] = (nobjects, 'number of objects detected by SExtractor')
    #log.info('number of objects detected by SExtractor: {}'.format(nobjects))
    
    # feed Astrometry.net only with brightest sources; N.B.:
    # keeping only objects with zero FLAGS does not work well in crowded fields
    # select stars for finding WCS solution
    #mask_use = (data_sexcat['FLAGS']<=1)
    # the above selection is not working well for crowded images; Danielle
    # found the following to work well for both non-crowded as crowded
    # images, i.e. use all objects that are not masked according to the
    # input mask (any pixel within the isophotal area of an object):

    if 'FLAGS_MASK' in data_sexcat.dtype.names:
        mask_use = (data_sexcat['FLAGS_MASK']==0)
    else:
        mask_use = (data_sexcat['FLAGS']<=3)

    # sort in brightness (E_FLUX_AUTO)
    if 'E_FLUX_AUTO' in data_sexcat.dtype.names:
        column_sort = 'E_FLUX_AUTO'
    elif 'E_FLUX_OPT' in data_sexcat.dtype.names:
        column_sort = 'E_FLUX_OPT'
    else:
        column_sort = 'E_FLUX_APER_R5xFWHM'

    index_sort = np.argsort(data_sexcat[column_sort][mask_use])

    # select the brightest objects
    nbright = get_par(set_zogy.ast_nbright,tel)
    sexcat_bright = '{}_cat_bright.fits'.format(base)
    #fits.writeto(sexcat_bright, data_sexcat[:][mask_use][index_sort][-nbright:],
    #             overwrite=True)
    fits.writeto(sexcat_bright, data_sexcat[mask_use][index_sort][-nbright:],
                 overwrite=True)

    # create ds9 regions text file to show the brightest stars
    if get_par(set_zogy.make_plots,tel):
        result = prep_ds9regions(
            '{}_cat_bright_ds9regions.txt'.format(base),
            data_sexcat['X_POS'][mask_use][index_sort][-nbright:],
            data_sexcat['Y_POS'][mask_use][index_sort][-nbright:],
            radius=5., width=2, color='green',
            value=np.arange(1,nbright+1))

    dir_out = '.'
    if '/' in base:
        dir_out = '/'.join(base.split('/')[:-1])

    dstep = 30
    depth_str = (str(list(range(dstep,nbright//2+1,dstep)))
                 .replace('[','').replace(']','').replace(' ',''))

    cmd = ['solve-field', '--no-plots',
           #'--no-fits2fits', cloud version of astrometry does not have this arg
           '--x-column', 'X_POS',
           '--y-column', 'Y_POS',
           '--sort-column', column_sort,
           '--no-remove-lines',
           '--uniformize', '0',
           # only work on brightest sources
           #'--objs', str(nbright),
           '--width', str(width),
           '--height', str(height),           
           #'--keep-xylist', sexcat,
           # ignore existing WCS headers in FITS input images
           #'--no-verify', 
           #'--verbose',
           #'--verbose',
           '--crpix-center',
           #'--code-tolerance', str(0.01), 
           #'--quad-size-min', str(0.1),
           # for KMTNet images restrict the max quad size:
           #'--quad-size-max', str(0.1),
           # number of field objects to look at:
           '--depth', depth_str,
           #'--scamp', scampcat,
           # give up solving after the specified number of seconds of CPU time
           '--cpulimit', '120',
           sexcat_bright,
           '--tweak-order', str(get_par(set_zogy.astronet_tweak_order,tel)),
           '--scale-low', str(scale_low),
           '--scale-high', str(scale_high),
           '--scale-units', 'app',
           '--ra', str(ra),
           '--dec', str(dec),
           '--radius', str(get_par(set_zogy.astronet_radius,tel)),
           '--new-fits', 'none',
           '--overwrite',
           '--out', base.split('/')[-1],
           '--dir', dir_out
    ]

    if tel=='ML1':
        cmd += ['--parity', 'neg']


    # log cmd executed
    cmd_str = ' '.join(cmd)
    log.info('Astrometry.net command executed:\n{}'.format(cmd_str))

    result = subprocess.run(cmd, capture_output=True)
    status = result.returncode
    log.info('stdout: {}'.format(result.stdout.decode('UTF-8')))
    log.info('stderr: {}'.format(result.stderr.decode('UTF-8')))
    log.info('status: {}'.format(status))


    if os.path.exists('{}.solved'.format(base)) and status==0:

        # read .match file, which describes the quad match that solved the
        # image, before it is deleted
        data_match = read_hdulist ('{}.match'.format(base))

    else:
        msg = ('solve-field (Astrometry.net) failed with exit code {}'
               .format(status))
        log.exception(msg)
        raise Exception(msg)

    
    if get_par(set_zogy.timing,tel): t2 = time.time()


    # read header saved in .wcs 
    wcsfile = '{}.wcs'.format(base)
    header_wcs = read_hdulist (wcsfile, get_data=False, get_header=True)

    # remove HISTORY, COMMENT and DATE fields from Astrometry.net header
    # they are still present in the base+'.wcs' file
    header_wcs.pop('HISTORY', None)
    header_wcs.pop('COMMENT', None)
    header_wcs.pop('DATE', None)


    # add specific keyword indicating index file of match
    if data_match['HEALPIX'][0]!=-1:
        anet_index = 'index-{}-{:02d}.fits'.format(data_match['INDEXID'][0],
                                                   data_match['HEALPIX'][0])
    else:
        anet_index = 'index-{}.fits'.format(data_match['INDEXID'][0])

    header_wcs['A-INDEX'] = (anet_index, 'name of index file WCS solution')

    # and pixelscale
    cd1_1 = header_wcs['CD1_1']  # CD1_1 = CDELT1 *  cos (CROTA2)
    cd1_2 = header_wcs['CD1_2']  # CD1_2 = CDELT2 *  sin (CROTA2)
    cd2_1 = header_wcs['CD2_1']  # CD2_1 = CDELT1 * -sin (CROTA2)
    cd2_2 = header_wcs['CD2_2']  # CD2_2 = CDELT2 *  cos (CROTA2)


    anet_pixscale_x = np.sqrt(cd1_1**2 + cd2_1**2) * 3600.
    anet_pixscale_y = np.sqrt(cd1_2**2 + cd2_2**2) * 3600.
    anet_pixscale = np.average([anet_pixscale_x, anet_pixscale_y])


    # and rotation with the angle between the North and the second
    # axis (y-axis) of the image, counted positive to the East
    anet_rot_x = np.arctan2(-cd2_1, cd1_1) * (180./np.pi)
    anet_rot_y = np.arctan2( cd1_2, cd2_2) * (180./np.pi)
    anet_rot = np.average([anet_rot_x, anet_rot_y])

    # add header keywords
    header_wcs['A-PSCALE'] = (anet_pixscale,
                              '[arcsec/pix] pixel scale WCS solution')
    header_wcs['A-PSCALX'] = (anet_pixscale_x,
                              '[arcsec/pix] X-axis pixel scale WCS solution')
    header_wcs['A-PSCALY'] = (anet_pixscale_y,
                              '[arcsec/pix] Y-axis pixel scale WCS solution')

    header_wcs['A-ROT'] = (anet_rot,
                           '[deg] rotation WCS solution (E of N for "up")')
    header_wcs['A-ROTX'] = (anet_rot_x,
                            '[deg] X-axis rotation WCS (E of N for "up")')
    header_wcs['A-ROTY'] = (anet_rot_y,
                            '[deg] Y-axis rotation WCS (E of N for "up")')

    # convert SIP header keywords from Astrometry.net to PV keywords
    # that swarp, scamp (and sextractor) understand using this module
    # from David Shupe: sip_to_pv

    # using the old version of sip_to_pv (before June 2017):
    #status = sip_to_pv(image_in, image_in, tpv_format=True)
    #if status == False:
    #    log.error('sip_to_pv failed.')
    #    return 'error'

    # new version (June 2017) of sip_to_pv works on image header
    # rather than header+image (see below); the header is modified in
    # place; compared to the old version this saves an image write
    result = sip_to_pv(header_wcs, tpv_format=True, preserve=False)

    # update input header with [header_wcs]
    #header += header_wcs
    for key in header_wcs:
        header[key] = (header_wcs[key], header_wcs.comments[key])


    # use astropy.WCS to find RA, DEC corresponding to X_POS,
    # Y_POS, based on WCS info saved by Astrometry.net in .wcs
    # file (wcsfile). The 3rd parameter to wcs.all_pix2world indicates
    # the pixel coordinate of the frame origin. Using astropy.WCS
    # avoids having to save the new RAs and DECs to file and read them
    # back into python arrays. It provides the same RA and DEC as
    # wcs-xy2rd and also as SExtractor run independently on the WCS-ed
    # image.
    # N.B.: WCS accepts header objects - gets rid of old warning about
    # axis mismatch, as .wcs files have NAXIS=0, while proper image
    # header files have NAXIS=2
    wcs = WCS(header)
    newra, newdec = wcs.all_pix2world(data_sexcat['X_POS'],
                                      data_sexcat['Y_POS'],
                                      1)

    # update catalog with new RA and DEC columns
    data_sexcat['RA'] = newra
    data_sexcat['DEC'] = newdec
    fits.writeto(sexcat, data_sexcat, overwrite=True)
    
    if get_par(set_zogy.timing,tel):
        t3 = time.time()

    # check how well the WCS solution just found, compares with an
    # external catalog defined in Constants module
    if os.path.isfile(get_par(set_zogy.cal_cat,tel)):

        # use .wcs file to get RA, DEC of central pixel
        xsize = width
        ysize = height
        ra_center, dec_center = wcs.all_pix2world(xsize/2+0.5, ysize/2+0.5, 1)
        log.info('ra_center: {}, dec_center: {}'.format(ra_center, dec_center))

        # determine cal_cat min and max declination zone of field
        # to determine fits extensions (=zone+1) to read
        fov_half_deg = np.amax([xsize, ysize]) * pixscale / 3600. / 2
        zone_indices = get_zone_indices (dec_center, fov_half_deg, zone_size=60.)
        log.info('declination zone indices read from calibration catalog: {}'
                 .format(zone_indices))

        # read specific extensions (=zone_indices+1) of calibration catalog
        data_cal = read_hdulist (get_par(set_zogy.cal_cat,tel),
                                 ext_name_indices=zone_indices+1)

        # use function [find_stars] to select stars in calibration
        # catalog that are within the current field-of-view
        index_field = find_stars (data_cal['ra'], data_cal['dec'], ra_center,
                                  dec_center, fov_half_deg)
        #index_field = np.where(mask_field)[0]
        # N.B.: this [data_cal] array is returned by this function
        # [run_wcs] and also by [sex_wcs] so that it can be re-used
        # for the photometric calibration in [prep_optimal_subtraction]
        data_cal = data_cal[index_field]
        ra_ast = data_cal['ra']
        dec_ast = data_cal['dec']
        mag_ast = data_cal[get_par(set_zogy.ast_filter,tel)]

        n_aststars = np.shape(index_field)[0]
        log.info('number of potential astrometric stars in FOV: {}'
                 .format(n_aststars))

        # add header keyword(s):
        cal_name = get_par(set_zogy.cal_cat,tel).split('/')[-1]
        header['A-CAT-F'] = (cal_name, 'astrometric catalog') 
        header['A-TNAST'] = (n_aststars,
                             'total number of astrometric stars in FOV')

        # Limit to brightest stars ([nbright] is defined above) in the field
        #index_sort_ast = np.argsort(mag_ast)
        #ra_ast_bright = ra_ast[index_sort_ast][0:nbright]
        #dec_ast_bright = dec_ast[index_sort_ast][0:nbright]

        # calculate array of offsets between astrometry comparison
        # stars and any non-saturated SExtractor source
        newra_bright = newra[mask_use][index_sort][-nbright:]
        newdec_bright = newdec[mask_use][index_sort][-nbright:]        
        __, indx_ast, __, dra_array, ddec_array = get_matches (
            newra_bright, newdec_bright, ra_ast, dec_ast, dist_max=2)


        n_aststars_used = np.shape(dra_array)[0]
        log.info('number of astrometric stars used: {}'.format(n_aststars_used))
        header['A-NAST'] = (n_aststars_used,
                            'number of brightest stars used for WCS')
        header['A-NAMAX'] = (get_par(set_zogy.ast_nbright,tel),
                             'input max. number of stars to use for WCS')

        # calculate means, stds and medians
        dra_mean, dra_median, dra_std = sigma_clipped_stats(
            dra_array, sigma=5, mask_value=0)
        ddec_mean, ddec_median, ddec_std = sigma_clipped_stats(
            ddec_array, sigma=5, mask_value=0)
        
        log.info('dra_mean [arcsec]: {:.3f}, dra_std: {:.3f}, dra_median: {:.3f}'
                 .format(dra_mean, dra_std, dra_median))
        log.info('ddec_mean [arcsec]: {:.3f}, ddec_std: {:.3f}, ddec_median: '
                 '{:.3f}'.format(ddec_mean, ddec_std, ddec_median))
        
        # add header keyword(s):
        header['A-DRA'] = (dra_median,
                           '[arcsec] dRA median offset to astrom. catalog')
        header['A-DRASTD'] = (dra_std, '[arcsec] dRA sigma (STD) offset')
        header['A-DDEC'] = (ddec_median,
                            '[arcsec] dDEC median offset to astrom. catalog')
        header['A-DDESTD'] = (ddec_std, '[arcsec] dDEC sigma (STD) offset')

        if True or get_par(set_zogy.make_plots,tel):

            # plot of dra vs. ddec in arcseconds, including histograms
            dr = np.sqrt(dra_std**2+ddec_std**2)
            dr = max(dr, 0.05)
            limits = 5. * np.array([-dr,dr,-dr,dr])
            mask_nonzero = ((dra_array!=0.) & (ddec_array!=0.))
            label1 = 'dRA={:.3f}$\pm${:.3f}"'.format(dra_median, dra_std)
            label2 = 'dDEC={:.3f}$\pm${:.3f}"'.format(ddec_median, ddec_std)
            result = plot_scatter_hist(
                dra_array[mask_nonzero], ddec_array[mask_nonzero], limits,
                xlabel='delta RA [arcsec]', ylabel='delta DEC [arcsec]',
                label=[label1,label2], labelpos=[(0.77,0.9),(0.77,0.85)],
                filename='{}_dRADEC.pdf'.format(base),
                title=base.split('/')[-1])
            
            # plot of dra vs. color and ddec vs. color, to check for
            # offset dependence on stellar color
            col = ['g', 'r']
            if (col[0] in data_cal.dtype.names and
                col[1] in data_cal.dtype.names):
                
                col_array = data_cal[col[0]][indx_ast]-data_cal[col[1]][indx_ast]
                
                fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, sharey=False)
                fig.subplots_adjust(hspace=0, wspace=0)
                fig.suptitle(base.split('/')[-1], fontsize=10)

                # dRA subplot
                ax0.plot(col_array, dra_array, 'o', color='tab:blue',
                         markersize=5, markeredgecolor='k')
                ax0.set_ylabel('delta RA [arcsec]')
                # 1st degree polynomial fit to determin slope
                p_ra, V_ra = np.polyfit(col_array, dra_array, 1, cov=True)
                label_ra = 'slope={:.3f}$\pm${:.3f}'.format(p_ra[0],
                                                            np.sqrt(V_ra[0,0]))
                ax0.annotate(label_ra, xy=(1,1), xycoords='axes fraction',
                             fontsize=10, ha='right', va='top')
                # plot fit
                xp = np.linspace(np.amin(col_array), np.amax(col_array))
                ax0.plot (xp, xp*p_ra[0]+p_ra[1], color='k')

                # dDEC subplot
                ax1.plot(col_array, ddec_array, 'o', color='tab:blue',
                         markersize=5, markeredgecolor='k')
                ax1.set_ylabel('delta DEC [arcsec]')
                ax1.set_xlabel('{} - {} colour'.format(col[0], col[1]))

                # 1st degree polynomial fit to determin slope
                p_dec, V_dec = np.polyfit(col_array, ddec_array, 1, cov=True)
                label_dec = 'slope={:.3f}$\pm${:.3f}'.format(p_dec[0],
                                                             np.sqrt(V_dec[0,0]))
                ax1.annotate(label_dec, xy=(1,1), xycoords='axes fraction',
                             fontsize=10, ha='right', va='top')
                # plot fit
                ax1.plot (xp, xp*p_dec[0]+p_dec[1], color='k')


                plt.savefig('{}_dRADEC_vs_color.pdf'.format(base))
                if get_par(set_zogy.show_plots,tel): plt.show()
                plt.close()

                # also save fits file with calibration info and
                # offsets for ADC analysis
                table_cal_offsets = Table(data_cal[indx_ast])
                table_cal_offsets['DRA'] = dra_array
                table_cal_offsets['DDEC'] = ddec_array
                table_cal_offsets.write ('{}_dRADEC_calcat.fits'.format(base),
                                         format='fits', overwrite=True)


    else:
        log.info('Warning: calibration catalog {} not found!'
                 .format(get_par(set_zogy.cal_cat,tel)))
        data_cal = None

        
    # add RA-CNTR and DEC-CNTR to header
    wcs = WCS(header)
    ra_center, dec_center = wcs.all_pix2world(width/2+0.5, height/2+0.5, 1)
    log.info('ra_center: {}, dec_center: {}'.format(ra_center, dec_center))
    header['RA-CNTR'] = (float(ra_center), 
                         'RA (ICRS) at image center (astrometry.net)')
    header['DEC-CNTR'] = (float(dec_center), 
                          'DEC (ICRS) at image center (astrometry.net)')


    # write data_cal with selection of calibration stars in this field
    # to fits table; one for each imtype as the shift may be
    # considerable
    if data_cal is not None:
        fits_calcat_field = '{}_calcat_field_{}.fits'.format(base, imtype)
        #data_cal.writeto (fits_calcat_field, overwrite=True)
        table_calcat_field = Table(data_cal)
        table_calcat_field.write (fits_calcat_field, format='fits',
                                  overwrite=True)


    # update header of input image
    header['DATEFILE'] = (Time.now().isot, 'UTC date of writing file')
    with fits.open(image_in, 'update', memmap=True) as hdulist:
        hdulist[0].header = header


    # remove file(s) if not keeping intermediate/temporary files
    if not get_par(set_zogy.keep_tmp,tel):

        if 'sexcat_bright' in locals():
            remove_files ([sexcat_bright], verbose=True)

        # astrometry.net output files
        list2remove = ['{}{}'.format(base, ext) for ext in
                       ['.solved', '.match', '.rdls', '.corr', '-indx.xyls',
                        '.axy', '.wcs']]
        remove_files (list2remove, verbose=True)



    if get_par(set_zogy.timing,tel):
        log_timing_memory (t0=t3, label='calculate offset wrt external catalog')
        
    if get_par(set_zogy.timing,tel):
        log_timing_memory (t0=t, label='run_wcs')

    return




from Funcoes import FuncoesReducao
import os
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import numpy as np
import time
import astropy.io.fits as fits

directory = os.getcwd()
imagem = '/Files/Imagens/2022-03-21/ImagensReduzidas/1037-293R_001.fits'
imagem2 = '/home/luis/Desktop/14mai28/1037-293R_002.fits'
raw = CCDData.read(directory + imagem, unit = 'adu') 





resultado = FuncoesReducao.sextractor(directory+imagem,directory,'test.cat')
print(resultado)
SexData = FuncoesReducao.ReadSextractorFile('test.cat')


# ~ result,status = FuncoesReducao.astrometry(imagem2,raw.header,directory + '/Astrometry','X','Y')

# ~ print(result)


# ~ print(FuncoesReducao.SeeFitsTables(directory + '/Astrometry/1037-293R_002.axy'))














vmin,vmax = None,None
plt.imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
for i in range(len(SexData['X_IMAGE'])):
	x,y,r = float(SexData['X_IMAGE'][i]),float(SexData['Y_IMAGE'][i]),float(SexData['FLUX_ISO'][i])/2000
	plt.scatter(x, y, r, marker='o', facecolors='none',color='Red')
plt.show()








# ~ if config["DARK"]:
	# ~ if os.path.isfile(diretoria + "/Master/MasterDark" + filtro + ".fits"):
		# ~ MasterDark = CCDData.read(diretoria + "/Master/MasterDark" + filtro + ".fits")
	# ~ elif os.path.isdir(diretoria + "/Master/MasterDark" + filtro + ".fits")
	
	
	
	
	# ~ if len(os.listdir(path=diretoria + "/Darks")) == 0:
		# ~ logging.warning("Não existem Darks, removendo este passo da reducao")
		# ~ config["DARK"] = False
	# ~ else:
		# ~ logging.debug("Criando os MasterDarks")
		# ~ for filtro in os.listdir(path=diretoria + "/Darks"):
			# ~ if not os.path.isfile(diretoria + "/Master/MasterDark" + filtro + ".fits"):
				# ~ logging.debug("Criando o MasterDark do filtro " + filtro)
				# ~ MasterDark = FuncoesReducao.CriarMaster(diretoria + "/Darks/" + filtro + '/', diretoria, filtro,'dark')
				# ~ logging.debug("MasterDark Criado do filtro " + filtro)
			# ~ else:
				# ~ MasterDark = CCDData.read(diretoria + "/Master/MasterDark" + filtro + ".fits") 
		# ~ logging.debug("MasterDarks Criados")
																
# ~ if config["FLAT"]:
	# ~ if len(os.listdir(path=diretoria + "/Flats")) == 0:
		# ~ logging.warning("Não existem Flats, removendo este passo da reducao")
		# ~ config["FLAT"] = False
	# ~ else:		
		# ~ logging.debug("Criando os MasterFlats")
		# ~ for filtro in os.listdir(path=diretoria + "/Flats"):		
			# ~ if not os.path.isfile(diretoria + "/Master/MasterFlat" + filtro + ".fits"):
				# ~ logging.debug("Criando o MasterFlat do filtro " + filtro)
				# ~ MasterFlat = FuncoesReducao.CriarMaster(diretoria + "/Flats/" + filtro + '/', diretoria, filtro,'flat')
				# ~ logging.debug("MasterFlat Criado do filtro " + filtro)
			# ~ else:
				# ~ MasterFlat = CCDData.read(diretoria + "/Master/MasterFlat" + filtro + ".fits") 
		# ~ logging.debug("MasterFlats Criados")

