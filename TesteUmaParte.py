from Funcoes import FuncoesReducao
import ccdproc 
from datetime import date
import os
from astropy.nddata import CCDData
import astropy.units as u
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
import numpy as np
import logging
import cv2
import time
import matplotlib.pyplot as plt
from pathlib import Path
import astroscrappy
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky




directory = os.getcwd()
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationDefesa.txt")
config['LOGGER'] = int(config['LOGGER'])
data = "2022-03-20"
diretoria = "%s/Files/Imagens/%s" %(directory,str(data))
MagLim = 14





# ~ ficheiro = '1037-293R_003.fits'
# ~ ficheiro = '1048-411R_004.fits'
ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fits'

DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		


nomeficheiro = ficheiro.split('/')[-1].split('.')[0]
#Gets the information from the sextractor file and the calculated positions
FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro	
Posicoes = fits.open(FileOutput)												# Gets positions calculated
Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'							# Saves RaDec positions
SexData = FuncoesReducao.ReadSextractorFileFit(SexFile)			# Reads sextractor file obtained in the beggining and saves desired values
X,Y,Mags,MagsErro,BackG,fluxo,FWHM= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE']
		
		

#Matches the objects found and the catalogue
stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
MatchTime = time.time()
FileWc = diretoria + '/AstrometryFiles/' + nomeficheiro + '.wcs'											# Gets directory : WCS file created from astrometry 
WCS = fits.open(FileWc)
scale = float(WCS[0].header['COMMENT'][-8].split()[1])
max_sep = float(config['max_sep'])*np.median(FWHM[FWHM!=0])*scale/3600
print(max_sep)

A = time.time()

catalog_path = directory + '/Files/Catalogos/UCAC4/'
catalogo= FuncoesReducao.ReadUCAC4(FuncoesReducao.ZonesUCAC4(min(Dec),max(Dec)),catalog_path,MagLim)
# ~ catalog_path = directory + '/Files/Catalogos/tycho2.dat'
# ~ catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)

tempoCatalogo = time.time()-A


idx,found,dx = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars
print(dx)
print(np.median(dx))
# ~ print(Dec[~found])
# ~ print(idx)
# ~ print(found)
# ~ print(dx)
# ~ print('------')
# ~ print(dx[found])
# ~ print(np.mean(dx[found]))
# ~ print(Mags[found])
# ~ print(np.mean(Mags[found]))
# ~ print('------')
# ~ print(dx[~found])
# ~ print(np.mean(dx[~found]))
# ~ print(Mags[~found])
# ~ print(np.mean(Mags[~found]))
MatchTime = time.time() - MatchTime
		






#Desenho
A = fits.open(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat')
Data = A[2].data
X1 = X[found]
Y1 = Y[found]
X2 = X[~found]
Y2 = Y[~found]
Imagem = CCDData.read(diretoria + "/ImagensReduzidas/" + ficheiro, unit = 'adu') 			# Reads Image

fig,axes = plt.subplots(2,2)
ax = axes.ravel()
ax[0].imshow(Imagem,cmap='gray')
ax[0].set_title('Todas')
ax[1].imshow(Imagem,cmap='gray')
ax[1].set_title('Encontradas')
ax[2].imshow(Imagem,cmap='gray')
ax[2].set_title('Faltam')
ax[3].imshow(Imagem,cmap='gray')
fig.suptitle('Mag: ' + str(MagLim)+ '  Tempo: ' + str(tempoCatalogo)+  '  NEstrelas: ' + str(len(catalogo.dec)))
for contador in range(len(Y)):
	c1 = plt.Circle((X[contador],Y[contador]),4, color = 'green', linewidth=2,fill=False)
	ax[0].add_patch(c1)

for contador in range(len(Y1)):
	c1 = plt.Circle((X1[contador],Y1[contador]),4, color = 'green', linewidth=2,fill=False)
	ax[1].add_patch(c1)
	
for contador in range(len(Y2)):
	c1 = plt.Circle((X2[contador],Y2[contador]),4, color = 'green', linewidth=2,fill=False)
	ax[2].add_patch(c1)
			
plt.show()

fig,axes = plt.subplots(2,3)
ax = axes.ravel()
ax[0].hist(dx[found].value,bins = 20)
ax[1].hist(dx[~found].value,bins = 20)
ax[2].hist(dx.value, bins = 40)

ax[0].set_title('Encontradas')
ax[1].set_title('Não Encontrados')
ax[2].set_title('Total')

ax[0].set_ylabel('Distâncias (graus)')
ax[3].set_ylabel('Magnitudes (instrumental)')


ax[3].hist(Mags[found],bins = 20)
ax[4].hist(Mags[~found], bins = 20)
ax[5].hist(Mags, bins = 40)

plt.show()

















# ~ def run_psfex (cat_in, file_config, cat_out, imtype, poldeg, nsnap=8,
               # ~ limit_ldac=False, nthreads=0):
    
    # ~ """Function that runs PSFEx on [cat_in] (which is a SExtractor output
       # ~ catalog in FITS_LDAC format) using the configuration file
       # ~ [file_config]"""


    # ~ if get_par(set_zogy.timing,tel):
        # ~ t = time.time()

    
    # ~ if imtype=='new':
        # ~ base = base_new
    # ~ else:
        # ~ base = base_ref

        
    # ~ # select a subset of entries in the input ldac catalog to speed up
    # ~ # psfex
    # ~ if limit_ldac:

        # ~ mem_use ('before limiting LDAC')
        
        # ~ with fits.open(cat_in, mode='update') as hdulist:
            # ~ data_ldac = hdulist[-1].data
            # ~ # SExtractor flags 0 or 1 and S/N larger than value
            # ~ # defined in settings file; note that this is the same as
            # ~ # the default rejection mask on SExtractor FLAGS defined
            # ~ # in the psfex.config file: SAMPLE_FLAGMASK = 0x00fe
            # ~ # (=254) or '0b11111110'
            # ~ s2n = get_par(set_zogy.psf_stars_s2n_min,tel)
            # ~ mask_ok = ((data_ldac['FLAGS']<=1) & 
                       # ~ (data_ldac['SNR_WIN']>=s2n) &
                       # ~ (data_ldac['IMAFLAGS_ISO']==0))


            # ~ mem_use ('after mask_ok filter')
            
            # ~ mask_sum = np.sum(mask_ok)
            # ~ log.info ('number of PSF stars available with FLAGS<=1, SNR>{} '
                      # ~ 'and IMAFLAGS_ISO==0: {}'.format(s2n, mask_sum))

            # ~ # to limit the number of PSF stars to a reasonable number
            # ~ # in crowded fields, pick a random set of [nlimit] stars
            # ~ # if there are at least that number available
            # ~ nlimit = 20000
            # ~ if mask_sum > nlimit:
                # ~ index_keep = (np.random.rand(nlimit) * mask_sum).astype(int)
                # ~ log.info ('using random subset of {} of these'
                          # ~ .format(len(index_keep)))
            # ~ else:
                # ~ index_keep = np.arange(mask_sum)

            # ~ hdulist[2].data = data_ldac[mask_ok][index_keep]
            
            
        # ~ if get_par(set_zogy.make_plots,tel):
            # ~ result = prep_ds9regions('{}_psfstars_ds9regions.txt'.format(base),
                                     # ~ data_ldac['XWIN_IMAGE'][mask_ok],
                                     # ~ data_ldac['YWIN_IMAGE'][mask_ok],
                                     # ~ radius=5., width=2, color='red')
                
        # ~ if get_par(set_zogy.timing,tel):
            # ~ log_timing_memory (
                # ~ t0=t,label='limiting entries in LDAC input catalog for PSFEx')

        
    # ~ # use function [get_samp_PSF_config_size] to determine [psf_samp]
    # ~ # and [psf_size_config] required to run PSFEx
    # ~ psf_samp, psf_size_config = get_samp_PSF_config_size(imtype)
    # ~ psf_size_config_str = '{},{}'.format(psf_size_config, psf_size_config)
    
    # ~ if get_par(set_zogy.verbose,tel):
        # ~ log.info('psf_size_config: {} in [run_psfex] for imtype {}'
                 # ~ .format(psf_size_config, imtype))

    # ~ # get FWHM and ELONGATION to limit the PSFex configuration
    # ~ # parameters SAMPLE_FWHMRANGE and SAMPLE_MAXELLIP
    # ~ #fwhm, fwhm_std, elong, elong_std = get_fwhm(cat_in, 0.05, class_sort=False,
    # ~ #                                                      get_elong=True)
    # ~ #print ('fwhm, fwhm_std, elong, elong_std', fwhm, fwhm_std, elong, elong_std)
    # ~ #sample_fwhmrange = str(fwhm-fwhm_std)+','+str(fwhm+fwhm_std)
    # ~ #print ('sample_fwhmrange', sample_fwhmrange)
    # ~ #maxellip = (elong+3.*elong_std-1)/(elong+3.*elong_std+1)
    # ~ #maxellip_str= str(np.amin([maxellip, 1.]))
    # ~ #print ('maxellip_str', maxellip_str)

    # ~ # run psfex from the unix command line
    # ~ cmd = ['psfex', cat_in, '-c', file_config,'-OUTCAT_NAME', cat_out,
           # ~ '-PSF_SIZE', psf_size_config_str, '-PSF_SAMPLING', str(psf_samp),
           # ~ '-SAMPLE_MINSN', str(get_par(set_zogy.psf_stars_s2n_min,tel)),
           # ~ '-NTHREADS', str(nthreads), '-PSFVAR_NSNAP', str(nsnap),
           # ~ '-PSFVAR_DEGREES', str(poldeg)]
    # ~ #       '-SAMPLE_FWHMRANGE', sample_fwhmrange,
    # ~ #       '-SAMPLE_MAXELLIP', maxellip_str]

    # ~ if get_par(set_zogy.make_plots,tel):
        # ~ cmd += ['-CHECKPLOT_TYPE', 'FWHM, ELLIPTICITY, COUNTS, COUNT_FRACTION, '
                # ~ 'CHI2, RESIDUALS',
                # ~ '-CHECKPLOT_DEV', 'PS',
                # ~ '-CHECKPLOT_ANTIALIAS', 'N',
                # ~ '-CHECKPLOT_NAME',
                # ~ 'psfex_fwhm, psfex_ellip, psfex_counts, psfex_countfrac, '
                # ~ 'psfex_chi2, psfex_resi']
        # ~ cmd += ['-CHECKIMAGE_TYPE', 'CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS,'
                # ~ 'BASIS,MOFFAT,-MOFFAT,-SYMMETRICAL',
                # ~ '-CHECKIMAGE_NAME',
                # ~ 'psfex_chi, psfex_proto, psfex_samp, psfex_resi, psfex_snap, '
                # ~ 'psfex_basis, psfex_moffat, psfex_minmoffat, psfex_minsymmetrical']

    # ~ # log cmd executed
    # ~ cmd_str = ' '.join(cmd)
    # ~ log.info('PSFEx command executed:\n{}'.format(cmd_str))

    # ~ result = subprocess.run(cmd, capture_output=True)
    # ~ status = result.returncode
    # ~ log.info('stdout: {}'.format(result.stdout.decode('UTF-8')))
    # ~ log.info('stderr: {}'.format(result.stderr.decode('UTF-8')))
    # ~ log.info('status: {}'.format(status))

    # ~ if status != 0:
        # ~ msg = 'PSFEx failed with exit code {}'.format(status)
        # ~ log.exception(msg)
        # ~ raise Exception(msg)


    # ~ # standard output of PSFEx is .psf; change this to _psf.fits
    # ~ psf_in = cat_in.replace('.fits', '.psf')
    # ~ psf_out = '{}_psf.fits'.format(base)
    # ~ os.rename (psf_in, psf_out)

    # ~ # if zero stars were found by PSFEx, raise exception
    # ~ header_psf = read_hdulist(psf_out, get_data=False, get_header=True)
    # ~ if header_psf['ACCEPTED']==0:
        # ~ msg = ('no appropriate source found by PSFEx in [run_psfex] for catalog '
               # ~ '{}'.format(cat_in))
        # ~ log.exception (msg)
        # ~ raise Exception (msg)


    # ~ # no need to limit disk space of LDAC catalog anymore, as it is no
    # ~ # longer saved
    # ~ if False:
        # ~ # for the reference image, limit the size of the ldac catalog fits
        # ~ # file to be saved with only those catalog entries used by PSFEx
        # ~ # --> why not for new image as well?
        # ~ if imtype=='ref' or imtype=='new':
            # ~ psfexcat = '{}_psfex.cat'.format(base)
            # ~ table = ascii.read(psfexcat, format='sextractor')
            # ~ # In PSFEx version 3.18.2 all objects from the input
            # ~ # SExtractor catalog are recorded, and in that case the
            # ~ # entries with FLAGS_PSF=0 need to be selected.
            # ~ if 'FLAGS_PSF' in table.colnames:
                # ~ mask_psfstars = (table['FLAGS_PSF']==0)
            # ~ # In PSFEx version 3.17.1 (last stable version), only stars
            # ~ # with zero flags are recorded in the output catalog, so need
            # ~ # to pick the source number of the input catalog minus 1.
            # ~ else:
                # ~ mask_psfstars = table['SOURCE_NUMBER']-1

            # ~ # overwrite the input catalog with these selected stars
            # ~ with fits.open(cat_in, mode='update') as hdulist:
                # ~ hdulist[2].data = hdulist[2].data[mask_psfstars]


    # ~ # record the PSFEx output check images defined above, into
    # ~ # extensions of a single fits file    
    # ~ if get_par(set_zogy.make_plots,tel):

        # ~ hdulist = fits.HDUList()
        # ~ hdulist.append(fits.PrimaryHDU())
    
        # ~ cwd = os.getcwd()
        # ~ file_list = glob.glob('{}/psfex*'.format(cwd))
        # ~ for name in file_list:
        
            # ~ short = name.split('/')[-1]
            # ~ prefix = '_'.join(short.split('_')[0:2])

            # ~ if 'fits' in name:
                # ~ # read image data and header
                # ~ data, header = read_hdulist (name, get_header=True)
                # ~ # add extension to hdulist with extname=prefix
                # ~ hdulist.append(fits.ImageHDU(data=data, header=header, 
                                             # ~ name=prefix))
                # ~ # remove individual image
                # ~ os.remove(name)
                # ~ log.info ('added fits image {} to multi-extension fits file'
                          # ~ .format(name))
                
            # ~ else:
                # ~ # if not a fits file, rename it
                # ~ name_new = short.replace('ldac', prefix)
                # ~ name_new = name_new.replace('{}_'.format(prefix),'')
                # ~ name_new = '{}/{}'.format(cwd, name_new)
                # ~ log.info ('name: {}, name_new: {}'.format(name, name_new))
                # ~ os.rename (name, name_new)
                
        # ~ # write image
        # ~ fits_output = '{}_psf_checkimages.fits'.format(base)
        # ~ hdulist.writeto(fits_output, overwrite=True)


    # ~ # remove file(s) if not keeping intermediate/temporary files
    # ~ if not get_par(set_zogy.keep_tmp,tel):
        # ~ remove_files ([cat_in], verbose=True)


    # ~ if get_par(set_zogy.timing,tel):
        # ~ log_timing_memory (t0=t, label='run_psfex')

    # ~ return
