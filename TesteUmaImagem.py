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
MagLim = 15





Imagens = ['MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fits']

# ~ catalog_path = directory + '/Files/Catalogos/tycho2.dat'
# ~ catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)

fig,axes = plt.subplots(2,3)
ax = axes.ravel()
n = 0
for ficheiro in Imagens:
	DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		
	nomeficheiro = ficheiro.split('/')[-1].split('.')[0]
	Data = CCDData.read(DirImage, unit = 'adu') 

	#Sextractor
	SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
	resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_LDAC',int(config['VIGNET']))  
	SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'
	resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_1.0')  
	FuncoesReducao.RemoveBorders(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat',2,int(config['VIGNET']))
	FuncoesReducao.RemoveBorders(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat',1,int(config['VIGNET']))	
			
	#Astrometry
	result,status = FuncoesReducao.astrometry(diretoria + '/SextractorCatalogue/' + ficheiro[:-5] + '.cat',Data.header,diretoria + '/AstrometryFiles',config) 

	# ~ #Transforms from XY to RaDec with WCS information from astrometry
	nomeficheiro = ficheiro.split('/')[-1].split('.')[0]														# Saves name of file (removes rest of extension)
	FileWc = diretoria + '/AstrometryFiles/' + nomeficheiro + '.wcs'											# Gets directory : WCS file created from astrometry 
	FileXY = diretoria + '/SextractorCatalogue/' + nomeficheiro + '.cat'										# Gets directory : catalogue created from sextractor
	FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro													# Output file for calculation of RaDec coordinate of all stars
	result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)														# Calculates RaDec positions for all stars obtained from sextractor


	#Gets the information from the sextractor file and the calculated positions
	FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro	
	Posicoes = fits.open(FileOutput)												# Gets positions calculated
	Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
	SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'							# Saves RaDec positions
	SexData = FuncoesReducao.ReadSextractorFileFit(SexFile)			# Reads sextractor file obtained in the beggining and saves desired values
	X,Y,Mags,MagsErro,BackG,fluxo,FWHM,SNR = SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE'],SexData['SNR_WIN']
	
	

	#Matches the objects found and the catalogue
	stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
	MatchTime = time.time()
	FileWc = diretoria + '/AstrometryFiles/' + nomeficheiro + '.wcs'											# Gets directory : WCS file created from astrometry 
	WCS = fits.open(FileWc)
	scale = float(WCS[0].header['COMMENT'][-8].split()[1])
	max_sep = 0.8*float(config['max_sep'])*np.median(FWHM[FWHM!=0])*scale/3600

	catalog_path = directory + '/Files/Catalogos/UCAC5/'
	catalogo,mag = FuncoesReducao.ReadUCAC5(FuncoesReducao.ZonesUCAC(min(Dec),max(Dec)),catalog_path,MagLim)
	idx,found,dx = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars
	
	#Run PSFEX
	Base = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5]
	CatIn = Base +'LDAC.cat'
	CatOut = Base + '_out.cat'
	ConfigFile =  directory + '/Files/Sextractor/defaultBRUNO.psfex' 
	FuncoesReducao.Rem0FWHM(CatIn)
	result, status = FuncoesReducao.psfex(CatIn,CatOut,ConfigFile)
	print(result)
	
	psf_in = CatIn.replace('.cat', '.psf')
	psf_out = '{}_psf.fits'.format(Base)
	os.rename (psf_in, psf_out)
	FuncoesReducao.SeePSF(psf_out)
	
	


	cores = ['r','g','c','y']
	Imagem = CCDData.read(diretoria + "/ImagensReduzidas/" + ficheiro, unit = 'adu') 			# Reads Image
	ax[n].imshow(Imagem,cmap='gray')
	for contador in range(len(X)):
		if FWHM[contador]!= 0:
			c1 = plt.Circle((X[contador],Y[contador]),3*(n+1), color = cores[n], linewidth=2,fill=False)
			ax[n].add_patch(c1)
			c1 = plt.Circle((X[contador],Y[contador]),3*(n+1), color = cores[n], linewidth=2,fill=False)
			ax[4].add_patch(c1)
	ax[n].set_title('Number of detections: '+str(len(Y)))
	n +=1
	

	# ~ #Desenho
	# ~ A = fits.open(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat')
	# ~ Data = A[2].data
	# ~ X1 = X[found]
	# ~ Y1 = Y[found]
	# ~ X2 = X[~found]
	# ~ Y2 = Y[~found]
	# ~ Imagem = CCDData.read(diretoria + "/ImagensReduzidas/" + ficheiro, unit = 'adu') 			# Reads Image

	# ~ fig,axes = plt.subplots(1,3)
	# ~ ax = axes.ravel()
	# ~ ax[0].imshow(Imagem,cmap='gray')
	# ~ ax[0].set_title('Objects Detected')
	# ~ ax[1].imshow(Imagem,cmap='gray')
	# ~ ax[1].set_title('Objects found in catalog')
	# ~ ax[2].imshow(Imagem,cmap='gray')
	# ~ ax[2].set_title('Objects not found in catalog')
	# ~ for contador in range(len(Y)):
		# ~ c1 = plt.Circle((X[contador],Y[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[0].add_patch(c1)

	# ~ for contador in range(len(Y1)):
		# ~ c1 = plt.Circle((X1[contador],Y1[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[1].add_patch(c1)
		
	# ~ for contador in range(len(Y2)):
		# ~ c1 = plt.Circle((X2[contador],Y2[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[2].add_patch(c1)
				
	# ~ plt.show()

	# ~ fig,axes = plt.subplots(2,3)
	# ~ ax = axes.ravel()
	# ~ ax[0].hist(dx[found].value,bins = 20)
	# ~ ax[1].hist(dx[~found].value,bins = 20)
	# ~ ax[2].hist(dx.value, bins = 40)

	# ~ ax[0].set_title('Encontradas')
	# ~ ax[1].set_title('Não Encontrados')
	# ~ ax[2].set_title('Total')

	# ~ ax[0].set_ylabel('Distâncias (graus)')
	# ~ ax[3].set_ylabel('Magnitudes (instrumental)')


	# ~ ax[3].hist(Mags[found],bins = 20)
	# ~ ax[4].hist(Mags[~found], bins = 20)
	# ~ ax[5].hist(Mags, bins = 40)

	# ~ plt.show()
ax[4].imshow(Imagem,cmap='gray')
plt.show()

