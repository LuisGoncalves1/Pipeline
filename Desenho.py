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
data = "2022-03-20"
diretoria = "%s/Files/Imagens/%s" %(directory,str(data))
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationDefesa.txt")

Imagens = ['MAPA_SURV_4_085283_220417_114004+060023_scan047_01',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_02',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_03',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_04']


catalog_path = directory + '/Files/Catalogos/tycho2.dat'
catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)











vmin = None
vmax = None
fig,axes = plt.subplots(2,3)
ax = axes.ravel()
for i in range(len(Imagens)):
	Imagem = CCDData.read(diretoria + "/ImagensReduzidas/" + Imagens[i] + '.fits', unit = 'adu') 				
	ax[i].imshow(Imagem,cmap='gray')
	
	A = fits.open(diretoria + '/SextractorCatalogue/'+ Imagens[i] +'.cat')
	Data = A[1].data
	X = Data['X_IMAGE']
	Y = Data['Y_IMAGE']
	FWHM = Data['FWHM_IMAGE']
	FWHM = FWHM[FWHM>0]
	
	FileOutput = diretoria + '/AstrometryFiles/' + Imagens[i]				# Output file for calculation of RaDec coordinate of all stars														
	Posicoes = fits.open(FileOutput)												# Gets positions calculated
	Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
	stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
	FileWc = diretoria + '/AstrometryFiles/' + Imagens[i] + '.wcs'
	
	WCS = fits.open(FileWc)
	scale = float(WCS[0].header['COMMENT'][-8].split()[1])
	max_sep = float(config['max_sep'])*np.median(FWHM)*scale/3600
	
	idx,found,dx = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars

	FileRDLS = diretoria + '/' + Imagens[i]
	RLDS = fits.open(FileRDLS)
	X2 = RLDS[1].data['X']
	Y2 = RLDS[1].data['Y']
	
	
	X = X[~found]
	Y = Y[~found]
	ax[i].set_title(str(max_sep*3600))
	for contador in range(len(Y)):
		color = 'green'
		# ~ if FWHM[contador] == 0:
			# ~ color = 'red'
		c1 = plt.Circle((X[contador],Y[contador]),4, color = color, linewidth=2,fill=False)
		ax[i].add_patch(c1)
	
	for contador in range(len(X2)):
		color = 'red'
		c1 = plt.Circle((X2[contador],Y2[contador]),2, color = color, linewidth=2,fill=False)
		ax[i].add_patch(c1)
plt.show()


