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

PipelineDirectory = os.getcwd()
data = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(data))
config = FuncoesReducao.ReadPipelineConfig(PipelineDirectory+"/Files/Config/ConfigurationDefesa.txt")



Base = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01'

# ~ A = fits.open('/home/luis/Desktop/Codigos/Pipeline/Files/Images/2022-03-20/SextractorCatalogue/MAPA_SURV_4_085283_220417_114004+060023_scan047_01_psfLDAC.cat')
# ~ with open('batata.txt','w') as ficheiro:
	# ~ X = A[2].data['X_IMAGE'] 
	# ~ Y = A[2].data['Y_IMAGE']
	# ~ for i in range(len(X)):
		# ~ linha = f'{X[i]:5},{Y[i]:5}\n'
		# ~ ficheiro.write(linha)



FileWc = DayDirectory + '/AstrometryFiles/' + Base + '.wcs'															# Gets directory : WCS file created from astrometry 
FileXY = '/home/luis/Desktop/Codigos/Pipeline/tal.fits'


FileOutput = DayDirectory + '/AstrometryFiles/' + Base	+'2'															# Output file for calculation of RaDec coordinate of all stars
FileOutput = PipelineDirectory + '/teste'
result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)																# Calculates RaDec positions for all stars obtained from sextractor

stars, Mags, catalogo, mag, idx, found = FuncoesReducao.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)
print(mag)
print(Mags)
# ~ print(found)
