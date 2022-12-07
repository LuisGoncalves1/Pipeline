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


ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits'
DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		


CatIn = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
ConfigFile =  directory + '/Files/Sextractor/defaultBRUNO.psfex' 

PSFFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] + '_psf.fits'
FuncoesReducao.SeePSF(PSFFile)
time.sleep(500)







FuncoesReducao.Rem0FWHM(CatIn)
# ~ result, status = FuncoesReducao.psfex(CatIn,ConfigFile)

Base = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5]
CatIn = Base +'LDAC.cat'
ConfigFile =  directory + '/Files/Sextractor/defaultBRUNO.psfex' 
psf_in = CatIn.replace('.cat', '.psf')
print('psf_in: ',psf_in)
psf_out = '{}_psf.fits'.format(Base)
print('psf_out: ',psf_out)
os.rename (psf_in, psf_out)










  
