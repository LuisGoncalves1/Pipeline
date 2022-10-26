from Funcoes import FuncoesReducao, FuncoesMaster
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
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationBruno.txt")
config['LOGGER'] = int(config['LOGGER'])
data = "2022-03-20"
diretoria = "%s/Files/Imagens/%s" %(directory,str(data))



ficheiro = '1037-293B_001.fits'
# ~ ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits'
DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		

SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_LDAC',VIGNET=25)  
SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'
resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_1.0')  
print('ola')

#Runs psfex
logging.debug("Starting PSFex on " + ficheiro)
PsfTime = time.time()
SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
PsfFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'psfex.cat'
PsfConfig = directory + '/Files/Sextractor/' + config['PSFexConfig']
resultado,status = FuncoesReducao.psfex(SexFile,PsfConfig)  
print(resultado)
print(status)
PsfTime = time.time() - PsfTime
