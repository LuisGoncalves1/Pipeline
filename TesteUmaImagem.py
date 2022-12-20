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
config = FuncoesReducao.ReadPipelineConfig(PipelineDirectory+"/Files/Config/ConfigurationDefesa.txt")
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day)) #diretoria
LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(DayDirectory,str(day))
logging.basicConfig(filename = log_file_name, level = int(config['LOGGER']),format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()
Masters = FuncoesReducao.SeeMasters(DayDirectory)




MagLim = 15
Imagens = DayDirectory + '/OriginalImages'
Imagens = os.listdir(Imagens)
# ~ Imagens = ['MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits',
# ~ 'MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fits',
# ~ 'MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fits',
# ~ 'MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fits']

print(Imagens)
n = 0
for ficheiro in Imagens:	
	Base = ficheiro.split('/')[-1].split('.')[0]
	FuncoesReducao.ImageProcess(Base,PipelineDirectory,DayDirectory,Masters,config,logging)
