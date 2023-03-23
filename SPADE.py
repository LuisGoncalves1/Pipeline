from Functions import SPADEFunctions
from datetime import date
import os
from astropy.nddata import CCDData
import ccdproc 
import astropy.units as u
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
import numpy as np
import cv2
import logging
import time
import matplotlib.pyplot as plt
from pathlib import Path
import astroscrappy
import concurrent.futures
import psutil

import sys

"""
The code is run by typing in terminal SPADE.py config_file
Headers: FILTER, EXPTIME, TIME,RA,DEC,NAXIS1,NAXIS2
Params: NUMBER,FLUX_ISOCOR,MAG_ISOCOR,MAGERR_ISOCOR,FLUX_APER(1),FLUXERR_APER(1),
X_IMAGE, Y_IMAGE, FLAGS, FWHM_IMAGE,ELONGATION, FLUX_RADIUS, SNR_WIN, VIGNET, MAG_PSF,

"""

def Analise(File):
	SPADEFunctions.ImageProcess(File,PipelineDirectory,DayDirectory,Masters,config,logging)

TimeBetweenImages = 10
PipelineDirectory = os.getcwd()
# ~ day = date.today()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
if not os.path.isdir(DayDirectory):
	print('Não existem dados ainda deste dia')
	os.mkdir(diretoria)
	os.mkdir(diretoria + "/Darks")
	os.mkdir(diretoria + "/Flats")
	os.mkdir(diretoria + "/Images") 
	os.mkdir(diretoria + "/ProcessedImages")
	os.mkdir(diretoria + "/ReducedImages")
	os.mkdir(diretoria + "/Masters")
	os.mkdir(diretoria + "/SextractorCatalogue")
	os.mkdir(diretoria + "/AstrometryFiles")
	os.mkdir(diretoria + "/OutputFiles")


config_file = sys.argv[1]
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
config['LOGGER'] = int(config['LOGGER'])


if type(config['LOGGER']) != int or 0>config['LOGGER'] or 50<config['LOGGER']:
	print('Logger fora de range')
LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(DayDirectory,str(day))
logging.basicConfig(filename = log_file_name, level = config['LOGGER'],format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()
logging.debug("Starting create Master Files")
A = time.time()
SPADEFunctions.CreateMasters(DayDirectory,config)
logging.info("Tempos para criarMasters : {:.2f}".format(time.time()-A) )
logging.debug("Master Files Created")
Masters = SPADEFunctions.SeeMasters(DayDirectory)

# ~ logging.debug("Starting to read catalogue")
# ~ A = time.time()
# ~ catalogo,mag,BTMag = SPADEFunctions.ReadTychoCatalogue(catalog_path)
# ~ logging.info("Tempos para ler catalogo : {:.2f}".format(time.time()-A) )
# ~ logging.debug("Catalogue Read")




threadpool = concurrent.futures.ThreadPoolExecutor(max_workers = 8)
A = time.time()
while True:
	# Enquanto n tem imagens
	while len(os.listdir(path=DayDirectory + "/Images")) == 0:
		print("actives: " + str(len(threadpool._threads)))
		print("queue: " + str(threadpool._work_queue.qsize()))
		if threadpool._work_queue.qsize() + len(threadpool._threads) == 0:	
			print('There is no more data for this day')
		time.sleep(TimeBetweenImages)
		
	
	Files = os.listdir(path=DayDirectory + "/Images")
	Files.sort()
	for File in Files:				
		Path(f"{DayDirectory}/Images/{File}").rename(f"{DayDirectory}/ProcessedImages/{File}")
		threadpool.submit(Analise,File)




