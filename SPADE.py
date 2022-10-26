from Funcoes import FuncoesReducao, FuncoesMaster
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

"""
Headers: FILTER, EXPTIME, TIME,RA,DEC,NAXIS1,NAXIS2


"""

def Analise(ficheiro,Masters):
	FuncoesReducao.AnalisarImagem(ficheiro,Masters,logging,config,diretoria,directory,catalogo,mag)



TempoEntreImagens = 1
directory = os.getcwd()
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationBruno.txt")
config['LOGGER'] = int(config['LOGGER'])
# ~ data = date.today()
data = "2022-03-20"

catalog_path = directory + '/Files/Catalogos/tycho2.dat'


diretoria = "%s/Files/Imagens/%s" %(directory,str(data))
if not os.path.isdir(diretoria):
	print('Não existem dados ainda deste dia')
	os.mkdir(diretoria)
	os.mkdir(diretoria + "/Darks")
	os.mkdir(diretoria + "/Flats")
	os.mkdir(diretoria + "/Imagens") 
	os.mkdir(diretoria + "/ImagensTratadas")
	os.mkdir(diretoria + "/ImagensReduzidas")
	os.mkdir(diretoria + "/Masters")
	os.mkdir(diretoria + "/SextractorCatalogue")
	os.mkdir(diretoria + "/AstrometryFiles")

if type(config['LOGGER']) != int or 0>config['LOGGER'] or 50<config['LOGGER']:
	print('Logger fora de range')
LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(diretoria,str(data))
logging.basicConfig(filename = log_file_name, level = config['LOGGER'],format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()
logging.debug("Starting create Master Files")
A = time.time()
FuncoesMaster.CreateMasters(diretoria,config)
logging.info("Tempos para criarMasters : {:.2f}".format(time.time()-A) )
logging.debug("Master Files Created")
Masters = FuncoesMaster.SeeMasters(diretoria)

logging.debug("Starting to read catalogue")
A = time.time()
catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)
logging.info("Tempos para ler catalogo : {:.2f}".format(time.time()-A) )
logging.debug("Catalogue Read")







threadpool = concurrent.futures.ThreadPoolExecutor(max_workers = 2)
A = time.time()
while True:
	# Enquanto n tem imagens
	while len(os.listdir(path=diretoria + "/Imagens")) == 0:
		print("ativos: " + str(len(threadpool._threads)))
		print("queue: " + str(threadpool._work_queue.qsize()))
			# ~ logging.info("Tempos para tratar todas imagens : {:.2f}".format(time.time()-A) )		
		if threadpool._work_queue.qsize() + len(threadpool._threads) == 0:	
			print('Não existem mais dados deste dia')
		time.sleep(TempoEntreImagens)
		
		
	for ficheiro in os.listdir(path=diretoria + "/Imagens"):				
		Path(diretoria + '/Imagens/' + ficheiro).rename(diretoria + '/ImagensTratadas/' + ficheiro)
		threadpool.submit(Analise,ficheiro,Masters)


