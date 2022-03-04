from Funcoes import FuncoesReducao
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
# ~ import astroscrappy
import threading

"""
Headers: FILTER, EXPOSURE, TIME,


"""

def Analise(ficheiro):
	FuncoesReducao.AnalisarImagem(ficheiro,logging,config,diretoria)


TempoEntreImagens = 5
directory = os.getcwd()
config = FuncoesReducao.LerFicheiro(directory+"/Ficheiros/Config/ConfigurationPASO.txt")
config['LOGGER'] = int(config['LOGGER'])
data = date.today()

diretoria = "%s/Ficheiros/Imagens/%s" %(directory,str(data))
if not os.path.isdir(diretoria):
	print('Não existem dados ainda deste dia')
	os.mkdir(diretoria)
	os.mkdir(diretoria + "/Darks")
	os.mkdir(diretoria + "/Flats")
	os.mkdir(diretoria + "/Imagens")
	os.mkdir(diretoria + "/ImagensTratadas")
	os.mkdir(diretoria + "/ImagensReduzidas")


if type(config['LOGGER']) != int or 0>config['LOGGER'] or 50<config['LOGGER']:
	print('Logger fora de range')

LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(diretoria,str(data))
logging.basicConfig(filename = log_file_name, level = config['LOGGER'],format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()

if config["DARK"]:
	if len(os.listdir(path=diretoria + "/Darks")) == 0:
		logging.warning("Não existem Darks, removendo este passo da reducao")
		config["DARK"] = False
	else:
		logging.debug("Criando os MasterDarks")
		for filtro in os.listdir(path=diretoria + "/Darks"):
			if not os.path.isfile(diretoria + "/MasterDark" + filtro + ".fits"):
				logging.debug("Criando o MasterDark do filtro " + filtro)
				MasterDark = FuncoesReducao.CriarMaster(diretoria + "/Darks/" + filtro + '/', diretoria, filtro,'dark')
				logging.debug("MasterDark Criado do filtro " + filtro)
			else:
				MasterDark = CCDData.read(diretoria + "/MasterDark" + filtro + ".fits") 
		logging.debug("MasterDarks Criados")
																
if config["FLAT"]:
	if len(os.listdir(path=diretoria + "/Flats")) == 0:
		logging.warning("Não existem Flats, removendo este passo da reducao")
		config["FLAT"] = False
	else:		
		logging.debug("Criando os MasterFlats")
		for filtro in os.listdir(path=diretoria + "/Flats"):		
			if not os.path.isfile(diretoria + "/MasterFlat" + filtro + ".fits"):
				logging.debug("Criando o MasterFlat do filtro " + filtro)
				MasterFlat = FuncoesReducao.CriarMaster(diretoria + "/Flats/" + filtro + '/', diretoria, filtro,'flat')
				logging.debug("MasterFlat Criado do filtro " + filtro)
			else:
				MasterFlat = CCDData.read(diretoria + "/MasterFlat" + filtro + ".fits") 
		logging.debug("MasterFlats Criados")

while True:
	# Enquanto n tem imagens
	while len(os.listdir(path=diretoria + "/Imagens")) == 0:
		print('Não existem mais dados deste dia')
		time.sleep(TempoEntreImagens)
	
	for ficheiro in os.listdir(path=diretoria + "/Imagens"):
		Path(diretoria + '/Imagens/' + ficheiro).rename(diretoria + '/ImagensTratadas/' + ficheiro)
		t = threading.Thread(target = Analise,args = (ficheiro,))
		t.start()
		
		
		
