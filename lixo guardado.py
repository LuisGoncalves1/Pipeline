from Funcoes import FuncoesReducao,FuncoesMaster
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
import threading

"""
Headers: FILTER, EXPOSURE, TIME,RA,DEC,NAXIS1,NAXIS2


"""

def Analise(ficheiro):
	FuncoesReducao.AnalisarImagem(ficheiro,logging,config,diretoria)


TempoEntreImagens = 5
directory = os.getcwd()
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationPASO.txt")
config['LOGGER'] = int(config['LOGGER'])
# ~ data = date.today()
data = "2022-03-20"


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


if type(config['LOGGER']) != int or 0>config['LOGGER'] or 50<config['LOGGER']:
	print('Logger fora de range')

LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(diretoria,str(data))
logging.basicConfig(filename = log_file_name, level = config['LOGGER'],format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()
FuncoesMaster.CreateMasters(diretoria)
