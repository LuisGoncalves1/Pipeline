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

ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01.cat'
FileName = diretoria  + '/SextractorCatalogue/' + ficheiro		


A = fits.open(FileName, mode = 'update')
A[1].data = A[1].data[A[1].data['FWHM_IMAGE'] != 0]
A.close()
		



