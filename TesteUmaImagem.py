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






MagLim = 15

Imagens = ['MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fits',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fits']

n = 0
for ficheiro in Imagens:	
	Base = ficheiro.split('/')[-1].split('.')[0]
	FuncoesReducao.ImageProcess(Base,PipelineDirectory,DayDirectory,config,logging)
	

	# ~ #Desenho
	# ~ A = fits.open(DayDirectory + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat')
	# ~ Data = A[2].data
	# ~ X1 = X[found]
	# ~ Y1 = Y[found]
	# ~ X2 = X[~found]
	# ~ Y2 = Y[~found]
	# ~ mags = np.array(mag['G'])[idx]
	# ~ Imagem = CCDData.read(DayDirectory + "/ImagensReduzidas/" + ficheiro, unit = 'adu') 			# Reads Image

	# ~ fig,axes = plt.subplots(1,3)
	# ~ ax = axes.ravel()
	# ~ ax[0].imshow(Imagem,cmap='gray')
	# ~ ax[0].set_title('Objects Detected')
	# ~ ax[1].imshow(Imagem,cmap='gray')
	# ~ ax[1].set_title('Objects found in catalog')
	# ~ ax[2].imshow(Imagem,cmap='gray')
	# ~ ax[2].set_title('Objects not found in catalog')
	# ~ for contador in range(len(Y)):
		# ~ c1 = plt.Circle((X[contador],Y[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[0].add_patch(c1)
		# ~ plt.text(X[contador], Y[contador], str(round(Mags[contador],1)), color="red", fontsize=6)

	# ~ for contador in range(len(Y1)):
		# ~ c1 = plt.Circle((X1[contador],Y1[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[1].add_patch(c1)
		
	# ~ for contador in range(len(Y2)):
		# ~ c1 = plt.Circle((X2[contador],Y2[contador]),4, color = 'green', linewidth=2,fill=False)
		# ~ ax[2].add_patch(c1)
	
				
	# ~ plt.show()
