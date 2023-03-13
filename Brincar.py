from Functions import SPADEFunctions
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
import json

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))

Images = ['MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits','MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fits','MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fits','MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fits'] # Data set1
# ~ Images = ['PAER_SURV_4_260711_230228_042347-031345_scan014_01.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_02.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_03.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_04.fit'] # Data set3
# ~ Images = ['PAER_SURV_4_263124_230302_122520+212638_89014Df_01.fit','PAER_SURV_4_263124_230302_122520+212638_89014Df_02.fit','PAER_SURV_4_263124_230302_122520+212638_89014Df_03.fit','PAER_SURV_4_263124_230302_122520+212638_89014Df_04.fit'] # Data set2
dic = SPADEFunctions.CompareImagesJson(DayDirectory,Images,scale = 4.7)

unicas = []
for i,key in enumerate(dic):
	if len(dic[key]['image']) == 1:
		unicas.append(i)
		print(dic[key]['image'])
		
print(unicas)
print(dic[unicas])

