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
data = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(data))
config = FuncoesReducao.ReadPipelineConfig(PipelineDirectory+"/Files/Config/ConfigurationDefesa.txt")
Base = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01'


A = FuncoesReducao.ReadOutput('/home/luis/Desktop/Codigos/Pipeline/Files/Images/2022-03-20/OutputFiles/MAPA_SURV_4_085283_220417_114004+060023_scan047_01')


mask = A['RaCatalogue']!= '-'
CalculatedValues = {'RA':np.asarray(A['RaCalculated'][mask],dtype=np.float64),'DEC':np.asarray(A['DecCalculated'][mask],dtype=np.float64)}
CatalogueValues  = {'RA':np.asarray(A['RaCatalogue'][mask],dtype=np.float64),'DEC':np.asarray(A['DecCatalogue'][mask],dtype=np.float64)}
Dists = FuncoesReducao.AngleDifferenceLists(CalculatedValues,CatalogueValues)
print(Dists)

plt.hist(Dists*3600,bins=20)
plt.ylabel('Number of counts')
plt.xlabel('Distance between Calculated and Catalogue (arcsec)')
plt.show()


time.sleep(500)
#Cenas
Detections , DetectionsMags, Catalogue, CatalogueMags, found = FuncoesReducao.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)



