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





#Cenas
Detections , DetectionsMags, Catalogue, CatalogueMags, found = FuncoesReducao.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)



