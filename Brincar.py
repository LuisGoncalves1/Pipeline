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
color = ['blue','green','red','black']
size = 10

Images =	 ['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_02.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_03.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit']
# ~ Images = ['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit']

# ~ Images = ['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit']
# ~ Images = ['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit']
A = time.time()
dic = SPADEFunctions.CompareImagesJsonSkyCoord2(DayDirectory,Images,scale = 4.7)
dic_unicas = SPADEFunctions.MakeDicUnique(dic)
objetos = SPADEFunctions.SeeVariation(dic_unicas,np.radians(130/3600))
objetos = SPADEFunctions.CleanObjects(objetos)




fig,axes = plt.subplots(1,2)
ax = axes.ravel()
vmin,vmax = None,None
Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Images[0]}", unit = 'adu')
raw = Data.data
ax[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))

print(objetos)

for obj in objetos:
	for Number in objetos[obj]['Number']:
		data1 = dic[Number]
		print(data1)
		image,X,Y = int(data1['image'][0])-1, data1['Xdetect'][0], data1['Ydetect'][0]
		c1 = plt.Circle((X,Y),size, color = color[image], linewidth=1,fill=False)
		ax[1].add_patch(c1)

plt.show()

unicas = []
for i,key in enumerate(dic):
	if len(dic[key]['image']) == 1:
		unicas.append(key)	
		
		
fig,axes = plt.subplots(2,3)
ax = axes.ravel()
vmin,vmax = None,None
fig.suptitle(f'Metodo 2:{Images[0]}')
for i in range(len(Images)):
		Image = Images[i]
		Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Image}", unit = 'adu')
		raw = Data.data
		ax[i].set_title(f'Image: {i} ', fontsize=16)
		ax[i].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
ax[4].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
ax[5].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))

for i in range(len(unicas)):
	obj = dic[unicas[i]]
	image,X,Y,FWHM = int(obj['image'][0])-1, obj['Xdetect'][0], obj['Ydetect'][0], obj['FWHM'][0]	
	c1 = plt.Circle((X,Y),size, color = color[image], linewidth=0.4,fill=False)
	ax[image].add_patch(c1)
	c1 = plt.Circle((X,Y),FWHM, color = color[image], linewidth=0.4,fill=False)
	ax[4].add_patch(c1)

for i in range(len(dic)):
	obj = dic[i]
	for i in range(len(obj['image'])):
		image,X,Y,FWHM = int(obj['image'][i])-1, obj['Xdetect'][i], obj['Ydetect'][i], obj['FWHM'][i]
		c1 = plt.Circle((X,Y),2*(image+1), color = color[image], linewidth=0.4,fill=False)
		ax[5].add_patch(c1)
	
# ~ fig,axes = plt.subplots(2,2)
# ~ ax = axes.ravel()
# ~ for i in range(len(Images)):
		# ~ Image = Images[i]
		# ~ Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Image}", unit = 'adu')
		# ~ raw = Data.data
		# ~ ax[i].set_title(f'Image: {i} ', fontsize=16)
		# ~ ax[i].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,3),vmax=vmax or np.nanpercentile(raw, 99))
		
# ~ for i in range(len(dic)):
	# ~ obj = dic[i]
	# ~ for j in range(len(obj['image'])):
		# ~ image,X,Y = int(obj['image'][j])-1, obj['Xdetect'][j], obj['Ydetect'][j]	
		# ~ c1 = plt.Circle((X,Y),size, color = 'red', linewidth=0.4,fill=False)
		# ~ ax[image].add_patch(c1)
	

plt.show()		
