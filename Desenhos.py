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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
config_file = 'ConfigurationDefesa.txt'
# ~ config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
color = ['blue','green','red','black']
size = 10

DesenhoPrecisao   	= False
DesenhoFotometria 	= False
DesenhoFotometria2	= False
VerMediaTempos 		= True
FicheiroTempos = '/home/luis/Desktop/ValoresTempos/Astrometry/VersaoAntiga/Defesa(10x) - FullTime:49.458'






Images = ['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_02.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_03.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit']
# ~ Images = ['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit']
# ~ Images = ['PAER_SURV_4_260711_230228_042347-031345_scan014_01.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_02.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_03.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_04.fit']
# ~ Images = ['PAER_SURV_4_260717_230228_042347+065500_scan014_01.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_02.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_03.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_04.fit']
# ~ Images = ['PAER_SURV_4_260723_230228_042347-031345_scan015_01.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_02.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_03.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_04.fit']
# ~ Images = ['PAER_SURV_4_260729_230228_042347+065500_scan015_01.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_02.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_03.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_04.fit']
# ~ Images = ['PAER_SURV_4_260735_230228_042347-004133_scan016_01.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_02.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_03.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_04.fit']
# ~ Images = ['PAER_SURV_4_260741_230228_042347+065500_scan016_01.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_02.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_03.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_04.fit']
# ~ Images = ['PAER_SURV_4_260747_230228_042347-004133_scan017_01.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_02.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_03.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_04.fit']
# ~ Images = ['PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_02.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_03.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_04.fit']
# ~ Images = ['PAER_SURV_4_260759_230228_042347+015037_scan018_01.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_02.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_03.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_04.fit']
# ~ Images = ['PAER_SURV_4_260765_230228_042347-081807_scan019_01.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_02.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_03.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_04.fit']
# ~ Images = ['PAER_SURV_4_260771_230228_042347+015037_scan019_01.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_02.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_03.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_04.fit']
# ~ Images = ['PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_02.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_03.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_04.fit']
# ~ Images = ['PAER_SURV_4_260783_230228_080700+045703_scan101_01.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_02.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_03.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_04.fit']


if DesenhoPrecisao:
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				RaCalc.append(float(Object['RaCalculated']))
				DecCalc.append(float(Object['DecCalculated']))
				RaCat.append(float(Object['RaCatalogue']))
				DecCat.append(float(Object['DecCatalogue']))
				



	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = SPADEFunctions.AngleDifferenceLists(L1,L2)
	plt.hist(Distance*3600,bins=2000)
	plt.xlim([0,20])
	plt.ylabel('Number of Point',fontsize = 20)
	plt.xlabel('Distance between Points Calculated and Catalogue (arcsec)',fontsize = 20)
	plt.show()
	
	fig,axes = plt.subplots(1,1)
	DiferencasRa  = np.array(RaCalc)-np.array(RaCat)
	print(np.std(3600*DiferencasRa))
	print(3600*np.median(DiferencasRa))
	print('Ra')
	DiferencasDec = np.array(DecCalc)-np.array(DecCat)
	print(np.std(3600*DiferencasDec))
	print(3600*np.median(DiferencasDec))
	print('Dec')

	for i in range(len(RaCalc)):
		c1 = plt.Circle((3600*(RaCalc[i]-RaCat[i]),3600*(DecCalc[i]-DecCat[i])),0.1, color = 'green', linewidth=2,fill=False)
		axes.add_patch(c1)
		
	axes.set_xlim([-20,20])
	axes.set_ylim([-20,20])
	axes.set_xlabel('Difference in RA (arcsec)')
	axes.set_ylabel('Difference in Dec (arcsec)')
	axes.set_aspect('equal', adjustable='box')
	axes.axhline(0, color='black',lw=1)
	axes.axvline(0, color='black',lw=1)
	plt.show()

if DesenhoFotometria:
	MagCalc,MagCat = [],[]
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				MagCalc.append(float(Object['MagInst ']))
				MagCat.append(float(Object['MagReal']))


	plt.plot(MagCalc,MagCat, 'o')
	plt.xlabel('Mag Instrumental')
	plt.ylabel('Mag Catalogo')
	plt.show()

if DesenhoFotometria2:
	MagCalc1,MagCat1 = [],[]
	MagCalc2,MagCat2 = [],[]
	Dec = []
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				RaCalc  = float(Object['RaCalculated'])
				DecCalc = float(Object['DecCalculated'])
				RaCat   = float(Object['RaCatalogue'])
				DecCat  = float(Object['DecCatalogue'])
				Dec.append(DecCalc)
				dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat])
				
				#5*4.7
				if dist > 0.1/3600:
					MagCalc1.append(float(Object['MagInst ']))
					MagCat1.append(float(Object['MagReal']))
				else:
					MagCalc2.append(float(Object['MagInst ']))
					MagCat2.append(float(Object['MagReal']))

	print(len(MagCalc1))
	print(len(MagCalc2))
	print(min(Dec))
	print(max(Dec))
	plt.plot(MagCalc1,MagCat1, 'or')
	plt.plot(MagCalc2,MagCat2, 'ob')
	plt.xlabel('Mag Instrumental')
	plt.ylabel('Mag Catalogo')
	plt.show()

if VerMediaTempos:
	tempos = {}
	with open(FicheiroTempos,'r') as ficheiro:
		Tipos = ficheiro.readline()
		Tipos = Tipos.strip().split(',')
		for Tipo in Tipos:
			tempos[Tipo] = []
		while True:
			linha = ficheiro.readline()
			if linha == '':
				break
			linha = linha.strip().split(',')
			for i,tempo in enumerate(linha):
				tempos[Tipos[i]].append(float(tempo))
	
	mediaTempos = {}
	for i in tempos:
		mediaTempos[i] = np.mean(tempos[i])
		
	print(mediaTempos)




















# ~ SexFile = '/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/SextractorCatalogue/PAER_SURV_4_260789_230228_080700-051142_scan102_01.cat'
# ~ SexData = fits.open(SexFile)
# ~ Data = SexData[1].data
# ~ X,Y,Number = Data['X_IMAGE'],Data['Y_IMAGE'],Data['ISOAREA_IMAGE']
# ~ fig,axes = plt.subplots(1,3)
# ~ ax = axes.ravel()
# ~ vmin,vmax = None,None
# ~ Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Images[0]}", unit = 'adu')
# ~ raw = Data.data
# ~ ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
# ~ ax[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
# ~ ax[2].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))

# ~ for i in range(len(X)):
	# ~ c1 = plt.Circle((X[i],Y[i]),5, color = 'green', linewidth=2,fill=False)
	# ~ ax[1].add_patch(c1)
	# ~ ax[1].text(X[i],Y[i],Number[i])



# ~ plt.show()

