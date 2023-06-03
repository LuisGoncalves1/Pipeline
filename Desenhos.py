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
fontsize = 25
ticksize = 15

DesenhoPrecisaoTodos   		= False		# Faz Desenho da diferença Ra e Dec para todas as imagens, junto em 1 gráfico
DesenhoPrecisaoTodosEspecial= False
DesenhoPrecisaoUma   		= False		# Faz Desenho da diferença Ra e Dec de um set num 1 gráfico
DesenhoPrecisaoTodosSep		= True		# Faz Desenho da diferença Ra e Dec para todas as imagens, separado em 15 gráficos
DesenhoPrecisaoTodosSep2	= False
DesenhoDiferentesGrupos		= False
DesenhoDiferentesGruposUnica= False
DesenhoSetsEspecificos		= False

DesenhoFotometriaUnica 		= False
DesenhoFotometriaUnica2		= False
DesenhoFotometriaUnica3		= False
FotometriaUnicaDeitada		= False
FotometriaUnicaDeitada2		= False
FotometriaUnicaDeitada3		= False
VerMediaTempos 				= False


DesenhoImage				= False
DesenhoImages				= False
DesenhoImages2				= False
DesenhoSet					= False
DesenhoEncontrados			= False
DesenhoMaster				= False
DesenhoReduzida				= False
DesenhoVariosSigmas			= False
FicheiroTempos = '/home/luis/Desktop/ValoresTempos/Detecao/5Sigma'


#apagar:
DesenhoFotometria 			= False
DesenhoFotometria2			= False

# ~ 3, 4, 7, 9, 11, 12, 13, 15 Sets especiais

# ~ Images = ['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_02.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_03.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit']
Images = ['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit']
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

Sets = [['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_02.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_03.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit'],
['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit'],
['PAER_SURV_4_260711_230228_042347-031345_scan014_01.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_02.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_03.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_04.fit'],
['PAER_SURV_4_260717_230228_042347+065500_scan014_01.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_02.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_03.fit','PAER_SURV_4_260717_230228_042347+065500_scan014_04.fit'],
['PAER_SURV_4_260723_230228_042347-031345_scan015_01.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_02.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_03.fit','PAER_SURV_4_260723_230228_042347-031345_scan015_04.fit'],
['PAER_SURV_4_260729_230228_042347+065500_scan015_01.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_02.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_03.fit','PAER_SURV_4_260729_230228_042347+065500_scan015_04.fit'],
['PAER_SURV_4_260735_230228_042347-004133_scan016_01.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_02.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_03.fit','PAER_SURV_4_260735_230228_042347-004133_scan016_04.fit'],
['PAER_SURV_4_260741_230228_042347+065500_scan016_01.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_02.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_03.fit','PAER_SURV_4_260741_230228_042347+065500_scan016_04.fit'],
['PAER_SURV_4_260747_230228_042347-004133_scan017_01.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_02.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_03.fit','PAER_SURV_4_260747_230228_042347-004133_scan017_04.fit'],
['PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_02.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_03.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_04.fit'],
['PAER_SURV_4_260759_230228_042347+015037_scan018_01.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_02.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_03.fit','PAER_SURV_4_260759_230228_042347+015037_scan018_04.fit'],
['PAER_SURV_4_260765_230228_042347-081807_scan019_01.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_02.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_03.fit','PAER_SURV_4_260765_230228_042347-081807_scan019_04.fit'],
['PAER_SURV_4_260771_230228_042347+015037_scan019_01.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_02.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_03.fit','PAER_SURV_4_260771_230228_042347+015037_scan019_04.fit'],
['PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_02.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_03.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_04.fit'],
['PAER_SURV_4_260783_230228_080700+045703_scan101_01.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_02.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_03.fit','PAER_SURV_4_260783_230228_080700+045703_scan101_04.fit']]



if DesenhoPrecisaoTodos:
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	method = 'normal'
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				RaCalc.append(float(Object['RaCalculated']))
				DecCalc.append(float(Object['DecCalculated'])-180)
				RaCat.append(float(Object['RaCatalogue']))
				DecCat.append(float(Object['DecCatalogue'])-180)
				


	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = 3600*SPADEFunctions.AngleDifferenceLists(L1,L2,method)
	plt.hist(Distance,bins=2000)
	mask = Distance<2
	print(sum(mask))
	print(len(mask))
	print(sum(mask)/len(mask))
	plt.xlim([0,12.5])
	plt.ylabel('Number of Points',fontsize = 20)
	plt.xlabel('Distance between Points Calculated and Catalogue (arcsec)',fontsize = 20)
	plt.xticks(fontsize=ticksize)
	plt.yticks(fontsize=ticksize)
	plt.axvline(2, color='orange',lw=1)
	
	percent = 100* sum(mask)/len(mask)
	col_labels = ['%']
	row_labels = ['>2%','<2%']
	table_vals = [[round(100-percent,1)],[round(percent,1)]]
	row_colors = ['red', 'blue', 'yellow']
	#plotting
	my_table = plt.table(cellText=table_vals,
						 colWidths=[0.1],
						 rowLabels=row_labels,
						 colLabels=col_labels,
						 loc='upper right',
						 cellLoc='center')
	my_table.auto_set_font_size(False)
	my_table.set_fontsize(ticksize)
	fig,axes = plt.subplots(1,1)
	DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
	DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))

	

	for i in range(len(RaCalc)):
		c1 = plt.Circle((DiferencasRa[i],DiferencasDec[i]),0.00005, color = 'green', linewidth=2,fill=False)
		axes.add_patch(c1)
	
	medDec,medRa = np.mean(DiferencasDec),np.mean(DiferencasRa)
	stdDec,stdRa = np.std(DiferencasDec),np.std(DiferencasRa)
	print('Valores Dec: ',medDec, ',', stdDec)
	print('Valores Ra: ',medRa, ',', stdRa)
	axes.set_xlim([-10,10])
	axes.set_ylim([-10,10])
	axes.set_xlabel('Difference in RA (arcsec)',fontsize=fontsize)
	axes.set_ylabel('Difference in Dec (arcsec)',fontsize=fontsize)
	axes.set_aspect('equal', adjustable='box')
	axes.axhline(medDec, color='black',lw=1)
	axes.axvline(medRa, color='black',lw=1)
	
	axes.axhline(medDec+stdDec, color='orange',lw=1)
	axes.axhline(medDec-stdDec, color='orange',lw=1)
	axes.axvline(medRa+stdRa, color='orange',lw=1)
	axes.axvline(medRa-stdRa, color='orange',lw=1)
	plt.xticks(fontsize=ticksize)
	plt.yticks(fontsize=ticksize)

if DesenhoPrecisaoTodosEspecial:
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	method = 'normal'
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				RaCalc.append(float(Object['RaCalculated']))
				DecCalc.append(float(Object['DecCalculated'])-180)
				RaCat.append(float(Object['RaCatalogue']))
				DecCat.append(float(Object['DecCatalogue'])-180)
				


	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = SPADEFunctions.AngleDifferenceLists(L1,L2,method)
	x  = 3600*(np.array(RaCalc)-np.array(RaCat))	
	y = 3600*(np.array(DecCalc)-np.array(DecCat))

	fig = plt.figure(constrained_layout=True,figsize=(6, 6))
	grid = plt.GridSpec(4, 4, hspace=0, wspace=0)

	main_ax = fig.add_subplot(grid[1:, :-1])
	main_ax.set_xlim([-10, 10])
	main_ax.set_ylim([-10, 10])
	main_ax.set_aspect('equal', adjustable='box')

	y_hist = fig.add_subplot(grid[1:, -1], xticklabels=[], sharey=main_ax)
	x_hist = fig.add_subplot(grid[0, :-1], yticklabels=[], sharex=main_ax)

	# scatter points on the main axes
	main_ax.plot(x, y, 'og', markersize=0.3, alpha=0.2)

	# histogram on the attached axes
	x_hist.hist(x, 2000, histtype='stepfilled',
				orientation='vertical', color='green')
	# ~ x_hist.invert_yaxis()

	y_hist.hist(y, 2000, histtype='stepfilled',
				orientation='horizontal', color='green')
	# ~ y_hist.invert_xaxis()
	
	# Remove ticks from the top histogram
	x_hist.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

	# Remove ticks from the right histogram
	y_hist.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

	main_ax.set_xlabel('Difference in Ra (arcsec)',fontsize=12)
	main_ax.set_ylabel('Difference in Dec (arcsec)',fontsize=12)
	
	medDec,medRa = np.mean(y),np.mean(x)
	stdDec,stdRa = np.std(y),np.std(x)
	print('Valores Dec: ',medDec, ',', stdDec)
	print('Valores Ra: ',medRa, ',', stdRa)
	print(max(x))
	print(max(y))
	main_ax.axhline(medDec, color='black',lw=1)
	main_ax.axvline(medRa, color='black',lw=1)
	
	main_ax.axhline(medDec+stdDec, color='orange',lw=1)
	main_ax.axhline(medDec-stdDec, color='orange',lw=1)
	main_ax.axvline(medRa+stdRa, color='orange',lw=1)
	main_ax.axvline(medRa-stdRa, color='orange',lw=1)
	plt.xticks(fontsize=ticksize)
	plt.yticks(fontsize=ticksize)
	
	plt.show()

if DesenhoPrecisaoUma:
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	method = 'normal'
	for Output in Images:
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				RaCalc.append(float(Object['RaCalculated']))
				DecCalc.append(float(Object['DecCalculated']))
				RaCat.append(float(Object['RaCatalogue']))
				DecCat.append(float(Object['DecCatalogue']))
				



	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = SPADEFunctions.AngleDifferenceLists(L1,L2,method)
	plt.hist(Distance*3600,bins=2000)
	plt.xlim([0,12.5])
	plt.ylabel('Number of Points',fontsize = 20)
	plt.xlabel('Distance between Points Calculated and Catalogue (arcsec)',fontsize = 20)
	
	fig,axes = plt.subplots(1,1)
	DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
	DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))
	
	

	for i in range(len(RaCalc)):
		c1 = plt.Circle((DiferencasRa[i],DiferencasDec[i]),0.005, color = 'green', linewidth=2,fill=False)
		axes.add_patch(c1)
	
	medDec,medRa = np.mean(DiferencasDec),np.mean(DiferencasRa)
	stdDec,stdRa = np.std(DiferencasDec),np.std(DiferencasRa)
	print('Valores Dec: ',medDec, ',', stdDec)
	print('Valores Ra: ',medRa, ',', stdRa)
	# ~ axes.set_xlim([-20,20])
	# ~ axes.set_ylim([-20,20])
	axes.set_xlabel('Difference in RA (arcsec)',fontsize=fontsize)
	axes.set_ylabel('Difference in Dec (arcsec)',fontsize=fontsize)
	axes.set_aspect('equal', adjustable='box')
	axes.axhline(medDec, color='black',lw=1)
	axes.axvline(medRa, color='black',lw=1)
	
	axes.axhline(medDec+stdDec, color='orange',lw=1)
	axes.axhline(medDec-stdDec, color='orange',lw=1)
	axes.axvline(medRa+stdRa, color='orange',lw=1)
	axes.axvline(medRa-stdRa, color='orange',lw=1)

if DesenhoPrecisaoTodosSep:
	fig1,axes1 = plt.subplots(4,4)
	ax1 = axes1.ravel()
	# ~ plt.xticks(fontsize=ticksize)
	# ~ plt.yticks(fontsize=ticksize)
	fig2,axes2 = plt.subplots(4,4)
	ax2 = axes2.ravel()
	# ~ plt.xticks(fontsize=ticksize)
	# ~ plt.yticks(fontsize=ticksize)
	for count,Images in enumerate(Sets):
		RaCalc,DecCalc = [],[]
		RaCat,DecCat   = [],[]
		method = 'normal'
		for Output in Images:
			with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
				ImageData = json.load(stream)
			for Object in ImageData:
				if Object['Match']:
					RaCalc.append(float(Object['RaCalculated']))
					DecCalc.append(float(Object['DecCalculated']))
					RaCat.append(float(Object['RaCatalogue']))
					DecCat.append(float(Object['DecCatalogue']))
					

		# Primeiro Desenho -> Histogramas
		L1 = {'RA':RaCalc,'DEC':DecCalc}
		L2 = {'RA':RaCat, 'DEC':DecCat}
		Distance = 3600*SPADEFunctions.AngleDifferenceLists(L1,L2,method)
		ax1[count].hist(Distance,bins=2000)
		if count != 13:
			ax1[count].set_xlim([0,15])
		fig1.supxlabel('Distance between Points Calculated and Catalogue (arcsec)',fontsize = fontsize/2)
		fig1.supylabel('Number of Points',fontsize = fontsize/2)
		
		NTotal 	= len(RaCalc)
		Ncentro = sum(Distance<3)
		print('Set: ',count+1)
		print(NTotal)
		print(Ncentro)
		print(Ncentro/NTotal)
		print('----')
		
		
		
		# Segundo Desenho -> Precisão
		DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
		DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))		
		for i in range(len(RaCalc)):
			c1 = plt.Circle((DiferencasRa[i],DiferencasDec[i]),0.005, color = 'green', linewidth=2,fill=False)
			ax2[count].add_patch(c1)
		
		medDec,medRa = np.mean(DiferencasDec),np.mean(DiferencasRa)
		stdDec,stdRa = np.std(DiferencasDec),np.std(DiferencasRa)
		print(stdDec)
		print(stdRa)
		
		# ~ axes.set_xlim([-20,20])
		# ~ axes.set_ylim([-20,20])
		fig2.supxlabel('Difference in RA (arcsec)',fontsize=fontsize/2)
		fig2.supylabel('Difference in Dec (arcsec)',fontsize=fontsize/2)

		ax2[count].set_aspect('equal', adjustable='box')
		ax2[count].axhline(medDec, color='black',lw=1)
		ax2[count].axvline(medRa, color='black',lw=1)
		
		ax2[count].axhline(medDec+stdDec, color='orange',lw=1)
		ax2[count].axhline(medDec-stdDec, color='orange',lw=1)
		ax2[count].axvline(medRa+stdRa, color='orange',lw=1)
		ax2[count].axvline(medRa-stdRa, color='orange',lw=1)
		
		if count != 13:
			ticks = [-10,-5,0,5,10]
			ax2[count].set_xlim([-10,10])
			ax2[count].set_ylim([-10,10])
			ax2[count].set_xticks(ticks)
		
		ax1[count].text(0.8, 0.7, str(count+1), fontsize=12, ha='center', va='center', transform=ax1[count].transAxes,
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square'))

		ax2[count].text(0.9, 0.9, str(count+1), fontsize=12, ha='center', va='center', transform=ax2[count].transAxes,
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square'))
	
	ax1[15].remove()
	ax2[15].remove()		

if DesenhoPrecisaoTodosSep2:
	Percentil = 0.85
	
	
	
	fig1,axes1 = plt.subplots(4,4)
	ax1 = axes1.ravel()
	fig2,axes2 = plt.subplots(4,4)
	ax2 = axes2.ravel()
	for count,Images in enumerate(Sets):
		RaCalc,DecCalc = [],[]
		RaCat,DecCat   = [],[]
		method = 'normal'
		for Output in Images:
			with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
				ImageData = json.load(stream)
			for Object in ImageData:
				if Object['Match']:
					RaCalc.append(float(Object['RaCalculated']))
					DecCalc.append(float(Object['DecCalculated']))
					RaCat.append(float(Object['RaCatalogue']))
					DecCat.append(float(Object['DecCatalogue']))
					

		# Primeiro Desenho -> Histogramas
		L1 = {'RA':RaCalc,'DEC':DecCalc}
		L2 = {'RA':RaCat, 'DEC':DecCat}
		Distance = SPADEFunctions.AngleDifferenceLists(L1,L2,method)
		ax1[count].hist(Distance*3600,bins=2000)
		if count != 13:
			ax1[count].set_xlim([0,15])
		fig1.supxlabel('Distance between Points Calculated and Catalogue (arcsec)',fontsize = fontsize)
		fig1.supylabel('Number of Points',fontsize = fontsize)

		Distance2 = Distance.copy()
		Distance2.sort()
		Distance2 = Distance2[0:int(Percentil*len(Distance2))]
		maxDistance = Distance2[-1]
		mask = Distance < maxDistance
		print('Set: ',count+1)
		print(3600*maxDistance)
		print('---')
		
		# Segundo Desenho -> Precisão
		DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
		DiferencasRa = DiferencasRa[mask]
		DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))
		DiferencasDec = DiferencasDec[mask]
		for i in range(len(DiferencasRa)):
			c1 = plt.Circle((DiferencasRa[i],DiferencasDec[i]),0.005, color = 'green', linewidth=2,fill=False)
			ax2[count].add_patch(c1)
		
		medDec,medRa = np.mean(DiferencasDec),np.mean(DiferencasRa)
		stdDec,stdRa = np.std(DiferencasDec),np.std(DiferencasRa)
		
		# ~ axes.set_xlim([-20,20])
		# ~ axes.set_ylim([-20,20])
		fig2.supxlabel('Difference in RA (arcsec)',fontsize=fontsize)
		fig2.supylabel('Difference in Dec (arcsec)',fontsize=fontsize)

		ax2[count].set_aspect('equal', adjustable='box')
		ax2[count].axhline(medDec, color='black',lw=1)
		ax2[count].axvline(medRa, color='black',lw=1)
		
		ax2[count].axhline(medDec+stdDec, color='orange',lw=1)
		ax2[count].axhline(medDec-stdDec, color='orange',lw=1)
		ax2[count].axvline(medRa+stdRa, color='orange',lw=1)
		ax2[count].axvline(medRa-stdRa, color='orange',lw=1)
		
		if count != 13:
			ticks = [-10,-5,0,5,10]
			ax2[count].set_xlim([-10,10])
			ax2[count].set_ylim([-10,10])
			ax2[count].set_xticks(ticks)
	
	ax1[15].remove()
	ax2[15].remove()		

if DesenhoSetsEspecificos:
	Percentil = 0.85
	Sets2 = [['PAER_SURV_4_260795_230228_080700+022452_scan102_01.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_02.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_03.fit','PAER_SURV_4_260795_230228_080700+022452_scan102_04.fit'],
	['PAER_SURV_4_260711_230228_042347-031345_scan014_01.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_02.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_03.fit','PAER_SURV_4_260711_230228_042347-031345_scan014_04.fit'],
	['PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_02.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_03.fit','PAER_SURV_4_260777_230228_080700-074353_scan101_04.fit']]	
	
	fig2,axes2 = plt.subplots(1,3)
	ax2 = axes2.ravel()
	for count,Images in enumerate(Sets2):
		RaCalc,DecCalc = [],[]
		RaCat,DecCat   = [],[]
		method = 'normal'
		for Output in Images:
			with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
				ImageData = json.load(stream)
			for Object in ImageData:
				if Object['Match']:
					RaCalc.append(float(Object['RaCalculated']))
					DecCalc.append(float(Object['DecCalculated']))
					RaCat.append(float(Object['RaCatalogue']))
					DecCat.append(float(Object['DecCatalogue']))
					
		L1 = {'RA':RaCalc,'DEC':DecCalc}
		L2 = {'RA':RaCat, 'DEC':DecCat}
		Distance = SPADEFunctions.AngleDifferenceLists(L1,L2,method)
		Distance2 = Distance.copy()
		Distance2.sort()
		Distance2 = Distance2[0:int(Percentil*len(Distance2))]
		maxDistance = Distance2[-1]
		mask = Distance < maxDistance
		print('Set: ',count+1)
		print(3600*maxDistance)
		print('---')
		
		# Segundo Desenho -> Precisão
		DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
		DiferencasRa = DiferencasRa[mask]
		DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))
		DiferencasDec = DiferencasDec[mask]
		for i in range(len(DiferencasRa)):
			c1 = plt.Circle((DiferencasRa[i],DiferencasDec[i]),0.005, color = 'green', linewidth=2,fill=False)
			ax2[count].add_patch(c1)
		
		medDec,medRa = np.mean(DiferencasDec),np.mean(DiferencasRa)
		stdDec,stdRa = np.std(DiferencasDec),np.std(DiferencasRa)
		
		# ~ axes.set_xlim([-20,20])
		# ~ axes.set_ylim([-20,20])
		suplabel_y = +0.001
		fig2.supxlabel('Difference in RA (arcsec)',fontsize=fontsize)
		fig2.supylabel('Difference in Dec (arcsec)',fontsize=fontsize)

		ax2[count].set_aspect('equal', adjustable='box')
		ax2[count].axhline(medDec, color='black',lw=1)
		ax2[count].axvline(medRa, color='black',lw=1)
		
		ax2[count].axhline(medDec+stdDec, color='orange',lw=1)
		ax2[count].axhline(medDec-stdDec, color='orange',lw=1)
		ax2[count].axvline(medRa+stdRa, color='orange',lw=1)
		ax2[count].axvline(medRa-stdRa, color='orange',lw=1)
		
		if count != 13:
			ticks = [-10,-5,0,5,10]
			ax2[count].set_xlim([-10,10])
			ax2[count].set_ylim([-10,10])
			ax2[count].set_xticks(ticks)
		ax2[count].tick_params(axis='both', labelsize=ticksize)
		
	ax2[0].text(0.9, 0.9, str(2), fontsize=12, ha='center', va='center', transform=ax2[0].transAxes,bbox=dict(facecolor='white', edgecolor='black', boxstyle='square'))
	ax2[1].text(0.9, 0.9, str(3), fontsize=12, ha='center', va='center', transform=ax2[1].transAxes,bbox=dict(facecolor='white', edgecolor='black', boxstyle='square'))
	ax2[2].text(0.9, 0.9, str(14), fontsize=12, ha='center', va='center', transform=ax2[2].transAxes,bbox=dict(facecolor='white', edgecolor='black', boxstyle='square'))	
	
if DesenhoDiferentesGruposUnica:
	#711,717,753,777,783
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	X,Y,Image	   = [],[],[]
	method = 'normal'
	fig,axes = plt.subplots(1,1)
	vmin,vmax = None,None
	Output = Images[1]
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
		ImageData = json.load(stream)
		Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Output}")
		raw = Data.data
		axes.imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
		axes.set_axis_off()
	for Object in ImageData:
		if Object['Match']:
			RaCalc.append(float(Object['RaCalculated']))
			DecCalc.append(float(Object['DecCalculated']))
			RaCat.append(float(Object['RaCatalogue']))
			DecCat.append(float(Object['DecCatalogue']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Image.append(int(Object['image']))



	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = 3600*SPADEFunctions.AngleDifferenceLists(L1,L2,method)
	DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
	DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))
	
	# ~ mask = Distance>5
	# ~ mask = DiferencasRa<2
	mask = DiferencasDec>3
	color = 'blue'
	for i in range(len(RaCalc)):
		if mask[i]:
			c1 = plt.Circle((X[i],Y[i]), 3, color = color, linewidth=2,fill=False)
			axes.add_patch(c1)
	
	mask = DiferencasDec<-3
	color = 'red'
	for i in range(len(RaCalc)):
		if mask[i]:
			c1 = plt.Circle((X[i],Y[i]), 3, color = color, linewidth=2,fill=False)
			axes.add_patch(c1)

if DesenhoDiferentesGrupos:
	#711,717,753,777,783
	RaCalc,DecCalc = [],[]
	RaCat,DecCat   = [],[]
	X,Y,Image	   = [],[],[]
	method = 'normal'
	fig,axes = plt.subplots(2,2)
	ax = axes.ravel()
	vmin,vmax = None,None
	for i,Output in enumerate(Images):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output.split('.')[0]+'.json'}", 'r') as stream:
			ImageData = json.load(stream)
			Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Output}")
			raw = Data.data
			ax[i].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
			ax[i].set_axis_off()
		for Object in ImageData:
			if Object['Match']:
				RaCalc.append(float(Object['RaCalculated']))
				DecCalc.append(float(Object['DecCalculated']))
				RaCat.append(float(Object['RaCatalogue']))
				DecCat.append(float(Object['DecCatalogue']))
				X.append(float(Object['Xdetect']))
				Y.append(float(Object['Ydetect']))
				Image.append(int(Object['image']))



	L1 = {'RA':RaCalc,'DEC':DecCalc}
	L2 = {'RA':RaCat, 'DEC':DecCat}
	Distance = 3600*SPADEFunctions.AngleDifferenceLists(L1,L2,method)
	DiferencasRa  = 3600*(np.array(RaCalc)-np.array(RaCat))	
	DiferencasDec = 3600*(np.array(DecCalc)-np.array(DecCat))
	
	# ~ mask = Distance>5
	# ~ mask = DiferencasRa<2
	mask = DiferencasDec>3
	color = 'blue'
	for i in range(len(RaCalc)):
		if mask[i]:
			c1 = plt.Circle((X[i],Y[i]), 1, color = color, linewidth=2,fill=False)
			ax[Image[i]-1].add_patch(c1)
	
	mask = DiferencasDec<-3
	color = 'red'
	for i in range(len(RaCalc)):
		if mask[i]:
			c1 = plt.Circle((X[i],Y[i]), 1, color = color, linewidth=2,fill=False)
			ax[Image[i]-1].add_patch(c1)
						
if DesenhoFotometria:
	MagCalc,MagCat = [],[]
	for Output in os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles'):
		with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Output}", 'r') as stream:
			ImageData = json.load(stream)
		for Object in ImageData:
			if Object['Match']:
				MagCalc.append(float(Object['MagInst']))
				MagCat.append(float(Object['MagReal']))


	plt.plot(MagCalc,MagCat, 'o')
	plt.xlabel('Mag Instrumental')
	plt.ylabel('Mag Catalogo')

if DesenhoFotometria2:
	Threshold = 0.25
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
				dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
				
				#5*4.7
				if dist > Threshold/3600:
					MagCalc1.append(float(Object['MagInst']))
					MagCat1.append(float(Object['MagReal']))
				else:
					MagCalc2.append(float(Object['MagInst']))
					MagCat2.append(float(Object['MagReal']))

	print(len(MagCalc2)/len(os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')))
	print(len(MagCat1)/len(os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')))
	plt.plot(MagCalc1,MagCat1, 'or')
	plt.plot(MagCalc2,MagCat2, 'ob')
	plt.xlabel('Mag Instrumental')
	plt.ylabel('Mag Catalogo')
	plt.title(f"Threshold: {Threshold}")

if DesenhoFotometriaUnica:
	MagCalc,MagCat = [],[]
	Imagens = os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')
	Imagem = Imagens[7]
	print(Imagem)
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))


	plt.plot(MagCalc,MagCat, 'o',markersize=4)
	plt.xlabel('Mag Instrumental',fontsize=fontsize)
	plt.ylabel('Mag Catalogo',fontsize=fontsize)

if DesenhoFotometriaUnica2:
	Threshold = 0.2
	MagCalc1,MagCat1 = [],[]
	MagCalc2,MagCat2 = [],[]
	Dec = []
	Imagens = os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')
	Imagem = Imagens[7]
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			Dec.append(DecCalc)
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			
			#5*4.7
			if dist > Threshold/3600:
				MagCalc1.append(float(Object['MagInst']))
				MagCat1.append(float(Object['MagReal']))
			else:
				MagCalc2.append(float(Object['MagInst']))
				MagCat2.append(float(Object['MagReal']))


	plt.plot(MagCalc1,MagCat1, 'or')
	plt.plot(MagCalc2,MagCat2, 'ob')
	plt.xlabel('Mag Instrumental')
	plt.ylabel('Mag Catalogo')
	plt.title(f"Threshold: {Threshold}")

if DesenhoFotometriaUnica3:
	MagCalc,MagCat,Distance = [],[],[]
	X,Y = [],[]
	Imagens = os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')
	Imagem = Imagens[26]
	print(Imagem)
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Distance.append(3600*dist)
			
	EscalaCores = np.log(Distance)

	plt.scatter(MagCalc,MagCat,cmap='hot',c=EscalaCores)
	plt.xlabel('Mag Instrumental',fontsize=fontsize)
	plt.ylabel('Mag Catalogo',fontsize=fontsize)
	cbar = plt.colorbar(ticks=[])
	cbar.set_label('Distance (arcsec)')
	
	fig,axes = plt.subplots(1,1)
	vmin,vmax = None,None
	Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Imagem.split('.')[0]+'.fit'}")
	raw = Data.data
	axes.imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
	axes.set_axis_off()
	for i in range(len(X)):
		if EscalaCores[i] > 2.5:
			c1 = plt.Circle((X[i],Y[i]),0.5, color = 'red', linewidth=2,fill=False)
			
		else:
			c1 = plt.Circle((X[i],Y[i]),0.5, color = 'green', linewidth=2,fill=False)
		axes.add_patch(c1)

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

if FotometriaUnicaDeitada:
	# Sets 1 e 3
	MagCalc,MagCat,Distance = [],[],[]
	X,Y = [],[]
	# ~ Images = os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')
	Imagem = Images[1]
	print(Imagem)
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem.split('.')[0]+'.json'}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Distance.append(3600*dist)
	
	EscalaCores = np.log(Distance)
	plt.scatter(MagCalc,np.array(MagCat)-np.array(MagCalc),cmap='hot',c=EscalaCores)
	plt.xlabel('Mag Instrumental',fontsize=fontsize)
	plt.ylabel('Mag Catalogo - Mag Instrumental',fontsize=fontsize)
	# ~ plt.xlim(-14,-7)
	# ~ plt.ylim(19,26)
	plt.tick_params(axis='x', labelsize=fontsize)
	plt.tick_params(axis='y', labelsize=fontsize)
	cbar = plt.colorbar(ticks=[])
	cbar.set_label('Distance',fontsize=fontsize)
	
if FotometriaUnicaDeitada2:
	Percentil = 0.70
	MagCalc,MagCat,Distance = [],[],[]
	X,Y = [],[]
	# ~ Images = os.listdir('/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles')
	Imagem = Images[1]
	print(Imagem)
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem.split('.')[0]+'.json'}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Distance.append(3600*dist)
	
	
	Distance2 = Distance.copy()
	Distance2.sort()
	Distance2 = Distance2[0:int(Percentil*len(Distance2))]
	maxDistance = Distance2[-1]
	mask = Distance < maxDistance
	
	print(len(Distance))
	print(len(mask))
	Distance 	= np.array(Distance)[mask]
	MagCat		= np.array(MagCat)[mask]
	MagCalc 	= np.array(MagCalc)[mask]	
	EscalaCores = np.log(Distance)

	plt.scatter(MagCalc,MagCat-MagCalc,cmap='hot',c=EscalaCores)
	plt.xlabel('Mag Instrumental',fontsize=fontsize)
	plt.ylabel('Mag Catalogo - Mag Instrumental',fontsize=fontsize)
	plt.xlim(-14,-7)
	plt.ylim(19,26)
	plt.tick_params(axis='x', labelsize=fontsize)
	plt.tick_params(axis='y', labelsize=fontsize)
	cbar = plt.colorbar(ticks=[])
	cbar.set_label('Distance',fontsize=fontsize)

if FotometriaUnicaDeitada3:
	Percentil = 0.40
	fig,axes = plt.subplots(2,2)
	ax = axes.ravel()
	Imagem1 = 'PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit' #Muitas
	Imagem2 = 'PAER_SURV_4_260771_230228_042347+015037_scan019_01.fit' #Poucas
	
	MagCalc,MagCat,Distance = [],[],[]
	X,Y = [],[]
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem1.split('.')[0]+'.json'}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Distance.append(3600*dist)
	
	
	Distance2 = Distance.copy()
	Distance2.sort()
	Distance2 = Distance2[0:int(Percentil*len(Distance2))]
	maxDistance = Distance2[-1]
	mask = Distance < maxDistance
	
	MagCat,MagCalc,Distance = np.array(MagCat),np.array(MagCalc),np.array(Distance)
	Distance2	= Distance[mask]
	MagCat2		= MagCat[mask]
	MagCalc2 	= MagCalc[mask]	
	EscalaCores = np.log(Distance)
	EscalaCores2 = np.log(Distance2)
	min_val=min(EscalaCores) 
	max_val=max(EscalaCores)

	sc0 = ax[0].scatter(MagCalc,MagCat-MagCalc,cmap='hot', vmin=min_val, vmax=max_val,c=EscalaCores)
	sc1 = ax[1].scatter(MagCalc2,MagCat2-MagCalc2,cmap='hot', vmin=min_val, vmax=max_val,c=EscalaCores2)
	fig.colorbar(sc0, ax=ax[0],ticks=[])
	fig.delaxes(fig.axes[4])
	fig.colorbar(sc1, ax=ax[1],ticks=[],label='Distance')
	ax[0].set_xlim(-14,-7)
	ax[0].set_ylim(19,26)
	ax[1].set_xlim(-14,-7)
	ax[1].set_ylim(19,26)
	
	with open(f'/home/luis/Desktop/Codigos/PipelineTrabalhar/{Imagem1}','a') as tempos:
		for i in range(len(MagCalc2)):
			linha = f'{MagCalc2[i]:8} {MagCat2[i]}\n'
			tempos.write(linha)
	
	MagCalc,MagCat,Distance = [],[],[]
	X,Y = [],[]
	with open(f"/home/luis/Desktop/Codigos/PipelineTrabalhar/Files/Images/2022-03-20/OutputFiles/{Imagem2.split('.')[0]+'.json'}", 'r') as stream:
		ImageData = json.load(stream)
	for Object in ImageData:
		if Object['Match']:
			RaCalc  = float(Object['RaCalculated'])
			DecCalc = float(Object['DecCalculated'])
			RaCat   = float(Object['RaCatalogue'])
			DecCat  = float(Object['DecCatalogue'])
			dist = SPADEFunctions.AngleDifferencePoints([RaCalc,DecCalc],[RaCat,DecCat],'normal')
			MagCalc.append(float(Object['MagInst']))
			MagCat.append(float(Object['MagReal']))
			X.append(float(Object['Xdetect']))
			Y.append(float(Object['Ydetect']))
			Distance.append(3600*dist)
	
	
	Distance2 = Distance.copy()
	Distance2.sort()
	Distance2 = Distance2[0:int(Percentil*len(Distance2))]
	maxDistance = Distance2[-1]
	mask = Distance < maxDistance
	
	MagCat,MagCalc,Distance = np.array(MagCat),np.array(MagCalc),np.array(Distance)
	Distance2	= Distance[mask]
	MagCat2		= MagCat[mask]
	MagCalc2 	= MagCalc[mask]	
	EscalaCores = np.log(Distance)
	EscalaCores2 = np.log(Distance2)
	min_val=min(EscalaCores) 
	max_val=max(EscalaCores)
	
	sc0 = ax[2].scatter(MagCalc,MagCat-MagCalc,cmap='hot', vmin=min_val, vmax=max_val,c=EscalaCores)
	sc1 = ax[3].scatter(MagCalc2,MagCat2-MagCalc2,cmap='hot', vmin=min_val, vmax=max_val,c=EscalaCores2)
	fig.colorbar(sc0, ax=ax[2],ticks=[])
	fig.delaxes(fig.axes[5])
	fig.colorbar(sc1, ax=ax[3],ticks=[],label='Distance')
	ax[2].set_xlim(-14,-7)
	ax[2].set_ylim(19,26)
	ax[3].set_xlim(-14,-7)
	ax[3].set_ylim(19,26)
	
	with open(f'/home/luis/Desktop/Codigos/PipelineTrabalhar/{Imagem2}','a') as tempos:
		for i in range(len(MagCalc2)):
			linha = f'{MagCalc2[i]:8} {MagCat2[i]}\n'
			tempos.write(linha)
			
	ax[0].set_xticks([])
	ax[1].set_xticks([])
	# ~ ax[1].set_yticks([])
	# ~ ax[3].set_yticks([])
	for ax in fig.axes:
		ax.tick_params(axis='both', which='both', labelsize=ticksize)  # Adjust the labelsize as desired
	fig.supxlabel('Mag Instrumental',fontsize=fontsize)
	fig.supylabel('Mag Catalogo - Mag Instrumental',fontsize=fontsize)
	
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

if DesenhoImage:
	Image = 'PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit'
	fig,axes = plt.subplots(1,2)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1= CCDData.read(f"{DayDirectory}/ReducedImages/{Image1}", unit = 'adu')
	raw1 	= Data1.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()

if DesenhoImages:
	Image1 = 'PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit'
	Image2 = 'PAER_SURV_4_260771_230228_042347+015037_scan019_01.fit'
	fig,axes = plt.subplots(1,2)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1,Data2 = CCDData.read(f"{DayDirectory}/ReducedImages/{Image1}", unit = 'adu'),CCDData.read(f"{DayDirectory}/ReducedImages/{Image2}", unit = 'adu')
	raw1,raw2 	= Data1.data,Data2.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()
	ax[1].imshow(raw2,cmap='gray',vmin=vmin or np.nanpercentile(raw2,1),vmax=vmax or np.nanpercentile(raw2, 99))
	ax[1].set_axis_off()
	
	fontprops = fm.FontProperties(size=8)
	Narcsec = 150
	
	scalebar = AnchoredSizeBar(ax[0].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[0].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[1].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[1].add_artist(scalebar)
	
if DesenhoImages2:
	Image1 = 'PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit'
	Image2 = 'PAER_SURV_4_260717_230228_042347+065500_scan014_01.fit'
	Image3 = 'PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit'
	Image4 = 'PAER_SURV_4_260777_230228_080700-074353_scan101_01.fit'
	Image5 = 'PAER_SURV_4_260783_230228_080700+045703_scan101_01.fit'
	
	
	
	fig,axes = plt.subplots(2,3)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1,Data2,Data3,Data4,Data5 = CCDData.read(f"{DayDirectory}/ReducedImages/{Image1}", unit = 'adu'), CCDData.read(f"{DayDirectory}/ReducedImages/{Image2}", unit = 'adu'), CCDData.read(f"{DayDirectory}/ReducedImages/{Image3}", unit = 'adu'), CCDData.read(f"{DayDirectory}/ReducedImages/{Image4}", unit = 'adu'), CCDData.read(f"{DayDirectory}/ReducedImages/{Image5}", unit = 'adu')
	raw1,raw2,raw3,raw4,raw5 	= Data1.data,Data2.data,Data3.data,Data4.data,Data5.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()
	ax[1].imshow(raw2,cmap='gray',vmin=vmin or np.nanpercentile(raw2,1),vmax=vmax or np.nanpercentile(raw2, 99))
	ax[1].set_axis_off()
	ax[2].imshow(raw3,cmap='gray',vmin=vmin or np.nanpercentile(raw3,1),vmax=vmax or np.nanpercentile(raw3, 99))
	ax[2].set_axis_off()
	ax[3].imshow(raw4,cmap='gray',vmin=vmin or np.nanpercentile(raw4,1),vmax=vmax or np.nanpercentile(raw4, 99))
	ax[3].set_axis_off()
	ax[4].imshow(raw5,cmap='gray',vmin=vmin or np.nanpercentile(raw5,1),vmax=vmax or np.nanpercentile(raw5, 99))
	ax[4].set_axis_off()
	
if DesenhoSet:
	fig,axes = plt.subplots(2,2,sharex=True,sharey=True)
	ax = axes.ravel()
	vmin,vmax = None,None
	for i,Image in enumerate(Images):
		Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Image}")
		raw = Data.data
		ax[i].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
		ax[i].set_axis_off()
	
	fontprops = fm.FontProperties(size=8)
	Narcsec=100
	scalebar = AnchoredSizeBar(ax[0].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[0].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[1].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[1].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[2].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[2].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[3].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[3].add_artist(scalebar)

if DesenhoEncontrados:
	for Images in Sets:
		A = time.time()
		dic = SPADEFunctions.CompareImagesJsonSkyCoord(DayDirectory,Images,scale = 4.7)
		dic_unicas = SPADEFunctions.MakeDicUnique(dic)
		objetos = SPADEFunctions.SeeVariation(dic_unicas,np.radians(float(config['VariationInImage'])/3600))
		objetos = SPADEFunctions.CleanObjects(objetos,config)
		Tempo = time.time() - A
		with open('/home/luis/Desktop/Codigos/PipelineTrabalhar/Tempos','a') as tempos:
			linha = f'{Tempo:8}\n'
			tempos.write(linha)
		fig,axes = plt.subplots(1,2)
		ax = axes.ravel()
		vmin,vmax = None,None
		Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Images[0]}", unit = 'adu')
		raw = Data.data
		ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,1),vmax=vmax or np.nanpercentile(raw, 99))


		for obj in objetos:
			for Number in objetos[obj]['Number']:
				data1 = dic[Number]
				image,X,Y = int(data1['image'][0])-1, data1['Xdetect'][0], data1['Ydetect'][0]
				c1 = plt.Circle((X,Y),size, color = color[image], linewidth=2,fill=False)
				ax[1].add_patch(c1)
		
		TLEOutput = f"{DayDirectory}/TLE/{Images[0][:-6]}Output.json"
		with open(TLEOutput, 'r') as stream:
			ImageData = json.load(stream)
		for obj in ImageData:
			X,Y,Name = obj['X (pixel) '],obj['Y (pixel) '], obj['Sat Catalog Number']
			c1 = plt.Circle((X,Y),size/2, color = 'orange', linewidth=1,fill=False)
			ax[1].text(X,Y,Name)
			ax[1].add_patch(c1)
			
		for image in dic_unicas:
			for i in range(len(dic_unicas[image]['X'])):
				X,Y = dic_unicas[image]['X'][i],dic_unicas[image]['Y'][i]
				c1 = plt.Circle((X,Y),size, color = color[image-1], linewidth=0.4,fill=False)
				ax[0].add_patch(c1)
		ax[1].set_title(f'Image: {Images[0]} ', fontsize=16)

if DesenhoMaster:
	MDark = 'MasterDark020.fits'
	MFlat = 'MasterFlatB.fits'
	fig,axes = plt.subplots(1,2)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1,Data2 = CCDData.read(f"{DayDirectory}/Masters/{MDark}", unit = 'adu'),CCDData.read(f"{DayDirectory}/Masters/{MFlat}", unit = 'adu')
	raw1,raw2 	= Data1.data,Data2.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()
	ax[1].imshow(raw2,cmap='gray',vmin=vmin or np.nanpercentile(raw2,1),vmax=vmax or np.nanpercentile(raw2, 99))
	ax[1].set_axis_off()

if DesenhoReduzida:
	Imagem = '1037-293B_001.fits'
	fig,axes = plt.subplots(1,2)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1,Data2 = CCDData.read(f"{DayDirectory}/ProcessedImages/{Imagem}", unit = 'adu'),CCDData.read(f"{DayDirectory}/ReducedImages/{Imagem}", unit = 'adu')
	raw1,raw2 	= Data1.data,Data2.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()
	ax[1].imshow(raw2,cmap='gray',vmin=vmin or np.nanpercentile(raw2,1),vmax=vmax or np.nanpercentile(raw2, 99))
	ax[1].set_axis_off()

if DesenhoVariosSigmas:
	Sigma2 = '/home/luis/Desktop/DadosSemCorrecaoCatalogo/Novos/2Sigma/'
	Sigma3 = '/home/luis/Desktop/DadosSemCorrecaoCatalogo/Novos/3Sigma/'
	Sigma5 = '/home/luis/Desktop/DadosSemCorrecaoCatalogo/Novos/5Sigma/'
	Imagem = 'PAER_SURV_4_260789_230228_080700-051142_scan102_01'
	
	fig,axes = plt.subplots(1,3,sharey=True,sharex=True)
	ax = axes.ravel()
	vmin,vmax = None,None
	Data1 = CCDData.read(f"{DayDirectory}/ProcessedImages/{Imagem+'.fit'}")
	raw1  = Data1.data
	ax[0].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[0].set_axis_off()
	ax[1].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[1].set_axis_off()
	ax[2].imshow(raw1,cmap='gray',vmin=vmin or np.nanpercentile(raw1,1),vmax=vmax or np.nanpercentile(raw1, 99))
	ax[2].set_axis_off()
	
	with open(Sigma2+Imagem+'.json', 'r') as stream:
		ImageData = json.load(stream)
		for Object in ImageData:
			X,Y = Object['Xdetect'],Object['Ydetect']
			c1 = plt.Circle((X,Y),size/2, color = 'red', linewidth=1,fill=False)
			ax[0].add_patch(c1)
			
	with open(Sigma3+Imagem+'.json', 'r') as stream:
		ImageData = json.load(stream)
		for Object in ImageData:
			X,Y = Object['Xdetect'],Object['Ydetect']
			c1 = plt.Circle((X,Y),size/2, color = 'red', linewidth=1,fill=False)
			ax[1].add_patch(c1)
			
	with open(Sigma5+Imagem+'.json', 'r') as stream:
		ImageData = json.load(stream)
		for Object in ImageData:
			X,Y = Object['Xdetect'],Object['Ydetect']
			c1 = plt.Circle((X,Y),size/2, color = 'red', linewidth=1,fill=False)
			ax[2].add_patch(c1)		

	fontprops = fm.FontProperties(size=8)
	Narcsec = 150
	
	scalebar = AnchoredSizeBar(ax[0].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[0].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[1].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[1].add_artist(scalebar)
	
	scalebar = AnchoredSizeBar(ax[2].transData,
                            Narcsec/4.7, f'({Narcsec} arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
	ax[2].add_artist(scalebar)
	
plt.show()	









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

