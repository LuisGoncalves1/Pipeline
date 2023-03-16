from Functions import SPADEFunctions
import os
import numpy as np
import math
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits
import json
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import EarthLocation
import datetime
import astropy.units as u
import requests
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84,utc
import statistics

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)

ficheiro='MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits'
tle_file='2022tle.json'
Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + ficheiro}", unit = 'adu') 
header=Data.header

#Testing fwhm 

def testing_fwhm(ficheiro,DayDirectory):
	fwhm=[]
	with open(f"{DayDirectory + '/OutputFiles/' + ficheiro}.json", 'r') as stream:
		jsonfile = json.load(stream)
		for i in jsonfile:
			fwhm.append(i['FWHM'])
	print(statistics.median(fwhm))
	plt.hist(fwhm,bins=1000)
	plt.xlim(0,20)
	plt.show()
	

def tle_txt_to_json(filename):
	
	items_2_store=[]
	tmp_dict={}
	with open(filename,'r') as fp:
		tle_data= fp.read().splitlines() 
		for i in tle_data:
			if i[0] != '1' and i[0] != '2':
				tmp_dict['satellite_name'] = i
			if i[0] == '1':
				tmp_dict['tle_1'] = i
			if i[0] == '2':
				tmp_dict['tle_2'] = i

			if "tle_1" in tmp_dict and "tle_2" in tmp_dict and "satellite_name" in tmp_dict:
				items_2_store.append(tmp_dict)
				tmp_dict={}
				
	with open(f"{filename[:-4]}.json", 'w') as f:
			json.dump(items_2_store, f, indent=3)


#convertion=SPADEFunctions.XYToRaDec(FileWc,FileXY,FileOutput)


def SatellitePosition(config,header,tle_file,DayDirectory,ficheiro,scale):
	'''
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	with open(f"{tle_file}", 'r') as stream:
		lista_ficheiro=json.load(stream)
		
	ts = load.timescale()
	
	Header=SPADEFunctions.ReadFilesHeader(DayDirectory,header,ficheiro)
	

	t = ts.utc(Header[0]['Year'],Header[0]['Month'],Header[0]['Day'],Header[0]['Hour'],Header[0]['Minute'],Header[0]['Second'])

	Alt_target=Header[0]['Alt']
	Az_target=Header[0]['Az']
	Dimensions=Header[0]['Dimensions']
	FOV=(Dimensions*scale)/3600
	
	Alt_window=[Alt_target-FOV/2, Alt_target+FOV/2]
	Az_window=[Az_target-FOV/2, Az_target+FOV/2]
		
	for i in lista_ficheiro:
		
		items={}
		line1=i['tle_1']
		line2=i['tle_2']
		nome=i['satellite_name']
	
		satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 
		bluffton = wgs84.latlon(float(config['Latitude']), float(config['Longitude']),float(config['Height']))
		
		difference = satellite - bluffton 
		topocentric = difference.at(t)
			
		alt, az, distance = topocentric.altaz()
		
		if alt.degrees >= Alt_window[0] and alt.degrees <= Alt_window[1] and az.degrees >= Az_window[0] and az.degrees <= Az_window[1]:
			print(f"Satellite {nome}, alt:{alt.degrees},az:{az.degrees}")
		

#SatellitePosition(config,header,tle_file,DayDirectory,ficheiro,scale=4.7)



def CreateRDLSFile(X,Y,FileName):
	''' 
	Function that changes the data from the RLDS file to a set of (Ra,Dec) coordinates.
	
	INPUTS:
	Ra 					-> array with Ra coordinates									(list)
	Dec 				-> array with Dec coordinates									(list)
	FileName 			-> Name of rdls file											(rdls)
	
	RETURNS:
	NONE
	
	FILES CHANGED:
	FileName			-> Changes Ra,Dec information to match with Ra,Dec given		(rdls)
	'''
	with fits.open(FileName, mode = 'update') as A:
		Dados = A[1].data
	
		x = np.array(X,Y,dtype = [('X_IMAGE', '>f4'),('Y_IMAGE', '>f4')])
		Dados = Dados[Dados['X_IMAGE'] > 2049] #elimina os dados menores 2049
		Dados = np.append(Dados,x)
		

		A.close()

X=np.array([1,1,2048,2048])
Y=np.array([1,2048,1,2048])
FileName='teste.cat'
CreateRDLSFile(X,Y,FileName)






def scattered_plot_parameters(lst):
	
	colors_detections=[]	
	for i in lst:
		if len(i)==1:
			black=(0.2,0.2,0.2)
			colors_detections.append(black)
		if len(i)==2:
			yellow=(1, 1, 0)
			colors_detections.append(yellow)
		if len(i)==3:
			blue=(0, 0.5, 1)
			colors_detections.append(blue)
		if len(i)==4:
			red=(1,0, 0)
			colors_detections.append(red)

	x=[]
	y=[]
	for i in lst:
		x.append(float(i[0]['RaCalculated']))

	for i in lst:
		y.append(float(i[0]['DecCalculated']))
		
	return x,y,colors_detections

def scatteredplot_all_detections(lst):
	
	alldect=scattered_plot_parameters(lst)
	
	
	plt.scatter(
	x=alldect[0],
	y=alldect[1],
	c=alldect[2],
	alpha=0.5, 
	)
	
	my_xticks=[]
	x=alldect[0]
	x.sort()
	print(x)
	for i in x:
		my_xticks.append(decdeg2dms(i))


	plt.xticks([x[0], x[-1]], [my_xticks[0], my_xticks[-1]],rotation=0,fontsize=8)
	
	my_yticks=[]
	y=alldect[1]
	y.sort()
	for i in y:
		my_yticks.append(decdeg2dms(i))
	plt.yticks([y[0], y[-1]], [my_yticks[0], my_yticks[-1]],rotation=0,fontsize=8)
	
	
	plt.title("All Objects")
	plt.ylabel('Dec (dms)')
	plt.xlabel('Ra (dms)')
	
	plt.show()
	
