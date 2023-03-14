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

ficheiro='MAPA_SURV_4_085283_220417_114004+060023_scan047_01'
tle_file='2022tle.json'
#Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + ficheiro}", unit = 'adu') 
#header=Data.header

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
testing_fwhm(ficheiro,DayDirectory)		




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

#TER EM ATENÇÃO AOS VALORES DE LAT E LONG QUE ESTOU A USAR NO CONFIG FILE QUANDO USO ESTA FUNÇÃO

def SatellitePosition(config,header,tle_file,DayDirectory,ficheiro):
	'''
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	with open(f"{tle_file}", 'r') as stream:
		lista_ficheiro=json.load(stream)
		
	ts = load.timescale()
	
	Header=SPADEFunctions.ReadFilesHeader(DayDirectory,header,ficheiro)
	
	second=str(Header[0]['Second'])+'.'+str(Header[0]['Microsecond'])
	t = ts.utc(Header[0]['Year'],Header[0]['Month'],Header[0]['Day'],Header[0]['Hour'],Header[0]['Minute'],float(second))
		
	Dimensions=Header[0]['Dimensions']
	x0=1; y0=1; xf=Dimensions; yf=Dimensions

	corner1=[x0,y0]
	corner2=[xf,y0]
	corner3=[x0,yf]
	corner4=[xf,yf]

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
		
		
		# if alt.degrees > Alt_window[0] and alt.degrees< Alt_window[1] and az.degrees > Az_window[0] and az.degrees< Az_window[1]:
			# print(f"Satellite {nome}, alt:{alt.degrees},az:{az.degrees}")
		# else:
			# print(f"Satellite {nome} is not observable in this conditions")
			
			


#SatellitePosition(config,header,tle_file,DayDirectory,ficheiro)





