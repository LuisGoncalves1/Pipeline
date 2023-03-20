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

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)

	
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

		# Define the field names and types for X and Y separately
		X_fields = [('X_IMAGE', '>f4')]
		Y_fields = [('Y_IMAGE', '>f4')]

		# Concatenate the field definitions
		fields = X_fields + Y_fields

		# Create the structured array using np.rec.fromarrays()
		x = np.rec.fromarrays([X, Y], dtype=fields)
		print(x)
		Dados = Dados[Dados['X_IMAGE'] > 2049] #elimina os dados menores 2049
		print(Dados)
		Dados = np.append(Dados,x)
		A.close()

X=[1,1,2048,2048]
Y=[1,2048,1,2048]
FileName='teste.cat'
CreateRDLSFile(X,Y,FileName)


def TLE_parameters(tle_catalogfile):
	
	data=[]
	earth_radius=6.3781*10**3 #in kms
	mu=3.986004418*10**5 #standard gravitacional parameter of the Earth in km³/s²
	geo_altitude=35786 #in kms
	
	with open(tle_catalogfile,'r') as fp:
		GP_elements= json.load(fp)
		
		for i in GP_elements:
			temp={}
			temp['tle_1']=i['tle_1']
			temp['tle_2']=i['tle_2']
			
			eccentricity=float('0.'+i['tle_2'][26:33])
			
			mean_motion=float(i['tle_2'][52:63]) #revolutions per day 

			convertion=((2*np.pi)/(24*3600))
			mean_motion=mean_motion*convertion #radians per second


			a=(mu**(1/3))/(mean_motion**(2/3)) # calculates the semi-major axis

			Perigee=(a*(1-eccentricity)-earth_radius)

			Apogee=(a*(1+eccentricity)-earth_radius) 

			if Perigee>geo_altitude-200 and Apogee<200+geo_altitude:
				data.append(temp)
	return data
									


ficheiros=['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit','PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit']


def SatellitePosition(config,DayDirectory,ficheiros):
	'''
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	
	tle_list=TLE_parameters('Full Catalog-20230228T0000.json')
	print(len(tle_list))
	
	ts = load.timescale()
	
	info=[]
	
	for image in ficheiros:
	
		Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + image}", unit = 'adu') 
		header=Data.header

		Header=SPADEFunctions.ReadFilesHeader(DayDirectory,header,image)
	
		info.append(Header)

	t0 = ts.utc(info[0]['Year'],info[0]['Month'],info[0]['Day'],info[0]['Hour'],info[0]['Minute'],info[0]['Second'])
	tf=  ts.utc(info[1]['Year'],info[1]['Month'],info[1]['Day'],info[1]['Hour'],info[1]['Minute'],info[1]['Second'])
	
	for i in tle_list:
			
		items={}
		line1=i['tle_1']
		line2=i['tle_2']
		nome='name'
	
		satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 

		bluffton = wgs84.latlon(float(config['Latitude']), float(config['Longitude']),float(config['Height']))
			
		difference = satellite - bluffton 
		for t in np.arange(t0,tf,0.00001):
			topocentric = difference.at(t)
				
			alt, az, distance = topocentric.altaz()
			print(f"alt:{alt.degrees},az:{az.degrees}")
			
			#if alt.degrees >= Alt_window[0] and alt.degrees <= Alt_window[1] and az.degrees >= Az_window[0] and az.degrees <= Az_window[1]:
				#print(f"Satellite {nome}, alt:{alt.degrees},az:{az.degrees}")



#SatellitePosition(config,DayDirectory,ficheiros)
