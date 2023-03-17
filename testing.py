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
ficheiro='MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits'
tle_file='2022tle.json'
# Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + ficheiro}", unit = 'adu') 
# header=Data.header
	
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
		
		Dados = Dados[Dados['X_IMAGE'] > 2049] #elimina os dados menores 2049
		
		Dados=np.concatenate([Dados,x])
		A.close()

# X=[1,1,2048,2048]
# Y=[1,2048,1,2048]
# FileName='teste.cat'
# CreateRDLSFile(X,Y,FileName)


def TLE_parameters():
	data=[]
	with open('Full Catalog-20230228T0000.tle','r') as fp:
		GP_elements= fp.read().splitlines() 
		print(len(GP_elements))
		for i in GP_elements:
			if i[0] == '2':
				
				temp={}
				temp['Catalog Number']=i[2:7]
	
				eccentricity=float('0.'+i[26:33])
				temp['Eccentricity']=eccentricity  

				mean_motion=float(i[52:63]) #revolutions per day 
	
				convertion=((2*np.pi)/(24*3600))
				mean_motion=mean_motion*convertion #radians per second

				temp['Mean Motion']=mean_motion 

				mu=3.986004418*10**5 #in km³/s², standard gravitacional parameter of the Earth
				a=(mu**(1/3))/(mean_motion**(2/3)) # calculates the semi-major axis

				temp['Semi-major Axis']=a #in kms
	
				Perigee=(a*(1-eccentricity)-6378.1)
				temp['Perigee']=Perigee #in kms
	
				Apogee=(a*(1+eccentricity)-6378.1)
				temp['Apogee']=Apogee #in kms
	
				if Perigee<200 and Apogee>200:
					temp['Perigee']=Perigee #in kms
					temp['Apogee']=Apogee #in kms
					data.append(temp)
									
	print(len(data))
	
TLE_parameters()

#verify calculations 

def verify_calculations():
	
	i='2 41729   0.0477  98.5584 0002767 195.8549 260.0049  1.00271002 24151'
	temp={}
	temp['Catalog Number']=i[2:7]
	
	eccentricity=float('0.'+i[26:33])
	temp['Eccentricity']=eccentricity  

	mean_motion=float(i[52:63]) #revolutions per day 
	
	convertion=((2*np.pi)/(24*3600))
	mean_motion=mean_motion*convertion #radians per second

	temp['Mean Motion']=mean_motion 

	mu=3.986004418*10**5 #in km³/s², standard gravitacional parameter of the Earth
	a=(mu)**(1/3)/(mean_motion**(2/3)) # calculates the semi-major axis
	
	
	temp['Semi-major Axis']=a #in kms
	print(a*(1-eccentricity))
	Perigee=(a*(1-eccentricity)-6378.1)*10**(-3)
	temp['Perigee']=Perigee #in kms
	
	Apogee=(a*(1+eccentricity)-6378.1)*10**(-3)
	temp['Apogee']=Apogee #in kms
	
	print(temp)


#verify_calculations()
def backwards_calculation():
	
	i='2 27952   7.5703  29.0833 7144309   6.0011 182.2870  2.25976661160559'
	temp={}
	temp['Catalog Number']=i[2:7]
	
	eccentricity=float('0.'+i[26:33])
	Apogee=35.676
	
	a=Apogee/(1+eccentricity)
	temp['Semi-major Axis']=a 
	
	miu=3.986004418*10**5 #in km³/s², standard gravitacional parameter of the Earth
	
	n=np.sqrt(miu/(a**3))
	temp['Mean Motion']=n
	
	convertion=n/(2.25976661)
	print(temp)
	print(convertion)
	
	
#backwards_calculation()
