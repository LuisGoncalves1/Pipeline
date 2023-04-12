import os
import shutil
from Functions import SPADEFunctions
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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
TLE_file='Full Catalog-20230228T0000.json'

#X=[1,1,2048,2048]
#Y=[1,2048,1,2048]
	
def CreateCatFile_XY(X,Y,FileName):
	''' 
	Function that changes the data from the cat file to a set of (x,y) coordinates. 
	
	INPUTS:
	Ra 					-> array with x coordinates									(list)
	Dec 				-> array with y coordinates									(list)
	FileName 			-> Name of cat file											 (cat)
	
	RETURNS:
	NONE
	
	FILES CHANGED:
	FileName			-> Changes x,y information to match with x,y given		    (cat)
	'''
	with fits.open(FileName, mode = 'update') as A:

		# Define the field names and types for X and Y separately
		X_fields = [('X_IMAGE', '>f4')]
		Y_fields = [('Y_IMAGE', '>f4')]

		# Concatenate the field definitions
		fields = X_fields + Y_fields

		# Create the structured array using np.rec.fromarrays()
		x = np.rec.fromarrays([X, Y], dtype=fields)

		A[1].data[A[1].data['X_IMAGE'] > 2049] #eliminate all data from the file
		
		A[1].data = np.array(x) # Convert the list of tuples to a NumPy array

		A.close()
		
		
def CreateRDLSFile_RaDEC(Ra,Dec,FileName):
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

		# Define the field names and types for X and Y separately
		
		X_fields = [('RA', '<f8')]
		Y_fields = [('DEC','<f8')]

		# Concatenate the field definitions
		fields = X_fields + Y_fields

		# Create the structured array using np.rec.fromarrays()
		x = np.rec.fromarrays([Ra, Dec], dtype=fields)
		
		A[1].data=A[1].data[A[1].data['RA'] > 180] #eliminate all data from the file
		
		
		A[1].data = np.array(x) # Convert the list of tuples to a NumPy array

		A.close()


def TLE_parameters(TLE_file):
	''' 
	Function that takes the TLE catalog previously stored in a json file and filters the satellites to those select in GeoProtectedZonePlus
	
	INPUTS:
	TLE_file		 					-> Tle Catalog in a json file									(json)

	RETURNS:
	data  								> List of satellites that meet the criteria						(list)
	
	FILES CREATED:
	NONE

	'''
	
	data=[]
	earth_radius=6.3781*10**3 #in kms
	mu=3.986004418*10**5 #standard gravitacional parameter of the Earth in km³/s²
	geo_altitude=35786 #in kms
	
	with open(TLE_file,'r') as fp:
		GP_elements= json.load(fp)
		
		for i in GP_elements:
			temp={}
			temp['Satellite Number']=i['tle_2'][2:8]
			temp['tle_1']=i['tle_1']
			temp['tle_2']=i['tle_2']
			
			eccentricity=float('0.'+i['tle_2'][26:33])
			
			mean_motion=float(i['tle_2'][52:63]) #revolutions per day 

			convertion=((2*np.pi)/(24*3600))
			mean_motion=mean_motion*convertion #radians per second


			a=(mu**(1/3))/(mean_motion**(2/3)) # calculates the semi-major axis

			Perigee=(a*(1-eccentricity)-earth_radius)

			Apogee=(a*(1+eccentricity)-earth_radius) 
			
			
			#The following condition filters the tles from the catalog and stores only those in GeoProtectedZonePlus orbit
			if Perigee>geo_altitude-200 and Apogee<200+geo_altitude:  
				data.append(temp)
	return data
									

FIT_Files=['PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit',
'PAER_SURV_4_260789_230228_080700-051142_scan102_02.fit',
'PAER_SURV_4_260789_230228_080700-051142_scan102_03.fit',
'PAER_SURV_4_260789_230228_080700-051142_scan102_04.fit']

FileXY='XY_window.cat'

def AltAz_window(config,info,image,TLE_Directory,DayDirectory):
	
	''' 
	Function that returns the maximum and minimum value of Alt,Az coordinates of the FOV of the image in question
	
	INPUTS:
	config				-> Configuration Dictionary										(dic)
	info 				-> Image header information										(dic)
	image				-> Image that's being analysed									(str)
	TLE_Directory       -> Directory of the TLE information								(str)

	RETURNS:
	Alt_window, Az_window  		-> Two lists, containing the maximum and minimum value of Alt,Az					(list)
	
	FILES CREATED:
	f'{image.split(".")[0]}'    > Fit file with RaDec coordinates from XY_window.cat file	

	'''
	
	#Transforms the XY coodinates [(1,1);(1,2048);(2048,1);(2048,2048)] from XY_window.cat file to RaDec for each image.

	FileWc=f'{DayDirectory}/AstrometryFiles/{image.split(".")[0]}.wcs'
	FileOutput= TLE_Directory + f'{image.split(".")[0]}'
	SPADEFunctions.XYToRaDec(FileWc,FileXY,FileOutput)

	
	location=EarthLocation(
			lat=float(config['Latitude'])*u.deg,
			lon=float(config['Longitude'])*u.deg,
			height=float(config['Height'])*u.m)
	
	Alt=[]
	Az=[]

	time = datetime.datetime(info['Year'],info['Month'],info['Day'],info['Hour'],info['Minute'],info['Second'],info['Microsecond'])

	A=fits.open(f'{TLE_Directory + image.split(".")[0]}')
		
	declination=A[1].data['DEC']
	right_ancention=A[1].data['RA']

	#Ra,Dec coodrinates are tranformed to Alt,Az coordinates 
	for i in range(0,len(declination)):
				
		dec=(declination[i])*u.deg
		ra=(right_ancention[i])*u.deg
		stars = SkyCoord(ra, dec, frame='icrs')
	
		stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))
		sky_coords = [stars_altaz.alt.degree,stars_altaz.az.degree]
	
		Alt.append(sky_coords[0])
		Az.append(sky_coords[1])
			
	Alt_window=[min(Alt),max(Alt)]
	Az_window=[min(Az),max(Az)]

	
	return Alt_window, Az_window

		
def SatellitePosition(config,DayDirectory,FileXY,FIT_Files,TLE_file):
	
	''' 
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	TLE_Directory = f"{PipelineDirectory}/Files/TLE/"
	tle_list=TLE_parameters(TLE_file)
	ts = load.timescale()
	
	satellites=[]
	
	for num in range(0,len(FIT_Files)):
		
	
		image=FIT_Files[num]
		Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + image}", unit = 'adu') 
		header=Data.header
		info=SPADEFunctions.ReadFilesHeader(DayDirectory,header,image)
		t = ts.utc(info['Year'],info['Month'],info['Day'],info['Hour'],info['Minute'],info['Second'])
	
		Alt,Az=AltAz_window(config,info,image,TLE_Directory,DayDirectory,PipelineDirectory)
		Ra=[]
		Dec=[]
	
		for i in tle_list:
			
			items={}
			line1=i['tle_1']
			line2=i['tle_2']
			nome=i['Satellite Number']
	
			satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 

			bluffton = wgs84.latlon(float(config['Latitude']), float(config['Longitude']),float(config['Height']))
			
			difference = satellite - bluffton 
		
			topocentric = difference.at(t)
				 
			alt, az, distance = topocentric.altaz()
			ra, dec, distance = topocentric.radec() 
			
			if alt.degrees >= Alt[0] and alt.degrees <= Alt[1] and az.degrees >= Az[0] and az.degrees <= Az[1]:
				item={}
				item['Satellite']=nome
				item['Ra']=math.degrees(ra.radians)
				item['Dec']=math.degrees(dec.radians)
				satellites.append(item)		
				
				
				Ra.append(math.degrees(ra.radians))
				Dec.append(math.degrees(dec.radians))
		

		FileName=TLE_Directory +f'Sat_RaDec_image{num+1}_dec=-5.rdls' 
		CreateRDLSFile_RaDEC(Ra,Dec,FileName)

		
		FileWc=f'{DayDirectory}/AstrometryFiles/{image.split(".")[0]}.wcs'
		FileOutput=TLE_Directory + f'{image.split(".")[0]}_radec'
		SPADEFunctions.RaDecToXY(FileWc,FileName,FileOutput)
		

	return satellites
				

satellites=SatellitePosition(config,DayDirectory,FileXY,FIT_Files,TLE_file,PipelineDirectory)

	
#image='PAER_SURV_4_260789_230228_080700-051142_scan102_01.fit'
		
def visualization_detections_and_satellites(TLE_Directory,image):
	
	'''
	Function created to plot the image being analysed and draw the detections on top of it, according to json output file
	
	INPUTS:
	all_files			->   A list of every sextrator file in str type
	DayDirectory		->   Directory 
	
	FILES CREATED:
	NONE
	'''
	
	with open(f"{TLE_Directory + '/OutputFiles/' + image.split('.')[0]}.json", 'r') as stream:
		jsonfile = json.load(stream)
		
	A=fits.open(f'{image.split(".")[0]}_radec')
		
	x=A[1].data['X']
	y=A[1].data['Y']

	
	fig = plt.figure()
	ax = fig.add_subplot(111) 
	vmin,vmax = None,None


	Imagem = CCDData.read(f"{TLE_Directory  + '/ProcessedImages/' + image}", unit = 'adu') 
	plt.imshow(Imagem,cmap='gray',vmin=vmin or np.nanpercentile(Imagem,3),vmax=vmax or np.nanpercentile(Imagem, 99))
	#plt.imshow(Imagem,cmap='gray')
	plt.title(f"{image.split('.')[0]}")


	#for i in jsonfile:
	#	c1 = plt.Circle((float(i['Xdetect']),float(i['Ydetect'])), color = 'red', linewidth=2,fill=False)
	#	ax.add_patch(c1)
	
	for i in range(0,len(x)):
		c2 = plt.Circle((x[i],y[i]),4, color = 'blue', linewidth=2,fill=False)
		ax.add_patch(c2)
		
	# fontprops = fm.FontProperties(size=8)

	# scalebar = AnchoredSizeBar(ax.transData,
                            # 75/4.7, '(75 arcsec)', 
                            # 'lower center', 
                           # pad=0.1,
                           # color='white',
                           # frameon=False,
                           # size_vertical=1,
                           # fontproperties=fontprops)

	# ax.add_artist(scalebar)
	
	plt.show()
	
	return True


#visualization_detections_and_satellites(DayDirectory,image)
 




