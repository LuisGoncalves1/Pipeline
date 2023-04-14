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
TLE_Directory = f"{DayDirectory}/TLE"
FileXY=f'{TLE_Directory}/XY_window.cat'

def ArrangeTLEFiles(DayDirectory,TLE_Directory,FIT_Files):

		for i in range(0,len(FIT_Files)):
			image=FIT_Files[i]
			
			#The output file for AltAz_window
			src_file_path = f'{DayDirectory}/AstrometryFiles/{image.split(".")[0]}'
			dst_file_path = f'{TLE_Directory}/{image.split(".")[0]}'
			shutil.copy(src_file_path, dst_file_path)
			new_dst_file_name = f'{image.split(".")[0]}_RaDec_Window'
			os.rename(dst_file_path, os.path.join(os.path.dirname(dst_file_path), new_dst_file_name))
			
		#The rdls file for satellite RaDec
		src_file_path = f'{DayDirectory}/AstrometryFiles/{FIT_Files[0].split(".")[0]}.rdls'
		dst_file_path = f'{TLE_Directory}/{image.split(".")[0]}.rdls'
		shutil.copy(src_file_path, dst_file_path)
		new_dst_file_name = f'{FIT_Files[0].split(".")[0][:-3]}_RaDec.rdls'
		os.rename(dst_file_path, os.path.join(os.path.dirname(dst_file_path), new_dst_file_name))
		
		#The xy file for satellite
		src_file_path = f'{DayDirectory}/AstrometryFiles/{FIT_Files[0].split(".")[0]}'
		dst_file_path = f'{TLE_Directory}/{image.split(".")[0]}'
		shutil.copy(src_file_path, dst_file_path)
		new_dst_file_name = f'{FIT_Files[0].split(".")[0][:-3]}_XY'
		os.rename(dst_file_path, os.path.join(os.path.dirname(dst_file_path), new_dst_file_name))
		
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

def TLE_parameters(TLE_Directory,TLE_file):
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
	
	with open(f'{TLE_Directory}/{TLE_file}','r') as fp:
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
			#if Perigee>geo_altitude-200 and Apogee<200+geo_altitude:  
			#	data.append(temp)
			data.append(temp)
	return data
									
def AltAz_window(config,info,image,TLE_Directory,DayDirectory,FileXY):
	
	''' 
	Function that returns the maximum and minimum value of Alt,Az coordinates of the FOV of each image expressed by FIT_Files.
	It's assumed that the image maximum value of Alt,Az corresponds to the positions [(1,1);(1,2048);(2048,1);(2048,2048)] in pixels.
	Therefore, this function creates a Fit file with RaDec coordinates from XY_window.cat file for each image, and then converts it to Alt,Az
	coordinates. 
	
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
	
	Alt=[]
	Az=[]
	
	#Transforms the XY coodinates [(1,1);(1,2048);(2048,1);(2048,2048)] from XY_window.cat file to RaDec for each image.
	FileWc=f'{DayDirectory}/AstrometryFiles/{image.split(".")[0]}.wcs'
	FileOutput=f'{TLE_Directory}/{image.split(".")[0]}_RaDec_Window'
	SPADEFunctions.XYToRaDec(FileWc,FileXY,FileOutput)
	
	#Opens fit file created with RaDec coordinates
	A=fits.open(f'{TLE_Directory}/{image.split(".")[0]}_RaDec_Window')
		
	declination=A[1].data['DEC']
	right_ancention=A[1].data['RA']
	
	#Location and date information is necessary to convert RaDec to Alt,Az
	location=EarthLocation(
			lat=float(config['Latitude'])*u.deg,
			lon=float(config['Longitude'])*u.deg,
			height=float(config['Height'])*u.m)
	
	time = datetime.datetime(info['Year'],info['Month'],info['Day'],info['Hour'],info['Minute'],info['Second'],info['Microsecond'])

	#Ra,Dec coordinates are tranformed to Alt,Az coordinates 
	for i in range(0,len(declination)):
				
		dec=(declination[i])*u.deg
		ra=(right_ancention[i])*u.deg
		stars = SkyCoord(ra, dec, frame='icrs')
	
		stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))
		sky_coords = [stars_altaz.alt.degree,stars_altaz.az.degree]
	
		Alt.append(sky_coords[0])
		Az.append(sky_coords[1])
		
	os.remove(f'{TLE_Directory}/{image.split(".")[0]}_RaDec_Window')
			
	Alt_window=[min(Alt),max(Alt)]
	Az_window=[min(Az),max(Az)]

	return Alt_window, Az_window
	
def SatellitePosition(config,DayDirectory,FileXY,FIT_Files,header_data,TLE_file,TLE_Directory):
	''' 
	This function loads TLE information from json file and returns a list of the satellites that were present at the time and place the images were obtained.
	In order to represent the satelites position in XY coordinates, cat files are created with RaDec coordinates and XY coordinates of the selected 
	satellites.
	
	INPUTS:
	config				-> Configuration Dictionary												(dic)
	DayDirectory 		-> Directory of present day												(str)
	FileXY 				-> Cat file with XYcoordinates [(1,1);(1,2048);(2048,1);(2048,2048)]	(str)
	FIT_Files 			-> List of Fit files of analysed images									(lst)
	TLE_file			-> Json file with TLE information										(str)
	
	RETURNS:
	satellites  		-> List of satellite names and Ra Dec coordinates that were present at a given time and date	(list)
	
	'''
	
	ArrangeTLEFiles(DayDirectory,TLE_Directory,FIT_Files)
	
	tle_list=TLE_parameters(TLE_Directory,TLE_file)
	ts = load.timescale()
	
	satellites=[]
	Ra=[]
	Dec=[]

	for num in range(0,len(FIT_Files)):
		
		image=FIT_Files[num]
		info=header_data[num]
		t = ts.utc(info['Year'],info['Month'],info['Day'],info['Hour'],info['Minute'],info['Second'])
		Alt,Az=AltAz_window(config,info,image,TLE_Directory,DayDirectory,FileXY)
	
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
				item['Image']=FIT_Files[num]
				item['Sat Catalog Number']=nome
				item['Ra (degrees)']=math.degrees(ra.radians)
				item['Dec (degrees)']=math.degrees(dec.radians)
				item['Alt (degrees) ']=alt.degrees
				item['Az (degrees) ']=az.degrees
				satellites.append(item)		
				
				
				Ra.append(math.degrees(ra.radians))
				Dec.append(math.degrees(dec.radians))
		
	FileName=f'{TLE_Directory}/{FIT_Files[0].split(".")[0][:-3]}_RaDec.rdls' 
	CreateRDLSFile_RaDEC(Ra,Dec,FileName)

	
	FileWc=f'{DayDirectory}/AstrometryFiles/{FIT_Files[0].split(".")[0]}.wcs'
	FileOutput=f'{TLE_Directory}/{FIT_Files[0].split(".")[0][:-3]}_XY'
	result, status = SPADEFunctions.RaDecToXY(FileWc,FileName,FileOutput)
	
	return satellites
	
def Write_TLE_Output_Json(FIT_Files,Satellite,TLE_Directory):
	
	A=fits.open(f'{TLE_Directory}/{FIT_Files[0].split(".")[0][:-3]}_XY')
	x=A[1].data['X']
	y=A[1].data['Y']
	data=[]
	for i in range(0,len(x)):
		if x[i]<0 or y[i]<0:
			continue
		else:
			Satellite[i]['X (pixel) ']=x[i]
			Satellite[i]['Y (pixel) ']=y[i]
			data.append(Satellite[i])
			
			
	with open(f'{TLE_Directory}/{FIT_Files[0].split(".")[0][:-3]}_Ouput.json', 'w') as f:
		json.dump(data, f, indent=3)
			
	return Satellite
	
def Drawing_Satellites_TLE_Position(DayDirectory,CatFilesDirectory,image,sat_xy_coordinates):
	
	'''
	Function created to plot the image being analysed and draw the position of the satellites on top of it, according to the files created in SatellitePosition
	function with TLES
	
	INPUTS:
	DayDirectory 		-> Directory of present DayDirectory							(str)	
	image				-> Image that's being analysed									(str)
	sat_xy_coordinates  -> Fit file with satellites XY coordinates					    (str)
	
	FILES CREATED:
	NONE
	'''
	
	
	A=fits.open(f'{CatFilesDirectory}/{sat_xy_coordinates}')
		
	x=A[1].data['X']
	y=A[1].data['Y']

	
	fig = plt.figure()
	ax = fig.add_subplot(111) 
	vmin,vmax = None,None


	Imagem = CCDData.read(f"{DayDirectory}/ProcessedImages/{image}", unit = 'adu') 
	plt.imshow(Imagem,cmap='gray',vmin=vmin or np.nanpercentile(Imagem,3),vmax=vmax or np.nanpercentile(Imagem, 99))
	plt.title(f"{image.split('.')[0]}")


	for i in range(0,len(x)):
		c2 = plt.Circle((x[i],y[i]),4, color = 'blue', linewidth=2,fill=False)
		ax.add_patch(c2)
		
	fontprops = fm.FontProperties(size=8)

	scalebar = AnchoredSizeBar(ax.transData,
                            75/4.7, '(75 arcsec)', 
                            'lower center', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)

	ax.add_artist(scalebar)
	
	plt.show()
	
	return True

def Run_TLE_Code(DayDirectory):
	
	files=os.listdir(f"{DayDirectory}/ProcessedImages")

	dic = {}
	for elem in files:
		if elem.split(".")[0][:-2] not in dic:
			dic[elem.split(".")[0][:-2]]=[]
		dic[elem.split(".")[0][:-2]].append(elem)
		
	result = list(dic.values())
	# result=[['PAER_SURV_4_263124_230302_122520+212638_89014Df_19.fit',
	# 'PAER_SURV_4_263124_230302_122520+212638_89014Df_20.fit',
	# 'PAER_SURV_4_263124_230302_122520+212638_89014Df_21.fit',
	# 'PAER_SURV_4_263124_230302_122520+212638_89014Df_22.fit']]
	
	for image_set in result:
		FIT_Files=image_set
		
		header_data=[]
		for img in range(0,len(FIT_Files)):
			Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + FIT_Files[img]}", unit = 'adu') 
			header=Data.header
			info=SPADEFunctions.ReadFilesHeader(DayDirectory,header,FIT_Files[img])
			header_data.append(info)
		
		for i in header_data:
			TLE=str(i['Year'])+str(i['Month'])+str(i['Day'])
		print(TLE)
		TLE_file=f'TLE Catalog-{TLE}.json'
		
		if os.path.exists(f'{TLE_Directory}/{TLE_file}'):
			print(FIT_Files)
			Satellite=SatellitePosition(config,DayDirectory,FileXY,FIT_Files,header_data,TLE_file,TLE_Directory)
			Write_TLE_Output_Json(FIT_Files,Satellite,TLE_Directory)
	
		else:
			print('TLE FILE missing')
			
		


Run_TLE_Code(DayDirectory)
