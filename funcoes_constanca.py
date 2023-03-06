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
import excel2json


##################################################################################################################################
#This information is here so I can sucessfully run the following functions. 
#When copied to the SPADEFunctions.py file, there's no need to include them
##################################################################################################################################

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
#Base='MAPA_SURV_4_085283_220417_114004+060023_scan047_03'
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
config['LOGGER'] = int(config['LOGGER'])

def AngleDifferencePoints(Coord1,Coord2):
	
	Ra1  = np.radians(Coord1[0])
	Dec1 = np.radians(Coord1[1])
	Ra2  = np.radians(Coord2[0])
	Dec2 = np.radians(Coord2[1])
	

	difference = np.arccos(np.sin(Dec1)*np.sin(Dec2)+np.cos(Dec1)*np.cos(Dec2)*np.cos(Ra1-Ra2)) 
	
	return np.degrees(difference)
	
def AngleDifferenceLists(L1,L2):
	
	''' Calculates Angular difference between two sets of data Keys: 'RA' and 'DEC'
	L -> {RA:[Ra1,Ra2],DEC:[Dec1,Dec2]}
	'''
	
	Dists = []

	for i in range(len(L1)):
		Coord1 = [L1[i],L1[i]]
		Coord2 = [L2[i],L2[i]]
		Dists.append(AngleDifferencePoints(Coord1,Coord2))
	
	
	return np.array(Dists)

##################################################################################################################################
#The following information needs to be included in the main code
##################################################################################################################################

Links=['https://celestrak.org/NORAD/elements/gp.php?SPECIAL=gpz-plus&FORMAT=tle']


##################################################################################################################################
# 								 OUTPUT FILES, HEADER INFORMATION, CONFIG FILE INFORMATION
##################################################################################################################################

def files(DayDirectory):
	'''
	Function created to store a list of the Bases of all the analysed images
	
	INPUTS:
	DayDirectory 		-> directory of the present day									(str)
	
	RETURNS:
	all_files     		-> sortted list of all analysed files					        (lst)
	
	FILES CREATED:
	None
	'''
	
	all_files=[]
	for ficheiro in os.listdir(path=DayDirectory + '/ProcessedImages/'):				
		ficheiro=ficheiro[:-5]
		all_files.append(ficheiro)
		
	all_files.sort()
		
	return all_files
	
def WriteSextractorFile_Json(DayDirectory):
	'''
	Function created to read catalogue created from sextractor in json type. It outputs 
	a list of dictionaries with the information related to each column as values of keys
	
	INPUTS:
	DayDirectory		->   directory of the present day									(str)
	
	FILES CREATED:
	SExtractorOutput.json
	'''
	
	all_files=files(DayDirectory)
	
	SEx_directory = DayDirectory + '/SextractorCatalogue/'
	
	lst_jsonfile=[]
	
	for ficheiro in all_files:
		
		filename= SEx_directory + f"{ficheiro}.cat"

		A=fits.open(filename)
	
		for i in A[1].data:
			tmp_dict = {}
			tmp_dict['image']= str(ficheiro[-2])+str(ficheiro[-1])
			tmp_dict['number'] = str(i[0])
			tmp_dict['Flux_isocor'] = str(i[1])
			tmp_dict['Mag_isocor'] = str(i[2])
			tmp_dict['MagErr_isocor'] = str(i[3])
			tmp_dict['Flux_aper'] = str(i[4])
			tmp_dict['FluxErr_aper'] = str(i[5])
			tmp_dict['x_image'] = str(i[6])
			tmp_dict['y_image'] = str(i[7])
			tmp_dict['Flags'] = str(i[8])
			tmp_dict['FWHM_Image'] = str(i[9])
			tmp_dict['Elongation'] = str(i[10])
			tmp_dict['Flux_Radius'] = str(i[11])
			tmp_dict['SNR_WIN'] = str(i[12])
			lst_jsonfile.append(tmp_dict)
	
	with open('SExtractorOutput.json', 'w') as f:
		json.dump(lst_jsonfile, f, indent=3)

def str2bool(BoolStr):
	''' 
	Function created to transform string to boolean
	
	INPUTS:
	BoolStr 			-> Boolean in string											(str)
	
	OUTPUTS:
	BoolStr				-> Boolean in bool												(bool)
	
	FILES CREATED:
	NONE
	'''
	if BoolStr == "True":
		return True
	if BoolStr == "False":
		return False
	else:
		return BoolStr

def ReadPipelineConfig(ficheiro):
	''' 
	Function to read pipeline config files. It skips over blank lines and lines which start with an #. 
	Returns in the form of a dictionary with properties as keys and True/False or numbers as values
	
	INPUTS:
	ficheiro 			-> Config file directory										(str)
	
	RETURNS:
	dic					-> Dictionary with config keys									(dict)
	
	FILES CREATED:
	NONE
	'''
	
	
	dic = {}
	with open(ficheiro,'r') as texto:
		while True:
			linha = texto.readline()
			if len(linha)!= 0 and linha[0] != '#' and linha[0]!= '\n':
				B = linha.strip().split()
				if B[1] == 'True' or 'False':
					dic[B[0]] = str2bool(B[1])
				else:
					dic[B[0]] = int(B[1])
			if linha == "":
				break
	return dic

def ReadFilesHeader(DayDirectory):
	
	'''
	Function created to return a list of dicts with the header information (time,date,scale,mean fwhm,dimensions) of all processed 
	files. The size of the list is associated with the number os processed files
	
	INPUTS:
	DayDirectory	->   directory of the present day									(str)
	
	RETURNS:
	header	->   header information of all processed images                  		(list of dicts)
	
	FILES CREATED:
	None
	'''
	
	header=[]
	all_files=files(DayDirectory)

	for ficheiro in all_files:
		
		Header_information={}
		Data = CCDData.read(f"{DayDirectory + '/ProcessedImages/' + ficheiro}.fits", unit = 'adu') 
		
		Header_information['Image'] =ficheiro[-2]+ficheiro[-1] #Stores the last two strings of file name, that is associated to the image number, so you know to which image the header information is associated to.
		Header_information['FWHM'] =Data.header['FWHM']
		Header_information['Dimensions'] =Data.header['NAXIS1']
		Header_information['Scale']=float(Data.header['SCALE'])
		Header_information['Alt']=float(Data.header['ALTUR'])
		Header_information['Az']=float(Data.header['AZIMUTH'])
		
		#In current file headers, date and time information is stored in the format: 21/04/23 and 13:09:12.2 (an example)
		#The code bellow splits the caracters from : and / and returns the strings that begin with 0 with only the second 
		#element (04 --> 4, 09 --> 9).This is helpful when using datetime module in further analysis, that can't interpret 
		#month 04 or hour 09. It also adds the string 20 to the caracter asscociated to the year for the same reason (23 --> 2023).
		
		Day=Data.header['Date'].split("/") 
		for i in Day: 
			if i[0]=='0':
				idx=Day.index(i)
				Day[idx]=i[1]

		Time_UT=Data.header['FULLTIME'].split(":")

		for i in Time_UT:
			if i[0]=='0':
				idx=Time_UT.index(i)
				Time_UT[idx]=i[1]
	
		Header_information['Year']=int('20'+Day[2])
		Header_information['Month']=int(Day[1])
		Header_information['Day']=int(Day[0])
		Header_information['Hour']=int(Time_UT[0])
		
		
		# This part of the algoritm adds half of the exposure time to the time when the image was taken.
		exp_time=float(Data.header['EXPTIME'])
		minute=int(Time_UT[1])
		second=float(Time_UT[2])+exp_time/2
		
		if second >= 60:
			Header_information['Minute']=minute+1
			Header_information['Second']=int(second-60) 
			Header_information['Microsecond']=int(abs(int(second)-second)*10**6)
			header.append(Header_information)
		else:
			Header_information['Minute']=minute
			Header_information['Second']=int(second)
			Header_information['Microsecond']=int(abs(int(second)-second)*10**6)
			header.append(Header_information)

	return header

def WriteAstrometryOutput_Json(Base,DayDirectory,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM):
	
	'''
	Function that creates a json file with the information about detections from Astrometry:
	X,Y,(Ra,Dec) calculated, (Ra,Dec) catalogue, Mag calculated, Mag catalogue, Altitude and Azimuthe. 
	
	INPUTS:
	Base 				-> File name without path and ending 						(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	Detections			-> SkyCoord Object with image detections (Ra/Dec)			(SkyCoord)
	DetectionsMags		-> Instrumental magnitudes corresponding to Detections		(np.array)
	Catalogue			-> SkyCoord Object with catalogue matches (Ra/Dec)			(SkyCoord)
	CatalogueMags		-> G magnitudes corresponding to Detections					(np.array)
	found				-> Boolean array with matches								(np.array)
	
	OUTPUTS:
	NONE
	
	FILES CREATED:
	{Base}.json --> the number of files created is related to the number of bases
	
	'''

	items_2_store=[]
	config_information = ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
	
	#The location where the images were taken is necessary to convert RA,DEC coordinates to Alt,Az coordinates
	location=EarthLocation(
			lat=float(config_information['Latitude'])*u.deg,
			lon=float(config_information['Longitude'])*u.deg,
			height=float(config_information['Height'])*u.m)
	
	#Header Information of proessed files
	header=ReadFilesHeader(DayDirectory)
	
	for c in range(0,len(header)):
		
		#The algoritm only continues when the image of which the header information is associated to is the same as the image of current Base
		if header[c]['Image'] == str(Base[-2])+str(Base[-1]):
			time = datetime.datetime(header[c]['Year'],header[c]['Month'],header[c]['Day'],header[c]['Hour'],header[c]['Minute'],header[c]['Second'],header[c]['Microsecond'])
	
			for i in range(len(Detections)):
			
					Match = found[i]
					items={}
					items['image']=str(Base[-2])+str(Base[-1])
					items['Xdetect']=str(X[i])
					items['Ydetect']= str(Y[i])
					items['RaCalculated'] = str(Detections.ra[i].value)
					items['DecCalculated'] = str( Detections.dec[i].value)
				
					if Match:
						items['Match']= '1'
						items['RaCatalogue'] = str(Catalogue.ra[i].value)
						items['DecCatalogue'] = str(Catalogue.dec[i].value)
						items['MagInst '] = str(DetectionsMags[i])
						items['MagReal'] =str(CatalogueMags[i])
						
					if Match == False:
						items['Match']= '0'
						items['MagInst '] = str(DetectionsMags[i])
						items['MagReal'] ='0'
					
		
					dec=(Detections.dec[i].value)*u.deg
					ra=(Detections.ra[i].value)*u.deg
					stars = SkyCoord(ra, dec, frame='icrs')
					
					#Converts Ra,Dec to Alt,Az
					stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))
					sky_coords = [stars_altaz.alt.degree,stars_altaz.az.degree]
					
					items['Alt']=str(sky_coords[0])
					items['Az']=str(sky_coords[1])
		
					items_2_store.append(items)
					
					with open(f"{DayDirectory + '/OutputFiles/' + Base}.json", 'w') as f: 
						json.dump(items_2_store, f, indent=3)
				
#X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM = SPADEFunctions.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)
#WriteAstrometryOutput_Json(Base,DayDirectory,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)

#################################################################################################################################
#        The purpose of the functions bellow is to find objects that are detected in more than  
#        one image using criteria that are explained in the function's comments, and
#        organize all detections in a list of dictionary lists. 
#################################################################################################################################


def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return round(mult*deg), round(mult*mnt), float("{:.1f}".format(mult*sec))
	
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
	
def criteria(pto, im_lst,fix_img,idx,header,scale):
	
	corr_pto=None
	
	FWHM_image=float(header[fix_img]['FWHM'])
	FWHM_next_image=float(header[idx]['FWHM'])
	
	for im_pto in im_lst:

		L1=[float(pto[0]['RaCalculated']),float(pto[0]['DecCalculated'])]
		L2=[float(im_pto['RaCalculated']),float(im_pto['DecCalculated'])]
		distance=AngleDifferencePoints(L1,L2)
	
		fwhm=[FWHM_image,FWHM_next_image]

		if float(distance) <= scale*max(fwhm)/3600:
			corr_pto=im_pto	
			continue

	return corr_pto

def corres_images(fix_list,image,fix_img,next_image,header,scale):  
	
	
	'''
	Function created to fixate a point in fix_list and compare it to all other points in image.
	
	INPUTS:
	fix_list	->          directory of the present day									(str)
	
	RETURNS:
	header	->   header information of all processed images                  		(list of dicts)
	
	FILES CREATED:
	None
	'''
	
	for c in fix_list:  #fixates a point an image to compares to all other points from different images
		  for i in range(next_image,len(image)): 
			  for im in image[i]: #goes through all points of an image
				  pt_2_find_in_next_im = c
				  new_pt= criteria(pt_2_find_in_next_im,im,fix_img,i,header,scale) 
				  if new_pt == None:
					  continue
				  else:
					  c.append(new_pt)
	return fix_list

def find_correspondences(DayDirectory):
	
	'''	
	This function is a list of dictionaries with all detections from all images. If the distance between two point 
	(calculated by AngleDifferencePoints function) is less than a specific value, it means that the detections correspond
	to the same object in different images, so they're stored in the same list.
	It's important to notice that the algorithm is currently treating every object has a potencial space debris detection,
	i.e the stars are not being distinguished from the rest.
	The following example clarifies the ouput of the function.

	Example of input:
	
	Obj=[[Image1],[Image2],[Image3],...,[ImageN]]
	
	Image 1={pt1A},{pt1B},{pt1C}
	Image 2={pt2A},{pt2B}
	Image 3={pt3A},{pt3C},{pt3E},{pt3F}
	...
	Image N={ptNA},{ptNC},{ptNF}

	Example of output:
	[[{pt1A},{pt2A},{pt3A},{ptNA}],[{pt1B},{pt2B}],[{pt1C},{pt3C},{ptNC}],[{pt3E}],[{pt3F},{ptNF}]]
	
	FILES CREATED:
	None
	'''
	all_files=files(DayDirectory)
	
	#Creates json Output files for astrometry data of all processed images. If such files are already created, the next three
	#lines can be commented
	# for ficheiro in all_files:
		# X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM = SPADEFunctions.CatalogComparison(ficheiro,PipelineDirectory,DayDirectory,config)
		# WriteAstrometryOutput_Json(ficheiro,DayDirectory,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)
	
	#Creates a list of N lists, (N = total number of images)
	obj = [[] for i in range(0,len(all_files))]
	print(f'[+] Start analysing data of {len(obj)} images.')
	
	#Appends data from image 1 to list of index [0], image 2 to index [1],...,image N to index[N-1]
	count=0
	for ficheiro in all_files:
		with open(f"{DayDirectory + '/OutputFiles/' + ficheiro}.json", 'r') as stream:
			jsonfile = json.load(stream)
			for i in jsonfile:
				obj[count].append([i])
			count+=1
			
	#Informs the number of detections in each image
	for i in range(0,len(all_files)):
		print(f'Image{all_files[i][-1]} has {len(obj[i])} objects.')

	header=ReadFilesHeader(DayDirectory)
	
	fix_img=0
	next_image=1
	
	#In this part, the algoritm creates a correspondences list, as a result of corres_images function.
	correspondences=corres_images(obj[fix_img],obj,fix_img,next_image,header,header[fix_img]['Scale'])
	
	while next_image<len(obj):		
		idx=[]
		img=obj[next_image]
		for d in img: 
			for c in correspondences:
				for i in range (0,len(c)):
					if c[i]['image'] == d[0]['image'] and c[i]['DecCalculated'] == d[0]['DecCalculated']: #significa que está
						index = img.index(d)
						idx.append(index)

		for index in sorted(idx, reverse=True): 
			del img[index]
		

		if next_image == len(obj)-1:
			correspondences.extend(img)
			next_image+=1
		else:
			fix_img+=1
			next_image+=1
			cors=corres_images(img,obj,fix_img,next_image,header,header[fix_img]['Scale'])
			correspondences.extend(cors)
			
	#scatteredplot_all_detections(correspondences)	
	return correspondences

#find_correspondences(DayDirectory)

################################################################################################################################
#														TLES   (ainda em construção)
################################################################################################################################

def ReadTLES_json(Links):
	
	'''
	ReadTLES(Links) function uploads tle information from a url and saves it in a .json file (a list of dicts).
	
	INPUTS:
	Links		->   Celestrak links of Satellites TLEs                               (lst)
	
	FILES CREATED:
	GEO_tles.json
	'''
	
	print('[+] Downloading TLE data...')
	tle_jsonfile=[]
	
	for i in range(len(Links)):
		response = requests.get(Links[i])

	if response.status_code == 200: #The request was successful and the server responded with the data you were requesting
		tmp_dict = {}
		tle_data=response.text.split('\n')
		count=0
		
		for i in range(0,len(tle_data)-1):
			
			if tle_data[i][0] == '1':
				tmp_dict['tle_1'] = tle_data[i].strip()
			if tle_data[i][0] == '2':
				tmp_dict['tle_2'] = tle_data[i].strip()
			if tle_data[i][0] != '1' and tle_data[i][0] != '2':
				tmp_dict['satellite_name'] = tle_data[i].strip()

			if "tle_1" in tmp_dict and "tle_2" in tmp_dict and "satellite_name" in tmp_dict:
				tle_jsonfile.append(tmp_dict)
				tmp_dict = {}
			else:
				pass
				
		with open('GEO_tles.json', 'w') as f:
			json.dump(tle_jsonfile, f, indent=3)
			print('[+] TLE data downloaded in GEO_tles.json')
	else:
		print('[+] ERROR in downloading TLE data')

def TLE_plot(Time,Alt,Az,time_ticks,nome):
	
	color_alt = "#FF0000"
	color_az = "#00FF00"

	fig, ax1 = plt.subplots()
	ax2 = ax1.twinx()
	
	ax1.plot(Time, Alt,color=color_alt)
	ax2.plot(Time, Az,color=color_az)
	
	ax1.set_xlabel("Time")
	ax1.set_ylabel("Elevantion (degrees)",color=color_alt, marker='o',fontsize=14)
	ax1.tick_params(axis="y", labelcolor=color_alt)
	
	ax2.set_ylabel("Azimuth (degrees)",color=color_az,marker='o', fontsize=14)
	ax2.tick_params(axis="y", labelcolor=color_az)
	
	fig.suptitle(f"Satellite {nome}", fontsize=20)
	plt.xticks(Time, time_ticks,rotation=0,fontsize=8)
	plt.show()
						
def SatellitePosition(DayDirectory,PipelineDirectory,config_file):
	'''
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	with open("Worksheet.json", 'r') as stream:
		lista_ficheiro=json.load(stream)
   
	
	header=ReadFilesHeader(DayDirectory)
	
	Alt_target=header[0]['Alt']
	Az_target=header[0]['Az']
	scale=header[0]['Scale']
	dimensions=header[0]['Dimensions']
	
		
	FOV=(scale*dimensions)/3600
	
	Alt_window=[Alt_target-FOV/2,Alt_target+FOV/2]
	Az_window=[Az_target-FOV/2,Az_target+FOV/2]
	
	 #The Timescale object returned by load.timescale() manages the conversions between different 
	 #time scales and is also how the programmer builds Time objects for specific dates. 
	 
	ts = load.timescale()

	location = ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
						
	for i in lista_ficheiro:
		items={}
		line1=i['tle_1']
		line2=i['tle_2']
		nome=i['satellite_name']
	
		satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 

		bluffton = wgs84.latlon(float(location['Latitude']), float(location['Longitude']),float(location['Height']))
		
		difference = satellite - bluffton 
		
		second_firstimage=float(str(header[0]['Second'])+'.'+(str(header[0]['Microsecond'])))
		second_lastimage=float(str(header[-1]['Second'])+'.'+(str(header[-1]['Microsecond'])))
		
		t0 = ts.utc(header[0]['Year'], header[0]['Month'], header[0]['Day'], header[0]['Hour'],header[0]['Minute'],second_firstimage) 
		
		t1 = ts.utc(header[-1]['Year'], header[-1]['Month'], header[-1]['Day'], header[-1]['Hour'],header[-1]['Minute'],second_lastimage)
		
		Alt=[]
		Az=[]
		Time=[]
		time_ticks=[]

		for t in np.arange(t0,t1,0.00001):

		
			topocentric = difference.at(t)
			
			alt, az, distance = topocentric.altaz()

			if alt.degrees > Alt_window[0] and alt.degrees< Alt_window[1]:
				if az.degrees > Az_window[0] and az.degrees< Az_window[1]:
					print(f"Satelite: {nome}, Alt: {alt.degrees}, Az: {az.degrees}")
					
SatellitePosition(DayDirectory,PipelineDirectory,config_file)		
					
					# Alt.append(alt.degrees)
					# Az.append(az.degrees)
					# Time.append(t.tdb)
					# strg=f"{t.utc.hour}:{t.utc.minute}:{round(t.utc.second)}"
					# time_ticks.append(strg)

					

	
					
		# if len(Alt)!= 0:
			# TLE_plot(Time,Alt,Az,time_ticks,nome)
			
		






def convertxlsx_to_json():
	xlsx_file='Oscar - Satellites 2023-03-02.xlsx'
	excel2json.convert_from_file(xlsx_file)

def ReadOscarJsonData():
	
	with open('Worksheet.json', 'r') as stream:
		lista_ficheiro=json.load(stream)
		
	print(f" There are {len(lista_ficheiro)} satellites")
	
	inactive_sat=[]
	for i in lista_ficheiro:
		if i['Sat status'] == "Inactive":
			inactive_sat.append(i["Acronym"])
	print(f"{ len(inactive_sat)} of those are inactive")
	print(inactive_sat)

def SatellitePosition(PipelineDirectory,config_file):
	'''
	SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
	and calculates their position in a given time and location.
	'''
	
	with open("GEO_tles.json", 'r') as stream:
		lista_ficheiro=json.load(stream)
		
	ts = load.timescale()
	
	config_information = ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
	
	location=EarthLocation(
			lat=float(config_information['Latitude'])*u.deg,
			lon=float(config_information['Longitude'])*u.deg,
			height=float(config_information['Height'])*u.m)
			
			
	sat=[]

	for i in lista_ficheiro:
		items={}
		line1=i['tle_1']
		line2=i['tle_2']
		nome=i['satellite_name']
	
		satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 
		sat.append(satellite)
	step = 1
	with open('data.txt', 'w') as f:
		for chunk in [sat[n:n+step] for n in range(0, len(sat), step)]:
			f.write(','.join(map(str, chunk)) + '\n')

#SatellitePosition(PipelineDirectory,config_file)



	
	







################################################################################################################################
#								RASCUNHOS : NÃO COLOCAR NO SPADEFunctions.py AINDA
################################################################################################################################


def separatedetections(lst1,lst2):	

	len1_str=[]
	len2_str=[]
	len3_str=[]
	len3_obj=[]
	
	for i in lst1:
		if len(i)==1:
			len1_str.append(i)
		if len(i)==2:
			len2_str.append(i)
		if len(i)==3:
			len3_str.append(i)

	for i in lst2:
		if len(i)==3:
			len3_obj.append(i)

	return len3_obj,len1_str,len2_str,len3_str




			
def scatterplotonedetection():
	
	detections=separatedetections(lst_cors_stars,lst_cors_nonidnt_objs)	
	fontprops = fm.FontProperties(size=8)
	fig, ax = plt.subplots()
	
	x=[]
	y=[]
	my_xticks=[]
	my_yticks=[]
	for i in detections[1]:
		for j in i:
			if j['image']=='1':

				plt.scatter(
					x=float(j['Ra_calculated']),
					y=float(j['Dec_calculated']),
					c='k',
					marker='o',
					alpha=0.25)
				plt.annotate('1',(float(j['Ra_calculated']),float(j['Dec_calculated'])))
			
			if j['image']=='2':

				plt.scatter(
					x=float(j['Ra_calculated']),
					y=float(j['Dec_calculated']),
					c='k',
					marker='o',
					alpha=0.5,
				)
				plt.annotate('2',(float(j['Ra_calculated']),float(j['Dec_calculated'])))
			
			if j['image']=='3':

				plt.scatter(
					x=float(j['Ra_calculated']),
					y=float(j['Dec_calculated']),
					c='k',
					marker='o',
					alpha=0.75,
				)
				plt.annotate('3',(float(j['Ra_calculated']),float(j['Dec_calculated'])))
			
			if j['image']=='4':

				plt.scatter(
					x=float(j['Ra_calculated']),
					y=float(j['Dec_calculated']),
					c='k',
					marker='o',
					alpha=1,
				)
				plt.annotate('4',(float(j['Ra_calculated']),float(j['Dec_calculated'])))
				
			x.append(float(j['Ra_calculated']))
		
			y.append(float(j['Dec_calculated']))
	
			x.sort()
			for i in x:
				my_xticks.append(decdeg2dms(i))


			plt.xticks([x[0], x[-1]], [my_xticks[0], my_xticks[-1]],rotation=0,fontsize=8)
	
			y.sort()
			for i in y:
				my_yticks.append(decdeg2dms(i))
			plt.yticks([y[0], y[-1]], [my_yticks[0], my_yticks[-1]],rotation=0,fontsize=8)
			
			
			
	
			scalebar = AnchoredSizeBar(ax.transData,
                            75/3600, '(75 arcsec)', 'lower left', 
                           pad=0.1,
                           color='black',
                           frameon=False,
                           #size_vertical=1,
                           fontproperties=fontprops)

			ax.add_artist(scalebar)
			plt.title("Objects detected in just 1 image")
			plt.ylabel('Dec (dms)')
			plt.xlabel('Ra (dms)')

	plt.show()
	
	
#scatteredplot_colors()
#scatterplotonedetection()


################################################################################################################################
#														Movement plots
################################################################################################################################

def plot_detection_movement(detections):
	
	x=[]
	y=[]
	x_dms=[]
	y_dms=[]
	
	image=0
	
	for i in detections:
		x.clear()
		y.clear()
		x_dms.clear()
		y_dms.clear()
		
		for n in i:
			
			point_x=(n['alt'])
			x.append(point_x)
			
			point_x_dms=(decdeg2dms(point_x))
			
			x_dms.append(str(point_x_dms))
			
		
			
			point_y=(n['az'])
			y.append(point_y)
			
			point_y_dms=(decdeg2dms(point_y))
			y_dms.append(str(point_y_dms))
			
		
		image+=1
		text = ["1", "2", "3", "4"]
		plt.plot(x, y,color='green', linestyle='dashed', linewidth = 3, marker='o', markerfacecolor='blue', markersize=4)
		for i in range(len(x)):
			plt.annotate(text[i], (x[i], y[i]))
		
		my_xticks = x_dms
		plt.xticks(x, my_xticks,rotation=45,fontsize=8)
		my_yticks = y_dms
		plt.yticks(y, my_yticks,rotation=0,fontsize=8)
		
	
		a=min(x)+(max(x)-min(x))/2
		b=min(y)+(max(y)-min(y))/2
		
		numx=(2*60)/(3600*2) #já estou a fazer metade aqui
		numy=(6*60)/(3600*2)  #já estou a fazer metade aqui 
		#num=(0.015*60)/(3600*2)
		plt.xlim(a-numx,a+numx)
		plt.ylim(b-numy,b+numy)

		
		plt.title(f'Star {image}')
		plt.ylabel(f'Az (degrees,min,sec)')
		plt.xlabel(f'Alt (degrees,min,sec)')
		plt.show()


################################################################################################################################
#														Image plots
################################################################################################################################

def plot_detections_in_image(Base,DayDirectory):
	
	'''
	Function created to plot the image being analysed and draw the detections on top of it, according to json output file
	
	INPUTS:
	all_files			->   A list of every sextrator file in str type
	DayDirectory		->   Directory 
	
	FILES CREATED:
	NONE
	'''
	
	
	with open(f"{DayDirectory + '/OutputFiles/' + Base}.json", 'r') as stream:
		jsonfile = json.load(stream)
	
	
	fig = plt.figure()
	ax = fig.add_subplot(111) 


	Imagem = CCDData.read(f"{DayDirectory  + '/ProcessedImages/' + Base}.fits", unit = 'adu') 
	plt.imshow(Imagem,cmap='gray')
	plt.title(f'Imagem{Base[-1]}')
	

	for i in jsonfile:
		c1 = plt.Circle((float(i['Xdetect']),float(i['Ydetect'])),4, color = 'red', linewidth=2,fill=False)
		ax.add_patch(c1)
		
	plt.show()
	
	return True

#plot_detections_in_image(Base,DayDirectory)

