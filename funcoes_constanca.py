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

##################################################################################################################################
#This information is here so I can sucessfully run the following functions. 
#When copied to the SPADEFunctions.py file, there's no need to include them
#################################################################################################################################


PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
Base='MAPA_SURV_4_085283_220417_114004+060023_scan047_01'
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
config['LOGGER'] = int(config['LOGGER'])
all_files=['MAPA_SURV_4_085283_220417_114004+060023_scan047_01','MAPA_SURV_4_085283_220417_114004+060023_scan047_02',
'MAPA_SURV_4_085283_220417_114004+060023_scan047_03','MAPA_SURV_4_085283_220417_114004+060023_scan047_04']

##################################################################################################################################
# These functions create json files for Sextractor and Astrometry data and add Alt,Az information
##################################################################################################################################


def WriteSextractorFile_Json(DayDirectory,all_files):
	
	'''
	Function created to read catalogue created from sextractor in json type. It outputs 
	a list of dictionaries with the information related to each column as values of keys
	
	INPUTS:
	all_files			->   A list of every sextrator file in str type
	DayDirectory		->   Directory 
	
	FILES CREATED:
	SExtractorOutput.json
	'''
	
	print('[+] Extracting SExtractor Output from .cat file')
	
	
	SEx_directory = DayDirectory + '/SextractorCatalogue/'
	
	lst_jsonfile=[]
	
	for ficheiro in all_files:
		
		filename= SEx_directory + f"{ficheiro}.cat"

		A=fits.open(filename)
	
		for i in A[1].data:
			tmp_dict = {}
			tmp_dict['image']= str(ficheiro[-2])+str(ficheiro[-1])
			tmp_dict['x_image'] = str(i[6])
			tmp_dict['y_image'] = str(i[7])
			tmp_dict['WIN_SNR'] = str(i[12])
			lst_jsonfile.append(tmp_dict)
	
	with open('SExtractorOutput.json', 'w') as f:
		json.dump(lst_jsonfile, f, indent=3)

	print('[+] SExtractor output sucessfully extracted to SExtractorOutput.json file')

#WriteSextractorFile_Json(DayDirectory,all_files)

def WriteAstrometryOutput_Json(Base,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM):
	
	'''
	Creates a json Output File with information about detections from Astrometry:
	X,Y,(Ra,Dec) calculated, (Ra,Dec) catalogue, Mag calculated, Mag catalogue, FWHM, Altitude and Azimuthe. 
	The FWHM is the mean FWHM associated to each image (This information is in the image header)
	
	To run this function:
	X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM = SPADEFunctions.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)
	WriteAstrometryOutput_Json(Base,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)
	
	
	FILES CREATED:
	MAPA_SURV_4_085283_220417_114004+060023_scan047_01.json
	MAPA_SURV_4_085283_220417_114004+060023_scan047_02.json
	MAPA_SURV_4_085283_220417_114004+060023_scan047_03.json
	MAPA_SURV_4_085283_220417_114004+060023_scan047_04.json
	
	'''
	
	print('[+] Extracting Output data from Astrometry')
	
	items_2_store=[]
	
	#The location where the images were taken is necessary to convert Ra,Dec coordinates to Alt,Az
	location=EarthLocation(
			lat=40.182524*u.deg,
			lon=-7.873811*u.deg,
			height=820*u.m)
	
	
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


			if Base[-1]=="1":
				items['FWHM'] = "2.53222818375E+000"
				time = datetime.datetime(2022, 4, 17,23,47,53)  
			if Base[-1]=="2":
				items['FWHM'] = "1.83843030334E+000"
				time = datetime.datetime(2022, 4, 17,23,47,58)  
			if Base[-1]=="3":
				items['FWHM'] = "1.57364977002E+000"
				time = datetime.datetime(2022, 4, 17,23,48,3)  
			if Base[-1]=="4":
				time = datetime.datetime(2022, 4, 17,23,48,8)  
				items['FWHM'] = "1.63440303802E+000"
				
			
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
				
	print('[+] Json file created in OutputFiles folder')



#X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM = SPADEFunctions.CatalogComparison(Base,PipelineDirectory,DayDirectory,config)
#WriteAstrometryOutput_Json(Base,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)


##################################################################################################################################
# The following 2 functions, AngleDifferencePoints and AngleDifferenceLists are already in SPADEFunctions.py. They are only
# here so I can sucessfully run other functions
##################################################################################################################################


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


#################################################################################################################################
#        The purpose of the functions bellow is to find objects that are detected in more than  
#        one image using criteria that are explained in the function's comments, and
#        organize all detections in a list of dictionary lists.
#################################################################################################################################


def find_corr_in_list(pto, im_lst):
	
	corr_pto=None
	
	for im_pto in im_lst:
		
		L1=[float(pto["DecCalculated"]),float(pto['RaCalculated'])]
		L2=[float(im_pto[0]['DecCalculated']),float(im_pto[0]['RaCalculated'])]
		distance=AngleDifferenceLists(L1,L2)
			
		fwhm=[float(pto["FWHM"]),float(im_pto[0]['FWHM'])]
			
		if float(distance[0]) <= (4.7/3600)*max(fwhm) and float(distance[1]) <= (4.7/3600)*max(fwhm):

			corr_pto=im_pto[0]	
			continue

				
	return corr_pto


def find_correspondence(DayDirectory,all_files):
	
	'''	
	This function is a list of dictionaries with all detections from all images. If the distance between two point 
	(calculated by AngleDifferencePoints function) is less than a specific value, it means that the detections correspond
	to the same object in different images, so they're stored in the same list.
	The following example clarifies the ouput of the function.

	Example of input:
	Image 1=[{pt1A},{pt1B},{pt1C}]
	Image 2=[{pt2A},{pt2B}]
	Image 3=[{pt3A},{pt3C},{pt3E},{pt3F}]
	Image 4=[{pt4A},{pt4C},{pt4F}]

	Example of output:
	[[{pt1A},{pt2A},{pt3A},{pt4A}],[{pt1B},{pt2B}],[{pt1C},{pt3C},{pt4C}],[{pt3E}],[{pt3F},{pt4F}]]
	
	FILES CREATED:
	None

	'''
	print('[+] Start analysing data.')
	
	
	# Stores the data in 4 different lists, associated with each image
	
	im1=[]
	im2=[]
	im3=[]
	im4=[]
	
	for ficheiro in all_files:
		with open(f"{DayDirectory + '/OutputFiles/' + ficheiro}.json", 'r') as stream:
			jsonfile = json.load(stream)

			for i in jsonfile:
				if i['image'] == '01':
					im1.append([i])
				if i['image']  == '02':
					im2.append([i])
				if i['image']  == '03':
					im3.append([i])
				if i['image']  == '04':
					im4.append([i])
					
	print(f'Image 1 has {len(im1)} objects.')
	print(f'Image 2 has {len(im2)} objects.')
	print(f'Image 3 has {len(im3)} objects.')
	print(f'Image 4 has {len(im4)} objects.')
	
	
	soma=len(im1)

	
	# In the following part, the algorithm fixates the elements of the list associated with image 1 and tries to find matches
	# with the objects in the image 2,3 and 4. If a match is found, the list associated with that element will 
	# extend with the data from other images.
	
	# EXAMPLE:
	# Before comparing im1 with im2,im3,im4 --->  im1=[[{ptA1}],[{ptB1}]]
	# After comparing im1 with im2,im3,im4 --->  im1=[[{ptA1},{ptA2},{ptA4}],[{ptB1}]]
	#This means that:
	# - Point A in image 1 is found also in image 2 and 4. 
	# - Point B only exists in image 1
	
	im_a_comparar=[]
	im_a_comparar.append(im2)
	im_a_comparar.append(im3)
	im_a_comparar.append(im4)

	for c in im1: 
		
		for im in im_a_comparar: 
			
			ii = -1
			while c[ii] == {}:
				ii += -1
				
			pt_2_find_in_next_im= c[ii]
			new_pt= find_corr_in_list(pt_2_find_in_next_im,im) #This function determines if these two points in different images are associated to the same object

			if new_pt == None:
				continue
			else:
				c.append(new_pt)
	
	print(len(im1))	
	# # You now have an updated im1 list, with all elements associated with the objects from image 1 and their match in image 2,3 and 4.
	# # It's important to note that an object may not be present in image 1 and appear in image 2,3 and 4.
	# # Therefore, it's important to repit the same analysis, only now fixating the elements from image 2 that weren't matched
	# # with any element from image 1, and compare them with image 3 and 4.
	# # For that, the following steps where taken:
	
	# # The algoritm stores the indexs of the elements of image 2 that where matched with detections from image 1
	
	# idx=[]
	# for d in im2: 
		# for c in im1:
			# for i in range (0,len(c)):
				# if c[i]['image'] == '2' and c[i]['Dec_calculated'] == d[0]['Dec_calculated']: #significa que está
					# index = im2.index(d)
					# idx.append(index)
	
	
	# # With the indexes, the algoritm delets all elements from image 2 that were already matched with elements from image 1.

	# for index in sorted(idx, reverse=True): 
		# del im2[index]
	# soma+=len(im2)
	
	# print(f'Image 2 has {len(im2)} objects that do not have a match in image 1.')
	
	
	# # As before, the algoritm fixates the elements of the list associated with image 2 and tries to find any matches
	# # with the objects in the image 3 and 4. If a match is found, im1 list is updated.
	
	# # EXAMPLE:
	# # im1 before im2 being compared with im3 and im4  --->  im1=[[{ptA1},{ptA2},{ptA4}],[{ptB1}]]
	# # After comparing im2 with im3,im4 --->  im1=[[{ptA1},{ptA2},{ptA4}],[{ptB1}],[{ptC2},{ptC3},{ptC4}],[{ptD2},{ptD3}]]
	
	# #This means that:
	# # - Point A in image 1 is found also in image 2 and 4. 
	# # - Point B only exists in image 1
	# # - Point C exists in image 2,3 and 4
	# # - Point D exists in image 2 and 3
	
	
	# im_a_comparar=[]
	# im_a_comparar.append(im3)
	# im_a_comparar.append(im4)

	
	# for c in im2: 
		
		# for im in im_a_comparar: 
			
			# ii = -1
			# while c[ii] == {}:
				# ii += -1
				
			# pt_2_find_in_next_im= c[ii]
			# new_pt= find_corr_in_list(pt_2_find_in_next_im,im)
			
			# if new_pt == None:
				# continue
			# else:
				# c.append(new_pt)
			
	# im1.extend(im2)
	
	
	# # As before, the algoritm stores the indexes of the elements of image 3 that where matched with detections from im1 list
	# idx=[]
	# for d in im3: 
		# for c in im1:
			# for i in range (0,len(c)):
				# if c[i]['image'] == '3' and c[i]['Dec_calculated'] == d[0]['Dec_calculated']:
					# index = im3.index(d)
					# idx.append(index)
	
	# # With the indexes, the algoritm delets all elements from image 3 that were already matched with elements from im1 list
	
	# for index in sorted(idx, reverse=True): 
		# del im3[index]
	# soma+=len(im3)
	# print(f'Image 3 has {len(im3)} objects that do not have a match in image 1 and 2.')
	
	# im_a_comparar=[]
	# im_a_comparar.append(im4)
	
	# for c in im3: 
		
		# for im in im_a_comparar: 
			
			# ii = -1
			# while c[ii] == {}:
				# ii += -1
				
			# pt_2_find_in_next_im= c[ii]
			# new_pt= find_corr_in_list(pt_2_find_in_next_im,im)
			
			# if new_pt == None:
				# continue
			# else:
				# c.append(new_pt)
				
	# im1.extend(im3)

	
	# idx=[]
	# for d in im4: 
		# for c in im1:
			# for i in range (0,len(c)):
				# if c[i]['image'] == '4' and c[i]['Dec_calculated'] == d[0]['Dec_calculated']:
					# index = im4.index(d)
					# idx.append(index)

	# for index in sorted(idx, reverse=True): 
		# del im4[index]
	# soma+=len(im4)	
	# print(f'Image 4 has {len(im4)} objects that do not have a match in image 1, 2 and 3.')
	
	# im1.extend(im4)
	
	# print(f'The total number of detections and their corresponding matches is {len(im1)}')
	# print(f'This result is correct if {len(im1)} is equal to {soma}')
	
	
	return im1


find_correspondence(DayDirectory,all_files)
