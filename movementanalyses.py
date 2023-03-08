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
		if ficheiro[0] == 'P':			
			ficheiro=ficheiro[:-4]
			all_files.append(ficheiro)
		
	all_files.sort()
		
	return all_files

def number_of_detections(DayDirectory):
	
	all_files = files(DayDirectory)
	for ficheiro in all_files:
		with open(f"{DayDirectory + '/OutputFiles/' + ficheiro}.json", 'r') as stream:
				jsonfile = json.load(stream)
		match=0
		nonmatch=0
		for i in jsonfile:
			if i['Match']=='1':
				match+=1
			if i['Match']=='0':
				nonmatch+=1
		print(f'Image {ficheiro[-2]+ficheiro[-1]}: Match={match},NonMatch={nonmatch}')

number_of_detections(DayDirectory)
	
def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return round(mult*deg), round(mult*mnt), float("{:.1f}".format(mult*sec))
	

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
