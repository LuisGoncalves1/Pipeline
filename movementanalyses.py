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

##################################################################################################################################
#Number of matches and non matches in each image
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

#number_of_detections(DayDirectory)

##################################################################################################################################
#Associate detections from all analysed images

def AngleDifferencePoints(Coord1,Coord2):
	
	Ra1  = np.radians(Coord1[0])
	Dec1 = np.radians(Coord1[1])
	Ra2  = np.radians(Coord2[0])
	Dec2 = np.radians(Coord2[1])
	print(np.sin(Dec1)*np.sin(Dec2)+np.cos(Dec1)*np.cos(Dec2)*np.cos(Ra1-Ra2))
	

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
	
def criteria(pto, im_pto):
	
	corr_pto=None

	L1=[float(pto[0]['RaCalculated']),float(pto[0]['DecCalculated'])]
	L2=[float(im_pto['RaCalculated']),float(im_pto['DecCalculated'])]
	
	distance=AngleDifferencePoints(L1,L2)
	
		
	if float(distance) <= 0.331:
		corr_pto=im_pto	
	
	return corr_pto

def corres_images(fix_list,image):  
	
	for c in fix_list:  
		for im in image:
			pt_2_find_in_next_im = c
			new_pt= criteria(pt_2_find_in_next_im,im) #criteria for correspondence between to points
			if new_pt == None: 
				fix_list.append([im])
				continue
			else:
				c.append(new_pt)
	return fix_list

def find_correspondences(DayDirectory):

	all_files=files(DayDirectory)
	corrs=[]

	
	with open(f"{DayDirectory}/OutputFiles/{all_files[0]}.json", 'r') as stream:
		jsonfile = json.load(stream)
		for i in jsonfile:
			if i['Match']=='0':
				corrs.append([i])
	
	del all_files[0]
	for ficheiro in all_files:
		obj=[]
		with open(f"{DayDirectory}/OutputFiles/{ficheiro}.json", 'r') as stream:
			jsonfile = json.load(stream)
			for i in jsonfile:
				if i['Match']=='0':
					obj.append(i)
					
		correspondences=corres_images(corrs,obj)
		
		
			
	return correspondences



find_correspondences(DayDirectory)



# def find_correspondences(DayDirectory):
	
	# '''	
	# This function is a list of dictionaries with all detections from all images. If the distance between two point 
	# (calculated by AngleDifferencePoints function) is less than a specific value, it means that the detections correspond
	# to the same object in different images, so they're stored in the same list.
	# It's important to notice that the algorithm is currently treating every object has a potencial space debris detection,
	# i.e the stars are not being distinguished from the rest.
	# The following example clarifies the ouput of the function.

	# Example of input:
	
	# Obj=[[Image1],[Image2],[Image3],...,[ImageN]]
	
	# Image 1={pt1A},{pt1B},{pt1C}
	# Image 2={pt2A},{pt2B}
	# Image 3={pt3A},{pt3C},{pt3E},{pt3F}
	# ...
	# Image N={ptNA},{ptNC},{ptNF}

	# Example of output:
	# [[{pt1A},{pt2A},{pt3A},{ptNA}],[{pt1B},{pt2B}],[{pt1C},{pt3C},{ptNC}],[{pt3E}],[{pt3F},{ptNF}]]
	
	# FILES CREATED:
	# None
	# '''
	# all_files=files(DayDirectory)
	
	# #Creates a list of N lists, (N = total number of images)
	# obj = [[] for i in range(0,len(all_files))]
	# header = [[] for i in range(0,len(all_files))]
	
	# #Appends data from image 1 to list of index [0], image 2 to index [1],...,image N to index[N-1] of obj list
	# count=0
	# for ficheiro in all_files:
		# with open(f"{DayDirectory}/OutputFiles/{ficheiro}.json", 'r') as stream:
			# jsonfile = json.load(stream)
			# for i in jsonfile:
				# if i['Match']=='0':
					# obj[count].append([i])
			# count+=1

	# #Appends header data from image 1 to list of index [0], image 2 to index [1],...,image N to index[N-1] of header list
	# count=0
	# for ficheiro in all_files:
		# with open(f"{DayDirectory + '/OutputFiles/' + ficheiro[-2]+ficheiro[-1]}_header.json", 'r') as stream:
			# jsonfile= json.load(stream)
			# for i in jsonfile:
				# header[count].append(i)
			# count+=1	

	# fix_img=0
	# next_image=1
	
	# #In this part, the algoritm creates a correspondences list, as a result of corres_images function. It starts by
	# # fixating the first image and compare it with the rest
	# correspondences=corres_images(obj[fix_img],obj,fix_img,next_image,header,float(header[fix_img][0]['Scale']))
	
	# #So far, you have a correspondences list with all the detections from all images that associate with a point in image 1
	# #Ex: corresponences=[[{pt1A},{pt2A},{pt3A},{ptNA}],[{pt1B},{pt2B}],[{pt1C},{pt3C},{ptNC}],[{pt1D}]]
	# #The next step is to fixate the next image (image 2 for example) and compare it with image 3,4,...N, to find correspondences 
	# #that weren't made with the data from image 1.
	
	# while next_image<len(obj):	
		# idx=[]
		# img=obj[next_image]
		
		# #The algoritm looks for data in the next image that isn't stored the the correspondences list.
		# #If a point in the next image is in the correspondence list, it means that that point already has a correspondence. So,
		# #the algorithm stores the index of that detection and delets it from img list. Therefore, when doing the 
		# #corres_images again, the algorithm only analyses data that hasn't been associated with previous images.
		# for detection in img: 
			# for associations in correspondences:
				# for i in range (0,len(associations)):
					# if associations[i]['image'] == detection[0]['image'] and associations[i]['DecCalculated'] == detection[0]['DecCalculated']: 
						# index = img.index(detection) 
						# idx.append(index)

		# for index in sorted(idx, reverse=True): 
			# del img[index]
		
		# #if the next image counter is the same size as the size of obj (a list of lists of all data), it means that
		# #the algoritm has already reached the last image, so the points that weren't deleted from img list are points that don't
		# #have a correspondence. Therefore, there are no more images to compare those points to, so the algorithm simply extends
		# #the correspondence list for those points
		# if next_image == len(obj)-1:
			# correspondences.extend(img)
			# next_image+=1
		# #If that's not the case, the algoritm find new correspondences with corres_images, only now fixating the next image
		# #and compare it with the following.
		# else:
			# fix_img+=1
			# next_image+=1
			# cors=corres_images(img,obj,fix_img,next_image,header,header[fix_img][0]['Scale'])
			# correspondences.extend(cors)
				
	# return correspondences


# def criteria(pto, im_lst,fix_img,idx,header,scale):
	
	# corr_pto=None
	# FWHM_image=float(header[fix_img][0]['FWHM'])
	# FWHM_next_image=float(header[idx][0]['FWHM'])
	# for im_pto in im_lst:

		# L1=[float(pto[0]['RaCalculated']),float(pto[0]['DecCalculated'])]
		# L2=[float(im_pto['RaCalculated']),float(im_pto['DecCalculated'])]
		# distance=AngleDifferencePoints(L1,L2)
	
		# fwhm=[FWHM_image,FWHM_next_image]

		# #if float(distance) <= config['max_sep']scale*max(fwhm)/3600:
		# #	corr_pto=im_pto	
		# #	continue

	# return corr_pto

# def corres_images(fix_list,image,fix_img,next_image,header,scale):  
	
	
	# '''
	# Function created to fixate a point in fix_list and compare it to all other points in image.
	
	# INPUTS:
	# fix_list	->          directory of the present day									(str)
	
	# RETURNS:
	# header	->   header information of all processed images                  		(list of dicts)
	
	# FILES CREATED:
	# None
	# '''
	# for c in fix_list:  #fixates a point an image to compares to all other points from different images
		  # for i in range(next_image,len(image)): 
			  # for im in image[i]: #goes through all points of an image
				  # pt_2_find_in_next_im = c
				  # new_pt= criteria(pt_2_find_in_next_im,im,fix_img,i,header,scale) #criteria for correspondence between to points
				  # if new_pt == None: 
					  # continue
				  # else:
					  # c.append(new_pt)
	# return fix_list	
	


#######################################################################################################################################
#Plots

def represent_correspondences(DayDirectory):
		
	lista=find_correspondences(DayDirectory)
	counts=[]
	for i in lista:
		counts.append(len(i))
	print(max(counts))
	print(min(counts))
	# plt.hist(counts)
	# plt.xlabel('Number of correspondences between detections')
	# plt.ylabel('counts')
	# plt.xlim([0,26])
		

def plot_detection_movement(DayDirectory):
	
	
	
	x=[]
	y=[]
	image=[]
	x_dms=[]
	y_dms=[]
	
	lista=find_correspondences(DayDirectory)
	detection=1
	
	for i in detections:
		x.clear()
		y.clear()
		x_dms.clear()
		y_dms.clear()
		image.clear()
		
		for n in i:
			
			image_number=(n['Image'])
			image.append(image_number)
			
			point_x=(n['Ra_calculated'])
			x.append(point_x)
			
			point_x_dms=(decdeg2dms(point_x))
			x_dms.append(str(point_x_dms))
			
		
			point_y=(n['Dec_calculated'])
			y.append(point_y)
			
			point_y_dms=(decdeg2dms(point_y))
			y_dms.append(str(point_y_dms))
			
		
	
		plt.plot(x, y,color='green', linestyle='dashed', linewidth = 3, marker='o', markerfacecolor='blue', markersize=4)
		for i in range(len(x)):
			plt.annotate(image[i], (x[i], y[i]))
		
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

		
		plt.title(f'Object {detection}')
		plt.ylabel(f'Az (degrees,min,sec)')
		plt.xlabel(f'Alt (degrees,min,sec)')
		detection+=1
		plt.show()
