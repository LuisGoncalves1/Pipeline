from Funcoes import FuncoesReducao
import os
import numpy as np
import math
import yaml
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import EarthLocation
import datetime
import astropy.units as u

###########################################################################################################################
### This function reads the output file, and stores the Ra Dec coordinates of the detections 
### that weren't identified as stars by the catalogue in a .yml file (list of dicts).
### This function works for now, but right now it only works for the 4 images of Defesa. In time I will improve
### the function and it will be automatized for whatever set of images


def convert_output_file():
	
	ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_0'
	directory = os.getcwd()
	data = "2022-03-20"
	diretoria = "%s/Files/Images/%s" %(directory,str(data))

	count=1
	while count<5:

		ficheiro= ficheiro+ str(count)
		filepath=diretoria + '/OutputFiles/' + ficheiro

		with open(filepath,'r') as fp:
			all_the_lines= fp.readlines()
			items=[]
			for i in all_the_lines:
				items.append(i)

		del(items[0])

		items_2_store = []
	
		for i in items:
			resultados=i.split(',')
			if float(resultados[2]) == 0:
				tmp = {}
				#tmp['x'] = float(resultados[0])
				#tmp['y'] = float(resultados[1])
				tmp['Ra'] = float(resultados[3])
				tmp['Dec'] = float(resultados[4])
				items_2_store.append(tmp)
	
		with open(f"data_image_{count}.yml", 'w') as outfile:
			yaml.dump(items_2_store, outfile, default_flow_style=False)
					
		ficheiro=ficheiro[:-1]
		count+=1
		
output=convert_output_file()

#########################################################################################################################################
### These two functions store the data from the .yml file of image A and compare with the data of image B,C,D .yml files, in order to
### identify the same point detected in different images.
### The output is a list of dictinaries, for example: [[{pt1A},{pt1B},{pt1C},{pt1D}],[{pt2A},{pt2B},{pt2C},{pt2D}],[{pt3A},{pt3B},{pt3D}]]
### The code is automatized for more than 4 images, however it only works if the number of detections in image A (the first image) is the 
### maximum number of detections for all images. This works for our case, since image A has 11 detections, image B 11 detections, 
### image C 10 detections and image D 11 detections. Currently, I am working on this function in order to correct this standback

def find_corr_in_list(pto, im_lst):
	corr_pto=None
	for im_pto in im_lst:
		if abs(pto['Ra']-im_pto['Ra'])< 0.01:
			corr_pto=im_pto
			continue
	return corr_pto



def find_correspondence(im1,im2):
	lst_cors= []
	im_a_comparar=[]
	
	
	with open(f"data_image_{str(im1)}.yml", 'r') as stream:
		load1 = yaml.safe_load(stream)
		for i in load1:
			lst_cors.append([i])
	
	while os.path.exists(f"data_image_{str(im2)}.yml") == True:

		with open(f"data_image_{str(im2)}.yml", 'r') as stream:
			load2 = yaml.safe_load(stream)
			im_a_comparar.append(load2)
					
				
		im2=int(im2)	
		im2+=1
	
	
	for c in lst_cors:
		
		for im in im_a_comparar:
			
			ii = -1
			while c[ii] == {}:
				ii += -1
				
			pt_2_find_in_next_im= c[ii]
			
			new_pt= find_corr_in_list(pt_2_find_in_next_im,im)
			
			if new_pt == None:
				continue
			else:
				c.append(new_pt)
				
	return lst_cors

lst_cors=find_correspondence(1,2)  #This doesn't need to be altered in order to analyse all images
#print(lst_cors)

################################################################################################################################################
### This function converts Ra, Dec coordinates from the .yml file to Alt,Az coordinates into a .yml file. This code works for now
### but it will be updated 

def radec_to_altaz(im):

	with open(f"data_image_{im}.yml", 'r') as stream:
		load= yaml.safe_load(stream)

	ra = []
	dec = []
	for detection in load:
		ra.append(detection['Ra']*u.deg)
		dec.append(detection['Dec']*u.deg)


	location=EarthLocation(
			lat=40.182524*u.deg,
			lon=-7.873811*u.deg,
			height=820*u.m)


	time = datetime.datetime(2022, 4, 17,23,47,53)
	stars = SkyCoord(ra, dec, frame='icrs') #o frame corresponde à tranformação de coordenadas 
	stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))


	sky_coords = []
	for i in range(len(stars_altaz)):
		sky_coords.append([stars_altaz[i].alt.degree,stars_altaz[i].az.degree])

	results=[]
	for c in sky_coords:
			tmp = {}
			tmp['az'] = float(c[0])
			tmp['alt'] = float(c[1])
			results.append(tmp)

	with open(f"altaz_image{im}.yml", 'w') as outfile:
		yaml.dump(results, outfile, default_flow_style=False)


#altazfile=radec_to_altaz('1')  # If you want to convert the data from image 2, change to radec_to_altaz('2')


###################################################################################################################################################
### This function plots a graph ilustrating the movement of the same detection in the images in which this point was observed. 
###
def plot_detection_movement(detections):
	x=[]
	y=[]
	image=0
	for i in detections:
		x.clear()
		y.clear()
		for n in i:
			point_x=(n['Ra'])
			x.append(point_x)
			point_y=(n['Dec'])
			y.append(point_y)
		
		image+=1
		plt.plot(x, y,color='green', linestyle='dashed', linewidth = 3, marker='o', markerfacecolor='blue', markersize=12)

		plt.title(f'Objeto {image}')
		plt.xlabel('Ra (degrees)')
		plt.ylabel('Dec (degrees)')
		plt.xlim([min(x)-0.0001,max(x)+0.0001])
		plt.ylim([min(y)-0.0001,max(y)+0.0001])
		plt.show()

#desenho=plot_detection_movement(lst_cors)
