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


################################################################################################################################
### This function reads the output file (stored in the OutputFiles folder), and stores the (x,y,Ra,Dec) coordinates of the 
### detections in a .yml file (list of dicts, in which the value is a float). Some parameters in this function depend on whether
### you want to export the data for the detections that were matched as a star (match=1) or as an unidentified object (match=0).
### There are comments in the function to tell you which parameters are these.
### This function is only automatized for 4 images.In time I will improve the function and it will be automatized for whatever 
### set of images.

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
			if float(resultados[2]) == 0:  #put 0 for the unidentified objects and 1 for stars
				tmp = {}
				tmp['x'] = float(resultados[0])
				tmp['y'] = float(resultados[1])
				tmp['Ra_calculated'] = float(resultados[3])
				tmp['Dec_calculated'] = float(resultados[4])
				#tmp['Ra_catalogue'] = float(resultados[5]) 			 #If you are analysing stars, you need the (Ra_catalogue, Dec_catalogue) parameters. However if you're analysing
				 #unidentified objects, comment the Ra_catalogue and  Dec_catalogue, since they only apply to the stars 
				 #from the catalogue.
				#tmp['Dec_catalogue'] = float(resultados[6])
				items_2_store.append(tmp)
	
		with open(f"output_matches=0_image{count}.yml", 'w') as outfile: #change the name of the file to matches=0 or matches=1
			yaml.dump(items_2_store, outfile, default_flow_style=False)
			
					
		ficheiro=ficheiro[:-1]
		count+=1
		
#output=convert_output_file() 

##################################################################################################################################
### This function exports the data given by the Sextractor in the .cat file to a .txt file. This data is associated with
### all detections (stars and unidentified objects), because its the data obtained before passing through astrometry. It gives
### you the values of (x,y,fwhm,SNR, among others) of all the detections in each image.
### The code is still not automatized, so you need to make changes in the function for each image.

def cattotxt():

	diretoria='./Files/Images/2022-03-20/SextractorCatalogue/'
	ficheiro='MAPA_SURV_4_085283_220417_114004+060023_scan047_04.cat' #change the name of the cat file, depending on the image you are analysing
	filename= diretoria + ficheiro


	from astropy.io import fits
	A=fits.open(filename)
	x=A[1].data['X_IMAGE']
	y=A[1].data['Y_IMAGE']
	fwhm=A[1].data['FWHM_IMAGE']
	SNR=A[1].data['SNR_WIN']

	output='SEx_Output_4.txt' #change the name of the ouput file depending on the image you are analysing

	hdulist = fits.open(filename)
	scidata = hdulist[1].data

	np.savetxt(output, scidata, fmt='%10.10f', delimiter=',')
	
#a=cattotxt()

###################################################################################################################################
### This function analyses the .txt file created by the cattxt() function,filters the data by a given SNR, and stores the result
### in a .yml file. Currently, it's storing the data in which 50<SNR<100. 
### This function is only automatized for 4 images.In time I will improve the function and it will be automatized for whatever 
### set of images.

def ficheiroestrelas():
	
		directory = os.getcwd()
		count=1
		
		while count<5:
			ficheiro=f'/SEx_Output_{count}.txt'

			filepath_2=directory + ficheiro

			with open(filepath_2,'r') as fp_2:
				all_the_lines_2= fp_2.readlines()
				items_sexdata=[]
				for i in all_the_lines_2:
					items_sexdata.append(i)

			items_2_store_sex=[]	
	
			for i in items_sexdata:
			
				resultados_2=i.split(',')
				
				if float(resultados_2[12])>100 and float(resultados_2[12])<150:
				
					tmp = {}
					tmp['x'] = float(resultados_2[6])
					tmp['y'] = float(resultados_2[7])
					items_2_store_sex.append(tmp)
		
			with open(f"SEx_data_SNRfiltered_image{count}.yml", 'w') as outfile:
				yaml.dump(items_2_store_sex, outfile, default_flow_style=False)
		
			count+=1

#a=ficheiroestrelas()


#def scatterstars():
#	x_list=[]
#
#	with open(f"SEx_data_SNRfiltered_image1.yml", 'r') as stream:
#		load = yaml.safe_load(stream)
#		for i in load:
#			x_list.append([i])
		
#final_lst=[]

#for i in range(0,len(x_list)-1):
#	if abs(x_list[i][0]['x']-x_list[i+1][0]['x'])>1000:
#		final_lst.append(x_list[i])
		
#print(final_lst)




################################################################################################################################
### This function compares the data from the sextractor file (SNR filtered)  with the out

def stardata():
	
	count=1
	while count<5:
	
	
		lst_sex=[]
		with open(f"SEx_data_SNRfiltered_image{count}.yml", 'r') as stream:
			load1 = yaml.safe_load(stream)
			for i in load1:
				lst_sex.append([i])
		lst_out=[]
		
		with open(f"output_matches=1_image{count}.yml", 'r') as stream:
			load2 = yaml.safe_load(stream)
			for i in load2:
				lst_out.append([i])

		result=[]	
		for a in lst_out:
			for b in lst_sex:
				if round(a[0]['x'],4)==round(b[0]['x'],4):
					items={}
					items['x']=a[0]['x']
					items['y']=a[0]['y']
					items['Ra_calculated']=a[0]['Ra_calculated']
					items['Dec_calculated']=a[0]['Dec_calculated']
					items['Ra_catalogue']=a[0]['Ra_catalogue']
					items['Dec_catalogue']=a[0]['Dec_catalogue']
					result.append(items)

			
		with open(f"star_data_image{count}.yml", 'w') as outfile:
			yaml.dump(result, outfile, default_flow_style=False)
		
		count+=1
#a=stardata()


###############################################################################################################################
### This function converts Ra, Dec coordinates from the .yml file (output_matches=1_image1.yml or output_matches=0_image1.yml
### depending on wheter you're analysing stars or unidentified objects). to Alt,Az coordinates into a .yml file. This code works for now
### but it will be updated 

def radec_to_altaz(im):

	with open(f"output_matches=0_image{im}.yml", 'r') as stream:
		load= yaml.safe_load(stream)

	ra = []
	dec = []
	for detection in load:
		ra.append(detection['Ra_calculated']*u.deg)
		dec.append(detection['Dec_calculated']*u.deg)
		
	print(ra)


	location=EarthLocation(
			lat=40.182524*u.deg,
			lon=-7.873811*u.deg,
			height=820*u.m)


	time = datetime.datetime(2022, 4, 17,23,48,8)  # you have to change the date, according to the image
	stars = SkyCoord(ra, dec, frame='icrs') #o frame corresponde à tranformação de coordenadas 
	stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))


	sky_coords = []
	for i in range(len(stars_altaz)):
		sky_coords.append([stars_altaz[i].alt.degree,stars_altaz[i].az.degree])

	results=[]
	for c in sky_coords:
			tmp = {}
			tmp['alt'] = float(c[0])
			tmp['az'] = float(c[1])
			results.append(tmp)

	with open(f"altaz_matches=0_image{im}.yml", 'w') as outfile:
		yaml.dump(results, outfile, default_flow_style=False)


#altazfile=radec_to_altaz('4')  # If you want to convert the data from image 2, change to radec_to_altaz('2')


##################################################################################################################################

def allcoordinates(count,match):

	radecdata=[]
	altazdata=[]
	alldata=[]
	
	
	with open(f"output_matches=0_image{count}.yml", 'r') as stream:
		load1 = yaml.safe_load(stream)
		for i in load1:
			radecdata.append([i])


	with open(f"altaz_matches={match}_image{count}.yml", 'r') as stream:
		load2 = yaml.safe_load(stream)
		for i in load2:
			altazdata.append([i])

		

	for i in range(0,len(radecdata)):
		radecdata[i].extend(altazdata[i])
		
	for element in radecdata:
		resultados={}
		resultados={**element[0],**element[1]}
		alldata.append(resultados)

	
	with open(f"alldata_matches={match}_image{count}.yml", 'w') as outfile:
			yaml.dump(alldata, outfile, default_flow_style=False)
			
	return alldata
			
#a=allcoordinates('4','0')
#print(a)
#########################################################################################################################################
### These two functions store the data from the .yml file of image A and compare with the data of image B,C,D .yml files, in order to
### identify the same point detected in different images.
### The output is a list of dictionaries, for example: [[{pt1A},{pt1B},{pt1C},{pt1D}],[{pt2A},{pt2B},{pt2C},{pt2D}],[{pt3A},{pt3B},{pt3D}]]
### The code is automatized for more than 4 images, however it only works if the number of detections in image A (the first image) is the 
### maximum number of detections for all images. This works for our case, since image A has 11 detections, image B 11 detections, 
### image C 10 detections and image D 11 detections. Currently, I am working on this function in order to correct this standback

def find_corr_in_list(pto, im_lst):
	corr_pto=None
	for im_pto in im_lst:		
		if pto['Dec_catalogue'] == im_pto['Dec_catalogue'] and pto['Ra_catalogue'] == im_pto['Ra_catalogue']:
		
		#if abs(pto['Dec_calculated']-im_pto['Dec_calculated'])<=(75/3600) and abs(pto['Ra_calculated']-im_pto['Ra_calculated'])<=(75/3600):
			corr_pto=im_pto
			continue
	return corr_pto



def find_correspondence(im1,im2,count):
	lst_cors= []
	im_a_comparar=[]
	
	
	with open(f"alldata_matches={count}_image{str(im1)}.yml", 'r') as stream:
		load1 = yaml.safe_load(stream)
		for i in load1:
			lst_cors.append([i])
	
	while os.path.exists(f"alldata_matches={count}_image{str(im2)}.yml") == True:

		with open(f"alldata_matches={count}_image{str(im2)}.yml", 'r') as stream:
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


#lst_cors=find_correspondence(1,2,1)  #This doesn't need to be altered in order to analyse all images
#print(len(lst_cors))

#correspondencias=[]
#for i in lst_cors:
#	if len(i)==4:
#		correspondencias.append(i)

#print(len(correspondencias))

############################################################################################################
#lst=[]
#for i in range(0,len(correspondencias)-1):
#	if abs(correspondencias[i][0]['x']-correspondencias[i+1][0]['x'])>500 and abs(correspondencias[i][0]['y']-correspondencias[i+1][0]['y'])>100:
#		lst.append(correspondencias[i])
		
		
		
#for i,a in enumerate(lst):
#	print(a)
#	print(i)
#	print("\n")
###############################################################################################################################
### This function converts degrees to (degrees,minutes,seconds)

def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return round(mult*deg), round(mult*mnt), float("{:.3f}".format(mult*sec))

dd = 45 + 30/60 + 1/3600


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

#az=[]	
#alt=[]
#for i in lst_cors:
#	for n in i:
#		az.append(float(n['az']))
#		alt.append(float(n['alt']))
		
#desenho=plot_detection_movement(lst)


##################################################################################################################################



#input: lista dos dicionarios dos pontos a plotar
#im1=1  #vai variar dependendo do ficheiro yml que queremos ir buscar


#with open(f"star_data_image{im1}.yml", 'r') as stream:
#	load = yaml.safe_load(stream)

def plot_point_dict_lst():
	
	PipelineDirectory = os.getcwd()
	day = "2022-03-20"
	DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day)) #diretoria
	ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits' #vai variar dependendo da imagem que queremos representar no background do plot
			
	fig = plt.figure()
	ax = fig.add_subplot(111) 


	Imagem = CCDData.read(DayDirectory  + "/ReducedImages/" + ficheiro, unit = 'adu') 
	plt.imshow(Imagem,cmap='gray')
	plt.title(f'Imagem{im1}')
	

	for i in lst:
		c1 = plt.Circle((i[0]['x'],i[0]['y']),4, color = 'red', linewidth=2,fill=False)
		ax.add_patch(c1)
		
	plt.show()
	
	return True

#plot_point_dict_lst()


