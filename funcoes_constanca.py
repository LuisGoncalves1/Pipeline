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
scale_arcsec=4.7
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
	The FWHM is the mean FWHM associated to each image 
	
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


def criteria(pto, im_lst,scale_arcsec):
	
	corr_pto=None
	
	for im_pto in im_lst:

		L1=[float(pto[0]['RaCalculated']),float(pto[0]['DecCalculated'])]
		L2=[float(im_pto['RaCalculated']),float(im_pto['DecCalculated'])]
		distance=AngleDifferencePoints(L1,L2)
	
		fwhm=[float(pto[0]["FWHM"]),float(im_pto['FWHM'])]

		if float(distance) <= (scale_arcsec*max(fwhm))/3600:
			corr_pto=im_pto	
			continue

	return corr_pto

def find_correspondence(fix_list,image,scale_arcsec,next_image):
	for c in fix_list:
		  for i in range(next_image,len(image)): 
			  for im in image[i]:
				  pt_2_find_in_next_im = c
				  new_pt= criteria(pt_2_find_in_next_im,im,scale_arcsec) 
				  if new_pt == None:
					  continue
				  else:
					  c.append(new_pt)
	return fix_list

def associate_images(DayDirectory,all_files,scale_arcsec):
	
	'''	
	This function is a list of dictionaries with all detections from all images. If the distance between two point 
	(calculated by AngleDifferencePoints function) is less than a specific value, it means that the detections correspond
	to the same object in different images, so they're stored in the same list.
	It's important to notice that the algorithm is currently treating every object has a potencial space debris detection,
	i.e the stars are not being distinguished from the rest.
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
	
	obj = [[] for i in range(0,len(all_files))]
	print(f'[+] Start analysing data of {len(obj)} images.')
	

	count=0
	for ficheiro in all_files:
		with open(f"{DayDirectory + '/OutputFiles/' + ficheiro}.json", 'r') as stream:
			jsonfile = json.load(stream)
			for i in jsonfile:
				obj[count].append([i])
			count+=1
	
	for i in range(0,len(all_files)):
		print(f'Image{all_files[i][-1]} has {len(obj[i])} objects.')
		
	#obj=sorted(obj, key=len,reverse=True)
	next_image=1
	correspondences=obj[0]
	print(len(correspondences))
	correspondences=find_correspondence(correspondences,obj,scale_arcsec,next_image)
	print(len(correspondences))

	
	while next_image<len(obj):
		
		idx=[]
		img=obj[next_image]
		for d in img: 
			for c in correspondences:
				for i in range (0,len(c)):
					if c[i]['image'] == d[0]['image'] and c[i]['DecCalculated'] == d[0]['DecCalculated']: #significa que está
						index = obj[next_image].index(d)
						idx.append(index)
					

		for index in sorted(idx, reverse=True): 
			del img[index]
		
		
		if next_image == len(obj)-1:
			correspondences.extend(img)
			next_image+=1
		else:
			next_image+=1
			cors=find_correspondence(img,obj,scale_arcsec,next_image)
			correspondences.extend(cors)
			
	return correspondences

		
associate_images(DayDirectory,all_files,scale_arcsec)
	


################################################################################################################################
#														Scattered plots
################################################################################################################################


def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return round(mult*deg), round(mult*mnt), float("{:.1f}".format(mult*sec))
	
	
def alldetections(lst1,lst2):

	number_detections_star=[]
	number_detections_nonidnt_objs=[]
	
	colors_stars=[]
	colors_nonidnt_objs=[]
	

	for i in lst1:
		number_detections_star.append(len(i))

	for i in number_detections_star:
			if i==1:
				black=(0.2,0.2,0.2)
				colors_stars.append(black)
			if i==2:
				yellow=(1, 1, 0)
				colors_stars.append(yellow)
			if i==3:
				blue=(0, 0.5, 1)
				colors_stars.append(blue)
			if i==4:
				red=(1,0, 0)
				colors_stars.append(red)

		
	
	
	for i in lst2:
		number_detections_nonidnt_objs.append(len(i))
	
	for i in number_detections_nonidnt_objs:
			if i==1:
				black=(0.2,0.2,0.2)
				colors_nonidnt_objs.append(black)
			if i==2:
				yellow=(1, 1, 0)
				colors_nonidnt_objs.append(yellow)
			if i==3:
				blue=(0, 0.5, 1)
				colors_nonidnt_objs.append(blue)
			if i==4:
				red=(1,0, 0)
				colors_nonidnt_objs.append(red)
			
	x_stars=[]
	y_stars=[]
	for i in lst1:
		x_stars.append(float(i[0]['Ra_calculated']))

	

	for i in lst1:
		y_stars.append(float(i[0]['Dec_calculated']))
	
	
	
	x_nonidnt_objs=[]
	y_nonidnt_objs=[]
	for i in lst2:
		x_nonidnt_objs.append(float(i[0]['Ra_calculated']))



	for i in lst2:
		y_nonidnt_objs.append(float(i[0]['Dec_calculated']))
		
	return x_stars,y_stars,x_nonidnt_objs,y_nonidnt_objs,colors_stars,colors_nonidnt_objs


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


def scatteredplot_colors():
	
	alldect=alldetections(lst_cors_stars,lst_cors_nonidnt_objs)
	
	
	#scatter=plt.subplot(1, 2, 1)
	plt.scatter(
	x=alldect[0],
	y=alldect[1],
	c=alldect[4],
	alpha=0.5, 
	)
	
	my_xticks=[]
	x=alldect[0]
	x.sort()
	for i in x:
		my_xticks.append(decdeg2dms(i))


	plt.xticks([x[0], x[-1]], [my_xticks[0], my_xticks[-1]],rotation=0,fontsize=8)
	
	my_yticks=[]
	y=alldect[1]
	y.sort()
	for i in y:
		my_yticks.append(decdeg2dms(i))
	plt.yticks([y[0], y[-1]], [my_yticks[0], my_yticks[-1]],rotation=0,fontsize=8)
	

	#plt.title("All Stars")
	#plt.ylabel('Dec (dms)')
	#plt.xlabel('Ra (dms)')

	#plt.subplot(1, 2, 2)
	plt.scatter(
		x=alldect[2],
		y=alldect[3],
		c=alldect[5],
		alpha=0.5,
	)
	
	my_xticks=[]
	x=alldect[2]
	x.sort()
	for i in x:
		my_xticks.append(decdeg2dms(i))


	plt.xticks([x[0], x[-1]], [my_xticks[0], my_xticks[-1]],rotation=0,fontsize=8)
	
	my_yticks=[]
	y=alldect[3]
	y.sort()
	for i in y:
		my_yticks.append(decdeg2dms(i))
	plt.yticks([y[0], y[-1]], [my_yticks[0], my_yticks[-1]],rotation=0,fontsize=8)
	
	#plt.title("Non-identified objects (black:1 , yellow:2 , blue:3 , red:4)")
	plt.title("All Objects")
	plt.ylabel('Dec (dms)')
	plt.xlabel('Ra (dms)')
	
	plt.show()
	

			
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


################################################################################################################################
#								RASCUNHOS : NÃO COLOCAR NO SPADEFunctions.py AINDA
################################################################################################################################
