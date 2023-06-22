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
import cv2

PipelineDirectory = os.getcwd()
day = "2022-03-20"
DayDirectory = "%s/Files/Images/%s" %(PipelineDirectory,str(day))
config_file='ConfigurationDefesa.txt'
config = SPADEFunctions.ReadPipelineConfig(PipelineDirectory + "/Files/Config/" + config_file)
TLE_Directory = f"{DayDirectory}/TLE"
FileXY=f'{TLE_Directory}/XY_window.cat'
File='PAER_SURV_4_260753_230228_042347+065500_scan017_02.fit'


Images=['PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_02.fit',
'PAER_SURV_4_260753_230228_042347+065500_scan017_03.fit','PAER_SURV_4_260753_230228_042347+065500_scan017_04.fit']

#SPADEFunctions.Run_TLE_Code(Images,config,PipelineDirectory,DayDirectory)

def Drawing_Satellites_TLE_Position(DayDirectory,TLE_Directory):
	
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
	
	

	image='PAER_SURV_4_260753_230228_042347+065500_scan017_02.fit'
	json_file='PAER_SURV_4_260753_230228_042347+065500_scan017_Output.json'
	

	f = open(f"{TLE_Directory}/TLEOutput/{json_file}")
	data = json.load(f)
	x=[]
	y=[]
	for i in data:
		if i["Image"]== image:
			x.append(i['X (pixel) '])
			y.append(i['Y (pixel) '])
	    
	    
	fig = plt.figure()
	ax = fig.add_subplot(111) 
	vmin,vmax = None,None


	Imagem = CCDData.read(f"{DayDirectory}/ProcessedImages/{image}", unit = 'adu') 
	plt.imshow(Imagem,cmap='gray',vmin=vmin or np.nanpercentile(Imagem,3),vmax=vmax or np.nanpercentile(Imagem, 99))
	plt.title(f"{image.split('.')[0]}")
	#plt.xlim(1666,2042)
	#plt.ylim(755,1122)


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

#--------------------------------------------------------------------------------------------------------
#----------------------------Probabilistic Hough Transform ----------------------------------------------
#--------------------------------------------------------------------------------------------------------

def Binarization(image):
	
	image = CCDData.read(f"{image}", unit = 'adu') 
	image = np.uint8(image.data/256) #convert uint16 to uint8
	
	#Image enhancement
	clahe = cv2.createCLAHE(clipLimit=5, tileGridSize=(20, 20))  
	enhanced_image = clahe.apply(image)
	
	gaussian_filter = cv2.GaussianBlur(enhanced_image, (3, 3), 0)
	binary_image = cv2.adaptiveThreshold(gaussian_filter,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,-5)
	
	#----Delete objets with an area with less than 10 pixels---
	
	# perform connected component analysis
	num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(binary_image)
	# create an array to store the filtered image
	filtered_image = np.zeros_like(binary_image)
	
	# iterate over the connected components
	for label in range(1, num_labels):
	    # check the area of the component
	    area = stats[label, cv2.CC_STAT_AREA]
	    # filter out small objects
	    if area >= 10:
	        # add the component to the filtered image
	        filtered_image[labels == label] = 255
	
	return binary_image,enhanced_image,filtered_image

	
def HoughLines(DayDirectory,File):
	
	image=f'{DayDirectory}/ProcessedImages/{File}'
	binary_image,enhanced_image,filtered_image=Binarization(image)
	
	#----------- Perform edge detection---------------
	
	# We need two threshold values, minVal and maxVal. Any edges with intensity gradient more than maxVal are sure to be edges 
	# and those below minVal are sure to be non-edges, so discarded.
	T_lower=100 #Lower threshold value in Hysteresis Thresholding
	T_upper=200	#Upper threshold value in Hysteresis Thresholding
	
	aperture_size=3 #Aperture size of the Sobel filter.
	#The default value is 3 and its value should be odd between 3 and 7. 
	#You can increase the Aperture size when you want to detect more detailed features.
	edges = cv2.Canny(filtered_image, T_lower,T_upper, aperture_size)
	
	# ---------- Perform Probabilistic Hough Transform-------------
	
	lines = cv2.HoughLinesP(edges,rho=1,theta=np.pi/180,threshold=20,minLineLength=100,maxLineGap=10)
	#number_of_detected_lines=len(lines)
	
	# Extract line segment endpoints from Hough lines
	line_segments = []
	if lines is not None:
		for line in lines:
			x1, y1, x2, y2 = line[0]
			line_segments.append((x1, y1, x2, y2))
	#line_segments=[(1913, 878, 1922, 981), (1939, 930, 1952, 1085)]
	# Draw the dected lines in enhanced_image_bgr_1 image. The enhanced image had to be converted to rgb in order to draw a green line
	enhanced_image_bgr_1 = cv2.cvtColor(enhanced_image, cv2.COLOR_GRAY2BGR)
	for segment in line_segments:
	    x1, y1, x2, y2 = segment
	    cv2.line(enhanced_image_bgr_1, (x1, y1), (x2, y2), (0, 255, 0), 2)  
	
	
	# ----------Conditions to verify if two lines correspond to the same streak -----------
	
	#Condition 1- if the difference of the slope values is less than 1, it could be the same line
	
	#Get the slope and b values (y=mx+b)
	
	if len(line_segments) != 0:
		slope_values=[]
		b_values=[]
		
		for line in line_segments:
				x1=line[0]
				y1=line[1]
				x2=line[2]
				y2=line[3]
				min_x = min(x1, x2)
				max_x = max(x1, x2)
				min_y = min(y1, y2)
				max_y = max(y1, y2)
				
				slope=(max_y-min_y)/(max_x-min_x)
				b=max_y-slope*max_x
				
				slope_values.append(slope)
				b_values.append(b)
				
		
		#if the difference of the slope values is <1, it saves the data given by the hough tranform function in the same key "streak"
		#isto dá aqui mas está mal feito, tenho que pensar um bocado
		streak=[]
		for i in range(0,len(slope_values)-1):
			my_dict={}
			for j in range(i+1,len(slope_values)):
				if abs(slope_values[i]-slope_values[j]) <= 1:
					if "streak" in my_dict:
						my_dict["streak"].append(line_segments[i])  
						my_dict["streak"].append(line_segments[j])  
					else:
						my_dict["streak"] = [line_segments[i], line_segments[j]]
			streak.append(my_dict)
	
		x1=[]
		x2=[]
		y1=[]
		y2=[]
		for i in streak:
			for j in range(0,len(i['streak'])):
				x1.append(i['streak'][j][0])
				y1.append(i['streak'][j][1])
				x2.append(i['streak'][j][2])
				y2.append(i['streak'][j][3])
		min_x = min(min(x1, x2))
		max_x = max(max(x1, x2))
		min_y = min(min(y1, y2))
		max_y = max(max(y1, y2))
				
		enhanced_image_bgr_2 = cv2.cvtColor(enhanced_image, cv2.COLOR_GRAY2BGR)
		cv2.line(enhanced_image_bgr_2, (min_x, min_y), (max_x, max_y), (0, 255, 0), 2)  
			
		x_metade=min_x+(max_x-min_x)/2
		y_metade=min_y+(max_y-min_y)/2
		XYToRaDec(FileWc,FileXY,FileOutput)
		
	    
		images=[enhanced_image,binary_image,filtered_image,edges,enhanced_image_bgr_1,enhanced_image_bgr_2]
		titles = ['Enhanced Image', 'Binary Image', 'Remove <10 px objects','Edge Detection',
	             f'0 Hough Lines Dectected', 'Streak Detection']
		
		for i in range(6):
			plt.subplot(2,3,i+1),plt.imshow(images[i],'gray')
			plt.title(titles[i])
			plt.xticks([]),plt.yticks([])
			plt.xlim(1800,1900)
			plt.ylim(0,150)
			
	else:
		images=[enhanced_image,binary_image,filtered_image,edges,enhanced_image_bgr_1]
		titles = ['Enhanced Image', 'Binary Image', 'Remove <10 px objects','Edge Detection',
	             f'Hough Lines Dectected']
		
		for i in range(5):
			plt.subplot(2,3,i+1),plt.imshow(images[i],'gray')
			plt.title(titles[i])
			plt.xticks([]),plt.yticks([])
			plt.xlim(1800,1900)
			plt.ylim(0,150)
			plt.ylim(0,150)
	plt.show()
	
streak_detection=HoughLines(DayDirectory,File)


#--------------------------------------------------------------------------------------------------------
#----------------------------Binarization Study ----------------------------------------------
#--------------------------------------------------------------------------------------------------------

def Binarization_Study(enhanced_image):
	

	gaussian_filter = cv2.GaussianBlur(enhanced_image, (3, 3), 0)
	
	#----------------------Global Thresh-----------------------------------------------
	# threshold_value, binary_image1 = cv2.threshold(gaussian_filter, 45, 1, cv2.THRESH_BINARY)
	#threshold_value, binary_image2 = cv2.threshold(gaussian_filter, 50,1, cv2.THRESH_BINARY)
	# threshold_value, binary_image3 = cv2.threshold(gaussian_filter, 55, 1, cv2.THRESH_BINARY)
	# titles = ['Global Thresh (c=45)',
             # 'Global Thresh (c=50)', 'Global Thresh (c=55)']
	# images = [ binary_image1 ,binary_image2, binary_image3]

	#Adaptative Guassian Thresh
	# thresh_gaussian1 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -2)
	#thresh_gaussian2 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -3)
	# thresh_gaussian3 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -4)
	# titles = ['Adaptive Gaussian Thresh (block=11,c=-2)',
             # 'Adaptive Gaussian Thresh (block=11,c=-3)', 'Adaptive Gaussian Thresh (block=11,c=-4)']
	# images = [ thresh_gaussian1 ,thresh_gaussian2, thresh_gaussian3]

	##Adaptative Thresh Mean 
	thresh_mean1 = cv2.adaptiveThreshold(gaussian_filter,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,-4)
	thresh_mean2 = cv2.adaptiveThreshold(gaussian_filter,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,-5)
	thresh_mean3 = cv2.adaptiveThreshold(gaussian_filter,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,-6)
	titles = ['Adaptive Mean Thresh (block=11,c=-4)',
             'Adaptive Mean Thresh (block=11,c=-5)', 'Adaptive Mean Thresh (block=11,c=-6)']
	images = [ thresh_mean1 ,thresh_mean2, thresh_mean3]

	#Image filtering
	#_, thresh_otsu_plusfilter = cv2.threshold(gaussian_filter,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
	 
	
	
	# titles = ['Global Thresh (c=50) ', 'Mean Adaptative Thresh',
            # 'Gaussian Adaptative Thresh', 'Otsu Thresh']
            
	# images = [binary_image2, thresh_mean2 , thresh_gaussian2,thresh_otsu_plusfilter]
	for i in range(3):
		plt.subplot(1,3,i+1),plt.imshow(images[i],'gray')
		plt.title(titles[i])
		plt.xticks([]),plt.yticks([])
		plt.xlim(1666,2042)
		plt.ylim(755,1122)
	plt.show()

#--------------------------------------------------------------------------------------------------------
#------------------------------------------ Grafico -----------------------------------------------------
#--------------------------------------------------------------------------------------------------------

TLE_file=f'TLE Catalog-2023228.json'
WindowFile=['PAER_SURV_4_260789_230228_080700-051142_scan102_01_RaDec_Window','PAER_SURV_4_260789_230228_080700-051142_scan102_02_RaDec_Window',
'PAER_SURV_4_260789_230228_080700-051142_scan102_03_RaDec_Window','PAER_SURV_4_260789_230228_080700-051142_scan102_04_RaDec_Window']

def HMS2deg(RA):
	''' 
	Function created to transform HMS to degrees
	
	INPUTS:
	Ra 					-> Ra value in HMS												(str)
	
	OUTPUTS:
	degrees				-> degrees value corresponding to input							(float)
	
	FILES CREATED:
	NONE
	'''
	
	degrees = []
	for value in RA:
		H = float(value[1:2])
		M = float(value[3:5])
		S = float(value[6:8])
		degree = 360/24 * H + 360/(24*60) * M + 360/(24*60*60) * S
		degrees.append(degree)
		
	return degrees

def DMS2deg(DEC):
	''' 
	Function created to transform DMS to degrees
	
	INPUTS:
	Dec 				-> Dec value in DMS												(str)
	
	OUTPUTS:
	degrees				-> degrees value corresponding to input							(float)
	
	FILES CREATED:
	NONE
	'''
	
	degrees = []
	for value in DEC:
		D = float(value[0:2])
		M = float(value[3:5])
		S = float(value[6:8])
		degree = D + M/60 + S/3600
		degrees.append(degree)
		
	return degrees
		
def convert_degrees_to_hours_minutes(degrees):
    total_minutes = degrees * 4  # Since 15 degrees = 1 hour, 1 degree = 4 minutes
    hours = int(total_minutes // 60)
    minutes = int(total_minutes % 60)
    string=str(hours)+':00'
    return string

def TLE_Graph(config,DayDirectory,TLE_file,TLE_Directory):
		
	file_contents=SPADEFunctions.TLE_parameters(TLE_Directory,TLE_file)
	print(len(file_contents))
	ts = load.timescale()
	t_inicial=ts.utc(2023,2,28,19,59)
	t_final=ts.utc(2023,2,28,20,21,49) 

	Ra=[]
	Dec=[]	
	print('aaa')
	count=0
	for i in file_contents:
			count=count+1
			print(f'{count}')
			line1=i['tle_1']
			line2=i['tle_2']
			nome='n' #o nome não interessa
			satellite = EarthSatellite(line1, line2, nome, ts)

			bluffton = wgs84.latlon(float(config['Latitude']), float(config['Longitude']),float(config['Height']))
			
			difference = satellite - bluffton 
			
			passo=(t_final-t_inicial)/60
			for t in np.arange(t_inicial,t_final,passo):
			
				topocentric = difference.at(t)
		
				ra, dec, distance = topocentric.radec()
				
				ra=math.degrees(ra.radians)
				Ra.append(ra)
		
				dec=math.degrees(dec.radians)
				Dec.append(dec)


	plt.scatter(Ra, Dec,alpha=0.5,s=5)
	print('ccc')

	# ra_target=['04 23 47.00','04 23 47.00','04 23 47.00','04 23 47.00','04 23 47.00','08 07 00.00','08 07 00.00','08 07 00.00','08 07 00.00']
	# dec_target=['-3 13 45.0','+6 55 00.0','-0 41 33.0','+1 50 37.0','-8 18 07.0','-7 43 53.0','+4 57 03.0','-5 11 42.0','+2 24 52.0']
	# ra_degrees= HMS2deg(ra_target) 
	# dec_degrees= DMS2deg(dec_target)
	#plt.scatter(ra_degrees, dec_degrees,marker='s',s=40)
	
	
	files=os.listdir(f"{DayDirectory}/TLE")
	print('aaa')
	for ficheiro in files:
		ra_window=[]
		dec_window=[]
		if ficheiro[-12:] == 'RaDec_Window':
			A=fits.open(f'{DayDirectory}/TLE/{ficheiro}')
			for i in A[1].data['RA']:
				ra_window.append(float(i))
			for i in A[1].data['DEC']:
				dec_window.append(float(i))
			
			min_x = min(ra_window)
			max_x = max(ra_window)
			min_y = min(dec_window)
			max_y = max(dec_window)

			square_x = [min_x, max_x, max_x, min_x, min_x]  # Vertices of the square
			square_y = [min_y, min_y, max_y, max_y, min_y]
			
			plt.plot(square_x, square_y, 'r-') 
	
	# ra_window=[]
	# dec_window=[]
	# for ficheiro in WindowFile:
		# A=fits.open(f'{TLE_Directory}/{ficheiro}')
		# ra_window.append(float(A[1].data['RA']))
		# dec_window.append(float(A[1].data['DEC']))
	
	# plt.scatter(ra_degrees, dec_degrees,marker='s',s=40)
	

	

	


	plt.xlabel('Ra (Hours) ',fontsize=20)
	plt.ylabel('Dec (Degrees)',fontsize=20)

	angle=[0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360]
	x=[]
	for i in angle:
		hms=convert_degrees_to_hours_minutes(i)
		x.append(hms)
	
	plt.xticks(angle, x,fontsize=15)
	
	plt.xticks(angle, x,fontsize=15)

	plt.xlim([0, 359])
	plt.ylim([-20, 20])
	plt.show()

#TLE_Graph(config,DayDirectory,TLE_file,TLE_Directory)
