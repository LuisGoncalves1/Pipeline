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


#--------------------------------------------------------------------------------------------------------
#----------------------------Probabilistic Hough Transform ----------------------------------------------
#--------------------------------------------------------------------------------------------------------
'''
The HoughLineP function in OpenCV performs probabilistic Hough line detection, which finds line segments in an image by considering only a subset of
the possible parameter space, rather than the entire space. It works as follows:

    Edge detection: The input image is typically pre-processed to extract edges using techniques such as Canny edge detection. 
					This is because the Hough transform is sensitive to edges, and edge points are used as the input for line detection.

    Line segment estimation: Random subsets of edge points are selected, and the corresponding lines in the parameter space are estimated. 
							 This is done by calculating the (ρ, θ) values for the lines passing through the selected edge points.

    Accumulation: The estimated lines are accumulated in a Hough accumulator, which is a 2D array representing the parameter space. 
				  Each cell in the accumulator represents a possible line in the image.

    Thresholding: After the accumulation, the accumulator is thresholded to identify the most prominent lines that have enough votes.
				  The threshold is set using the threshold parameter in the HoughLineP function.

    Line extraction: Line segments are extracted from the accumulator based on the thresholded votes.
					 For each line segment, the corresponding (ρ, θ) values are converted back to Cartesian coordinates (x, y) to obtain the start and end points of the line segment.

    Optional line filtering: Additional filtering can be applied, such as discarding lines that are too short or have large gaps between segments, 
							using the minLineLength and maxLineGap parameters in the HoughLineP function.


Parameters of HoughLineP Function

image -  the input image on whiche the hough line detection is perfomed. It should be a single-channel grayscale image

rho - This is the distance resolution of the accumulator in pixels. It determines the accuracy of the detection. 
	  Smaller values of rho result in more accurate detection but may increase computation time.

theta- This is the angular resolution of the accumulator in radians. It determines the angular accuracy of the detection.
	   Smaller values of theta result in more accurate detection but may increase computation time.

threshold - This is the minimum number of votes (or intersections in the accumulator) required for a line to be considered as a valid line. 
			It is a threshold value that determines the line detection sensitivity. Smaller values of threshold result in more lines being detected,
			but may also result in more false positives.
			
minLineLength - This is the minimum length of a line segment to be considered as a valid line. Lines shorter than this length will be discarded. 
				It is specified in pixels.
				
maxLineGap - This is the maximum gap allowed between segments to be merged into a single line. Lines with gaps greater than this value will not be merged. 
			It is specified in pixels.
'''

from skimage import io, filters, morphology

def Histogram_Equalization():

	image='PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit'
	image = CCDData.read(f"{image}", unit = 'adu') 
	image = np.uint8(image.data/256) #convert uint16 to uint8
	
	
	#Image enhancement
	clahe = cv2.createCLAHE(clipLimit=5, tileGridSize=(20, 20))  
	enhanced_image = clahe.apply(image)
	#enhanced_image=cv2.equalizeHist(image)
	
	# Calculate the histogram
	hist_image = cv2.calcHist([image], [0], None, [256], [0, 256])
	hist_enhanced_image = cv2.calcHist([enhanced_image], [0], None, [256], [0, 256])
	
	
	plt.subplot(2, 2, 1)
	plt.imshow(image, cmap='gray')
	plt.title('Parte of Original Image')
	plt.xlim(1666,2042)
	plt.ylim(755,1122)
	
	plt.subplot(2, 2, 2)
	plt.title('Grayscale Histogram')
	plt.xlabel('Pixel Intensity')
	plt.ylabel('Frequency')
	plt.plot(hist_image)
	plt.xlim(0,255)
	plt.ylim(0,1100000)
	
	plt.subplot(2, 2, 3)
	plt.imshow(enhanced_image, cmap='gray')
	plt.title('Part of Enhanced Image')
	plt.xlim(1666,2042)
	plt.ylim(755,1122)
	
	
	plt.subplot(2, 2, 4)
	plt.title('Grayscale Histogram (Block=20, ClipLimit=5)')
	plt.xlabel('Pixel Intensity')
	plt.ylabel('Frequency')
	plt.plot(hist_enhanced_image)
	plt.xlim(0,255)
	plt.ylim(0,1100000)
	

	#plt.show()
	
	return enhanced_image	

def Binarization_Study(gaussian_filter):
	
	#----------------------Global Thresh-----------------------------------------------
	threshold_value, binary_image1 = cv2.threshold(gaussian_filter, 45, 1, cv2.THRESH_BINARY)
	threshold_value, binary_image2 = cv2.threshold(gaussian_filter, 50,1, cv2.THRESH_BINARY)
	threshold_value, binary_image3 = cv2.threshold(gaussian_filter, 55, 1, cv2.THRESH_BINARY)
	titles = ['Global Thresh (c=45)',
             'Global Thresh (c=50)', 'Global Thresh (c=55)']
	images = [ binary_image1 ,binary_image2, binary_image3]

	
	#Adaptative Guassian Thresh
	# thresh_gaussian1 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -2)
	# thresh_gaussian2 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -3)
	# thresh_gaussian3 = cv2.adaptiveThreshold(gaussian_filter, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,11, -4)
	# titles = ['Adaptive Gaussian Thresh (block=11,c=-2)',
             # 'Adaptive Gaussian Thresh (block=11,c=-3)', 'Adaptive Gaussian Thresh (block=11,c=-4)']
	# images = [ thresh_gaussian1 ,thresh_gaussian2, thresh_gaussian3]

	
	##Adaptative Thresh Mean 
	thresh_mean = cv2.adaptiveThreshold(gaussian_filter,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,11,-3)

	
	
	
	#Image filtering
	#_, thresh_otsu = cv2.threshold(enhanced_image,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
	_, thresh_otsu_plusfilter = cv2.threshold(gaussian_filter,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
	 
	
	
	# titles = ['Global Thresh (c=50) + G.Filter ', 'Mean Adaptative Thresh+ G.Filter',
            # 'Gaussian Adaptative Thresh+ G.Filter', 'Otsu Thresh+ G.Filter']
	# images = [binary_image2, thresh_mean , thresh_gaussian,thresh_otsu_plusfilter]
	for i in range(3):
		plt.subplot(1,3,i+1),plt.imshow(images[i],'gray')
		plt.title(titles[i])
		plt.xticks([]),plt.yticks([])
		plt.xlim(1666,2042)
		plt.ylim(755,1122)
	plt.show()

def Binarization():
	
	enhanced_image	= Histogram_Equalization()
	gaussian_filter = cv2.GaussianBlur(enhanced_image, (3, 3), 0)
	hist_image_filtered = cv2.calcHist([gaussian_filter], [0], None, [256], [0, 256])
	hist_image_enhanced = cv2.calcHist([enhanced_image], [0], None, [256], [0, 256])
	
	plt.subplot(2, 2, 1)
	plt.imshow(enhanced_image, cmap='gray')
	plt.title('Part of Enhanced Image')
	plt.xlim(1666,2042)
	plt.ylim(755,1122)
	
	plt.subplot(2, 2, 2)
	plt.title('Grayscale Histogram of the Full Image')
	plt.xlabel('Pixel Intensity')
	plt.ylabel('Frequency')
	plt.plot(hist_image_enhanced)
	plt.xlim(0,255)
	plt.ylim(0,1100000)
	
	plt.subplot(2, 2, 3)
	plt.imshow(gaussian_filter, cmap='gray')
	plt.title('Part of Enhanced Image + Gaussian Filter')
	plt.xlim(1666,2042)
	plt.ylim(755,1122)
	
	plt.subplot(2, 2, 4)
	plt.title('Grayscale Histogram of the Full Image')
	plt.xlabel('Pixel Intensity')
	plt.ylabel('Frequency')
	plt.plot(hist_image_filtered)
	plt.xlim(0,255)
	plt.ylim(0,1100000)
	#plt.show()
	
	binary_image=Binarization_Study(gaussian_filter)
	
Binarization()

	# cv2.namedWindow('Gaussian', cv2.WINDOW_NORMAL)  
	# cv2.imshow('Gaussian', thresh_gaussian)
	# screen_width, screen_height = 1090, 1080  
	# cv2.resizeWindow('Gaussian', screen_width, screen_height)
	# cv2.waitKey(0)
	# cv2.destroyAllWindows()
	
	
	# cv2.namedWindow('Mean', cv2.WINDOW_NORMAL)  
	# cv2.imshow('Mean', thresh_mean)
	# screen_width, screen_height = 1090, 1080  
	# cv2.resizeWindow('Mean', screen_width, screen_height)
	# cv2.waitKey(0)
	# cv2.destroyAllWindows()
	
	
	
	# cv2.namedWindow('otsu', cv2.WINDOW_NORMAL)  
	# cv2.imshow('otsu', thresh_otsu)
	# screen_width, screen_height = 1090, 1080  
	# cv2.resizeWindow('otsu', screen_width, screen_height)
	# cv2.waitKey(0)
	# cv2.destroyAllWindows()
		
# def binarization(image):
	
	# #Image enhancement
	# clahe = cv2.createCLAHE(clipLimit=5, tileGridSize=(20, 20))  
	# enhanced_image = clahe.apply(image)
	
	
	
	# # fig, axs = plt.subplots(1, 2,figsize=(10, 5))
	# # axs[0].imshow(image, cmap='gray')
	# # axs[0].set_title('Original Image')
	# # axs[0].axis('off')
	# # axs[1].imshow(enhanced_image, cmap='gray')
	# # axs[1].set_title('Enhanced Image')
	# # axs[1].axis('off')
	# # plt.show()
	
	# #Image filtering
	# gaussian_filter = cv2.GaussianBlur(enhanced_image, (3, 3), 0)
	
	# #Binarization
	# threshold_value, binary_image = cv2.threshold(gaussian_filter, 50, 255, cv2.THRESH_BINARY)
	
	
	# #Delete objets with an area with less than 60 pixels
	# # perform connected component analysis
	# num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(binary_image)
	# # create an array to store the filtered image
	# filtered_image = np.zeros_like(binary_image)
	
	# # iterate over the connected components
	# for label in range(1, num_labels):
	    # # check the area of the component
	    # area = stats[label, cv2.CC_STAT_AREA]
	    # # filter out small objects
	    # if area >= 50:
	        # # add the component to the filtered image
	        # filtered_image[labels == label] = 255
	        
	# # fig, axs = plt.subplots(1, 2, figsize=(10, 5))
	# # axs[0].imshow(binary_image, cmap='gray')
	# # axs[0].set_title('thresh Image')
	# # axs[0].axis('off')
	# # axs[1].imshow(filtered_image, cmap='gray')
	# # axs[1].set_title('filtered Image')
	# # axs[1].axis('off')
	# # plt.show()
    
	# return enhanced_image,filtered_image
	

def houghlines():

	image='PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit'
	image = CCDData.read(f"{image}", unit = 'adu') 
	image = np.uint8(image.data/256) #convert uint16 to uint8
	
	enhanced_image,binary_image=binarization(image)
	
	# Perform edge detection
	# We need two threshold values, minVal and maxVal. Any edges with intensity gradient more than maxVal are sure to be edges 
	# and those below minVal are sure to be non-edges, so discarded.
	T_lower=100 #Lower threshold value in Hysteresis Thresholding
	T_upper=200	#Upper threshold value in Hysteresis Thresholding
	
	aperture_size=3 #Aperture size of the Sobel filter.
	#The default value is 3 and its value should be odd between 3 and 7. 
	#You can increase the Aperture size when you want to detect more detailed features.
	edges = cv2.Canny(binary_image, T_lower,T_upper, aperture_size)
	
	cv2.namedWindow('edge', cv2.WINDOW_NORMAL)  
	cv2.imshow('edge', edges)
	screen_width, screen_height = 1090, 1080  
	cv2.resizeWindow('edge', screen_width, screen_height)
	cv2.waitKey(0)
	cv2.destroyAllWindows()

	#Perform Probabilistic Hough Transform
	lines = cv2.HoughLinesP(edges,rho=1,theta=np.pi/180,threshold=20,minLineLength=100,maxLineGap=10)
	
	# Extract line segment endpoints from Hough lines
	line_segments = []
	if lines is not None:
		for line in lines:
			x1, y1, x2, y2 = line[0]
			line_segments.append((x1, y1, x2, y2))
	
	for line in line_segments:
			x1=line[0]
			y1=line[1]
			x2=line[2]
			y2=line[3]
			
			min_x = min(x1, x2)
			max_x = max(x1, x2)
			min_y = min(y1, y2)
			max_y = max(y1, y2)
	       
			
	print(f'min_x={min_x},max_x={max_x},min_y={min_y},max_y={max_y}')
	for segment in line_segments:
	    x1, y1, x2, y2 = segment
	    cv2.line(enhanced_image, (x1, y1), (x2, y2), (0, 2, 255), 2)  
	
	# cv2.circle(enhanced_image, (min_x, min_y), 5, (0, 255, 0), 2)
	# cv2.circle(enhanced_image, (max_x, max_y), 5, (0, 255, 0), 2)

	
	cv2.namedWindow('Probabilistic Hough Transform', cv2.WINDOW_NORMAL)  
	cv2.imshow('Probabilistic Hough Transform', enhanced_image)
	screen_width, screen_height = 1090, 1080  
	cv2.resizeWindow('Probabilistic Hough Transform', screen_width, screen_height)
	cv2.waitKey(0)
	cv2.destroyAllWindows()
	
	
#houghlines()
	



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
