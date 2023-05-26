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

#SPADEFunctions.Run_TLE_Code(config,PipelineDirectory,DayDirectory,FileXY,TLE_Directory)
	
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


#-----------------Probabilistic Hough Transform--------------------------------------------------------------

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

def houghlines():

	#import cv2
	image='PAER_SURV_4_260753_230228_042347+065500_scan017_01.fit'
	image = CCDData.read(f"{image}", unit = 'adu') 
	img_uint8 = np.uint8(image.data/256) #convert uint16 to uint8
	
	# Perform edge detection
	#We need two threshold values, minVal and maxVal. Any edges with intensity gradient more than maxVal are sure to be edges 
	#and those below minVal are sure to be non-edges, so discarded.
	T_lower=15  #Lower threshold value in Hysteresis Thresholding
	T_upper=17	#Upper threshold value in Hysteresis Thresholding
	
	aperture_size=7 #Aperture size of the Sobel filter.
	#The default value is 3 and its value should be odd between 3 and 7. 
	#You can increase the Aperture size when you want to detect more detailed features.
	edges = cv2.Canny(img_uint8, T_lower,T_upper, aperture_size)
	
	
	# cv2.imshow('Canny Edges', edges)
	# cv2.waitKey(0)
	# cv2.destroyAllWindows()
	
	
	
	#Perform Probabilistic Hough Transform
	
	lines = cv2.HoughLinesP(edges,rho=1,theta=np.pi/180,threshold=10,minLineLength=30,maxLineGap=30)
	
	# Extract line segment endpoints from Hough lines
	line_segments = []
	if lines is not None:
		for line in lines:
			x1, y1, x2, y2 = line[0]
			line_segments.append((x1, y1, x2, y2))
	#print(len(line_segments))
	
	for segment in line_segments:
	    x1, y1, x2, y2 = segment
	    cv2.line(img_uint8, (x1, y1), (x2, y2), (0, 0, 255), 2)  # Draw line segments on the image
	
	cv2.imshow('Probabilistic Hough Transform', img_uint8)
	cv2.waitKey(0)
	cv2.destroyAllWindows()


#---------------------------- Grafico --------------------------------------------------------------------


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
	ts = load.timescale()
	t_inicial=ts.utc(2023,2,28,19,59)
	t_final=ts.utc(2023,2,28,20,21,49) 

	Ra=[]
	Dec=[]	

	for i in file_contents:

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
	

	# ra_target=['04 23 47.00','04 23 47.00','04 23 47.00','04 23 47.00','04 23 47.00','08 07 00.00','08 07 00.00','08 07 00.00','08 07 00.00']
	# dec_target=['-3 13 45.0','+6 55 00.0','-0 41 33.0','+1 50 37.0','-8 18 07.0','-7 43 53.0','+4 57 03.0','-5 11 42.0','+2 24 52.0']
	# ra_degrees= HMS2deg(ra_target) 
	# dec_degrees= DMS2deg(dec_target)
	#plt.scatter(ra_degrees, dec_degrees,marker='s',s=40)
	
	
	files=os.listdir(f"{DayDirectory}/TLE")
	
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
	

	

	


	plt.xlabel('Ra (Hours) ',fontsize=12)
	plt.ylabel('Dec (Degrees)',fontsize=12)

	angle=[0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360]
	x=[]
	for i in angle:
		hms=convert_degrees_to_hours_minutes(i)
		x.append(hms)
	
	plt.xticks(angle, x)
	
	plt.xticks(angle, x)

	plt.xlim([0, 359])
	plt.ylim([-20, 20])
	plt.show()
TLE_Graph(config,DayDirectory,TLE_file,TLE_Directory)

WindowFile=['PAER_SURV_4_260789_230228_080700-051142_scan102_01_RaDec_Window','PAER_SURV_4_260789_230228_080700-051142_scan102_02_RaDec_Window',
'PAER_SURV_4_260789_230228_080700-051142_scan102_03_RaDec_Window','PAER_SURV_4_260789_230228_080700-051142_scan102_04_RaDec_Window']
