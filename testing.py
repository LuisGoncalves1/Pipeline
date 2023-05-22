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


#-----------------Probabilistic Hough Transform--------------

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


#---------------------------- Grafico ----------------------------------


TLE_file=f'TLE Catalog-2022417.json'

def convert_degrees_to_hms(degrees):
    total_seconds = degrees * 3600  # 1 degree = 3600 seconds (60 minutes * 60 seconds)
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = int(total_seconds % 60)
    hms=str(hours)+':'+str(minutes)+':'+str(seconds)
    return hms
    
def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    dms=str(round(mult*deg))+':'+str(round(mult*mnt))+':'+str(float("{:.1f}".format(mult*sec)))
    return dms


def TLE_Graph(config,DayDirectory,TLE_file,TLE_Directory):

	
	with open(f'{TLE_Directory}/{TLE_file}') as user_file:
		file_contents = json.load(user_file)
		
	ts = load.timescale()
	t = ts.utc(2023,2,28,20,20,12)
		
	Ra=[]
	Dec=[]	
	Ra_ticks=[]
	Dec_ticks=[]
	for i in file_contents:
			
			line1=i['tle_1']
			line2=i['tle_2']
			nome='n'
	
			satellite = EarthSatellite(line1, line2, nome, ts)

			bluffton = wgs84.latlon(float(config['Latitude']), float(config['Longitude']),float(config['Height']))
			
			difference = satellite - bluffton 
		
			topocentric = difference.at(t)
	
			ra, dec, distance = topocentric.radec()
			
			ra=math.degrees(ra.radians)
			hms=convert_degrees_to_hms(ra)
			Ra.append(ra)
			Ra_ticks.append(hms)
			
			dec=math.degrees(dec.radians)
			dms=decdeg2dms(dec)
			Dec.append(dec)
			Dec_ticks.append(dms)

	
	return Ra,Dec,Ra_ticks,Dec_ticks

import matplotlib.dates as mdates
from datetime import datetime
import matplotlib.ticker as ticker

Ra,Dec,Ra_ticks,Dec_ticks,=TLE_Graph(config,DayDirectory,TLE_file,TLE_Directory)
plt.scatter(Ra, Dec,alpha=0.5)

ax = plt.gca()
plt.xlabel('Ra (hms) ')
#plt.xticks(Ra, Ra_ticks)
plt.ylabel('Dec (dms)')
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=2))
plt.yticks(Dec, Dec_ticks)
plt.xlim([0, 359])
plt.ylim([-20, 20])
plt.show()


