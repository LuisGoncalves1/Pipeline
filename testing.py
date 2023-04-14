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

SPADEFunctions.Run_TLE_Code(config,PipelineDirectory,DayDirectory,FileXY,TLE_Directory)
	
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


			
		



