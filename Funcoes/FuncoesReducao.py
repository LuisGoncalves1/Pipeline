import ccdproc 
from datetime import date
import os
from astropy.nddata import CCDData
import numpy as np
import astropy.units as u
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
import numpy as np
import logging
import matplotlib.pyplot as plt
import astroscrappy
from pathlib import Path
import cv2
import time
import subprocess
import itertools
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky
from .FuncoesMaster import BestDark



def sextractor(image,directory,catout):
	''' image -> image name
		directory -> Pipeline folder directory
		catout -> directory of output catalog
	
	
	'''
	SexDir = directory + '/Files/Sextractor/'
	cmd = {}
	cmd['source-extractor'] = image
	cmd['-c'] = SexDir + 'default.sex'
	# ~ cmd['-CATALOG_NAME'] = SexDir + catout
	cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items())))
	result = subprocess.run(list_cmd, capture_output=True)
	return result
	
def astrometry(DirImagem,header,DirAstro):
	cmd = ['solve-field', '--no-plots',
           # ~ '--x-column', 'X_IMAGE',
           # ~ '--y-column', 'Y_IMAGE',
           '--no-remove-lines',
           '--uniformize', '0',
           '--width', str(header['NAXIS1']),
           '--height', str(header['NAXIS2']),           
           '--crpix-center',
           # ~ '--depth', depth_str,
           '--cpulimit', '120',
           '--scale-low', '0.1',
           '--scale-high', '1',
           '--scale-units', 'app',
           '--ra', header['RA'],
           '--dec', header['DEC'],
           '--tweak-order','2',
           '--radius', '0.5',
           # ~ '--new-fits', 'none',
           '--overwrite',
           '--dir', DirAstro,
           # ~ '--out', 'Ola.fits',
            DirImagem
    ]
    # log cmd executed
	cmd_str = ' '.join(cmd)
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode
	return result, status
	      		
def ReadSextractorFile(SexFile):
	dic = {}
	config = {}
	with open(SexFile,'r') as text:
		while True:
			line = text.readline()
			A = line.split()
			if line != '':
				if line[0] == '#':
					dic[int(A[1])] = A[2]
					config[A[2]] = []
				else:
					for i in range(len(A)):
						config[dic[i+1]].append(A[i])
			else:
				break
			
	return config		
				
def str2bool(valor):
	if valor == "True":
		return True
	if valor == "False":
		return False
	else:
		return valor

def ReadPipelineConfig(ficheiro):
	dic = {}
	with open(ficheiro,'r') as texto:
		while True:
			linha = texto.readline()
			if len(linha)!= 0 and linha[0] != '#' and linha[0]!= '\n':
				B = linha.strip().split()
				if B[1] == 'True' or 'False':
					dic[B[0]] = str2bool(B[1])
				else:
					dic[B[0]] = int(B[1])
			if linha == "":
				break
	return dic

def AnalisarImagem(ficheiro,Masters,logging,config,diretoria):
	logging.debug("Startins to analise " + ficheiro)
	Data = CCDData.read(diretoria + "/ImagensTratadas/" + ficheiro, unit = 'adu') 
	raw = Data.data.copy()
	Filter = Data.header['FILTER']
	exposure = int(float(Data.header['EXPTIME'].replace(',','')))
	
	
	if config["DARK"]:
		logging.debug("Starting Dark correction of file " + ficheiro)
		Master = BestDark(diretoria + '/Masters',exposure,Filter)
		if Master == []:
			logging.warning('There are no Darks available of type %s' (Filter))
		else:
			Master = diretoria + '/Masters/' + Master
			MasterDark = CCDData.read(Master) 
			Data = ccdproc.ccd_process(Data, dark_frame=MasterDark,dark_exposure = float(MasterDark.header["EXPOSURE"].replace(',','.')) * u.s,data_exposure = float(Data.header["EXPOSURE"].replace(',','.')) *u.s)
			Data.header['DARK'] = True
			logging.debug("Finished Dark correction of file " + ficheiro)
							
	if config["FLAT"]:
		logging.debug("Starting Fits correction of file " + ficheiro)
		if Filter not in Masters['Flat']:
			logging.warning('There are no Flats available of type %s' (Filter))
		else:
			Master = diretoria + '/Masters/MasterFlat' + Filter + '.fits'
			MasterFlat = CCDData.read(Master)	
			Data = ccdproc.ccd_process(Data, master_flat = MasterFlat)
			Data.header['FLAT'] = True
			logging.debug("Finished Fits correction of file " + ficheiro)	
				
	if config["Background"]:
		logging.debug("Eliminating Background of " + ficheiro)
		bkg = MeanBackground(SigmaClip(sigma=3.))
		bkg_value = bkg.calc_background(Data)
		Data = Data.subtract(bkg_value*np.ones(np.shape(Data))*u.adu,propagate_uncertainties=True,handle_meta='first_found')
		logging.debug("Background removed from " + ficheiro)
					
	if config["CrossTalk"]:
		logging.debug("Starting CrossTalk correction of " + ficheiro)
		alfa1 = 0.015
		alfa2 = 0.008
		convolution_kernel = np.array([[0, alfa1, 0], 
									   [alfa1, 1-4*alfa1, alfa1], 
									   [0, alfa1, 0]])

		Data.data = cv2.filter2D(np.array(Data.data), -1, convolution_kernel)
		logging.debug("Crosstalk corrected from " + ficheiro)

	if config["CosmicRays"]:
		logging.debug("Starting CosmicRays Correction of " + ficheiro)
		crmask, Data.data = astroscrappy.detect_cosmics(Data.data, inmask = None, sepmed = False, cleantype='medmask', fsmode = 'median', verbose=False)  #deteta os raios cósmicos numa matriz numpy
		logging.debug("CosmicRays removed from " + ficheiro)
		
	
						
	Data.write(diretoria  + '/ImagensReduzidas/' + ficheiro)


	if config['GRAPH']:
		fig,axes = plt.subplots(1,2)
		ax = axes.ravel()
		vmin,vmax = None,None
		ax[0].set_title('Imagem Inicial')
		ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].set_title('Imagem Final')
		ax[1].imshow(Data.data,cmap='gray',vmin=vmin or np.nanpercentile(Data.data, 0.1),vmax=vmax or np.nanpercentile(Data.data, 99))
		plt.show() 

def SeeFitsTables(filedir):
	hdul = fits.open(filedir)
	return hdul[1].data

def ReadCatalogue(catalog_path,max_magnitude=None,min_magnitude=None,max_variability=None):
	catalog = Table.read(catalog_path)
	mask = np.isfinite(catalog['ra']) & np.isfinite(catalog['dec'])
	catalog = catalog[mask]

	if max_magnitude is not None:
		catalog = catalog[catalog['v_mag'] <= max_magnitude]
        
	if min_magnitude is not None:
		catalog = catalog[catalog['v_mag'] >= min_magnitude]

	if max_variability is not None:
		catalog = catalog[catalog['variability'] == max_variability]
	
	
	# ~ return SkyCoord(ra=Angle(catalog['ra']), dec=Angle(catalog['dec']))
	return SkyCoord(ra = Longitude(catalog['ra'], unit = u.deg), dec = np.array(catalog['dec'])*u.deg,frame = 'cirs'), catalog['v_mag']

def find_matching_stars(catalog_stars, image_stars, max_sep=0.5 * u.deg):
    if len(image_stars) == 0:
        return np.array([], dtype=int), np.zeros(len(catalog_stars), dtype=bool)

    idx, d2d, d3d = match_coordinates_sky(catalog_stars, image_stars,storekdtree = True)
    mask = d2d < max_sep

    return idx[mask], mask
