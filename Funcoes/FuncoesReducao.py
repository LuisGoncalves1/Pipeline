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
	''' Function to run sextractor on an image and output catalogue file
	image -> image directory to give to sextractor
	directory -> directory of the pipeline folder
	catout -> directory and name of the output catalog with all information chosen
	
	'''
	
	
	
	SexDir = directory + '/Files/Sextractor/' 						# The sextractor directory is inside the pipeline folder
	cmd = {}		
	cmd['source-extractor'] = image									# Chooses the image directory as the image to analyse
	cmd['-c'] = SexDir + 'default.sex'								# Chooses the sextractor catalogue as the default in the pipeline folder
	cmd['-CATALOG_NAME'] =  catout									# Chooses the catalogue output as the catout given in the function
	cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'				# Chooses the needed files from the sextractor folder
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items()))) # Compacts the cmd dictionary created so it can be run through terminal
	result = subprocess.run(list_cmd, capture_output=True) 			  # Runs the program in terminal and outputs the result and status (0 = good, 255 = bad)
	status = result.returncode
	return result
	    		
def astrometry(DirCatalogo,header,DirAstro,config):
	''' Function to run astrometry.net on an image and output created files
	DirCatalogo -> directory where the image to be analysed is, in this case its an sextractor catalogue with X,Y positions
	Header -> Image header, to know information such as width,height and center-pixel Ra and Dec
	DirAstro -> Directory where we are going to save the output files
	config -> Config file for the chosen telescope, with informations regarding the parameters to be used to run astrometry
	'''
	
	
	cmd = ['solve-field', '--no-plots',
	       '--x-column', 'X_IMAGE',					# Columns in sextractor 
           '--y-column', 'Y_IMAGE',					# Columns in sextractor
           '--no-remove-lines',
           # ~ '--uniformize', '0',
           '--width', str(header['NAXIS1']),    	# Width of image
           '--height', str(header['NAXIS2']),       # Height of image
           '--crpix-center',						# The reference pixel is in the center
           '--no-verify',
           '--cpulimit', config['cpulimit'],		# Max time for computation, stops after this time
           '--scale-low', config['scale-low'],		# Minimum scale to be used by astrometry
           '--scale-high', config['scale-high'],	# Maximum scale to be used by astrometry
           '--scale-units', config['scale-units'], 	# In what units are given each scale
           '--ra', header['RA'],					# Center-pixel RA coordinate in HH:MM:SS or degrees
           '--dec', header['DEC'],					# Center-pixel Dec coordinate in DD:MM:SS or degrees
           '--pixel-error','2',
           '--tweak-order','2',
           '--rdls', 'none',						# none in all files which we do not want, if you want
           '--solved','none',						# any of the files, comment the line
           '--corr','none',
           '--axy','none',
           '--match','none',
           # ~ '--new-fits', 'none',
           '--radius', config['radius'],			# Radius of the image in degrees (in the diagonal of the image)
           '--overwrite',							# Overwrite if there is already a file with the same name
           '--dir', DirAstro,						# Output directory
            DirCatalogo
    ]
    # log cmd executed
	cmd_str = ' '.join(cmd) 						# Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True) # Runs astrometry.net with given parameters
	status = result.returncode						  # Returns result, status = 0 if good and 255 if bad
	return result, status	      		
	      			      		      		
def XYToRaDec(FileWc,FileXY,FileOutput):
	''' Function created to transform XY positions of an image into RaDec coordinates
	using WCS created by previous astrometry run. IT uses an astrometry function wcs-xy2rd.
	
	FileWc -> directory of WCS file of astrometry for corresponding image
	FileXY -> directory of XY coordinates file given by sextractor
	FileOutput -> directory of output file to write with Ra Dec coordinates of corresponding XY coordinates
	'''
	
	
	cmd = ['wcs-xy2rd', 
		'-w',FileWc, 							# WCS file directory
		'-i', FileXY,							# Sextractor catalogue with XY positions file directory
		'-o', FileOutput,						# Output file directory
		'-X','X_IMAGE',
		'-Y','Y_IMAGE']
		
	cmd_str = ' '.join(cmd)						# Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode  						# Returns result, status = 0 if good and 255 if bad
	return result, status	           		
	      		
def ReadSextractorFileTxt(SexFile):
	''' Function created to read catalogue created from sextractor in txt type. It outputs 
	a dictionary with the information related to each column as values of keys
	SexFile -> Directory of sextractor fits catalogue
	
	'''
	
	
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
	
def ReadSextractorFileFit(SexFile):
	''' Function created to read catalogue created from sextractor in fits type. It outputs 
	a dictionary with the information related to each column as values of keys
	SexFile -> Directory of sextractor fits catalogue
	
	'''
	
	
	config = {}
	Ficheiro = fits.open(SexFile)
	Dados = Ficheiro[1]
	Colunas = Dados.columns
	for i in Colunas:
		Coluna = i.name
		config[Coluna] = Dados.data[Coluna]
		
		
	return config
							
def str2bool(valor):
	''' Function to transform a string boolean into a boolean, if possible. 
	valor -> string to be converted 
	'''
	if valor == "True":
		return True
	if valor == "False":
		return False
	else:
		return valor

def ReadPipelineConfig(ficheiro):
	''' Function to read pipeline config files. It skips over blank lines and lines which start with an #. 
	Returns in the form of a dictionary with properties as keys and True/False or numbers as values
	ficheiro -> Config file directory
	
	'''
	
	
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

def AnalisarImagem(ficheiro,Masters,logging,config,diretoria,directory,catalogo,mag):
	''' Main function to analise an image. 
	ficheiro -> Directory of image to analise
	Masters -> Dictionary if available masters
	logging -> Logger
	config -> Telescope config dictionary
	diretoria -> Images directory for that day
	directory -> Pipeline directory
	catalogo -> Catalogue to be used for astrometry in SkyCoord type
	mag -> Corrrespoinding magnitudes for catalogo positions
	
	Writes in header: DARK (True or False); 
	FLAT (True or False);
	BG (True or False);
	CT (True of False);
	CR (True or False);
	
	'''
	
	
	logging.debug("Startins to analise " + ficheiro)										# Logger information
	Data = CCDData.read(diretoria + "/ImagensTratadas/" + ficheiro, unit = 'adu') 			# Reads Image
	raw = Data.data.copy()																	# Copys data so it conserves initial image
	Filter = Data.header['FILTER']															# Gets filter from header
	exposure = int(float(Data.header['EXPTIME'].replace(',','')))
	
	A = time.time()
	if config["DARK"]:																		# If DARK is True does dark correction
		logging.debug("Starting Dark correction of file " + ficheiro)						# Logger information
		Master = BestDark(diretoria + '/Masters',exposure,Filter)							# Sees what is the best masterdark to apply
		if Master == []:
			logging.warning('There are no Darks available of type ' + Filter)				# IF there is none, it skips over dark correction
			Data.header['DARK'] = 'False'													# Writes in header
		else:
			Master = diretoria + '/Masters/' + Master										# Gets master directory
			MasterDark = CCDData.read(Master) 												# Reads and then applies master
			Data = ccdproc.ccd_process(Data, dark_frame=MasterDark,dark_exposure = float(MasterDark.header["EXPOSURE"].replace(',','.')) * u.s,data_exposure = float(Data.header["EXPOSURE"].replace(',','.')) *u.s)
			Data.header['DARK'] = 'True'													# Writes in header
			logging.debug("Finished Dark correction of file " + ficheiro)					# Logger information
	else:	
		Data.header['DARK'] = 'False'														# Writes in header
		
	A = time.time() - A
	
	B = time.time()
	if config["FLAT"]:
		logging.debug("Starting Fits correction of file " + ficheiro)						# Logger information
		if Filter not in Masters['Flat']:
			logging.warning('There are no Flats available of type ' + Filter)				# Logger information if there is no flat
			Data.header['FLAT'] = 'False'													# Writes in header
		else:
			Master = diretoria + '/Masters/MasterFlat' + Filter + '.fits'					# Gets master directory
			MasterFlat = CCDData.read(Master)												# Reads master
			Data = ccdproc.ccd_process(Data, master_flat = MasterFlat)						# Applies Master
			Data.header['FLAT'] = 'True'													# Writes in header
			logging.debug("Finished Fits correction of file " + ficheiro)					# Logger information
	else:
		Data.header['FLAT'] = 'False'														# Writes in header
	
	B = time.time()-B 													
				
	C = time.time()
	if config["Background"]:
		logging.debug("Eliminating Background of " + ficheiro)								# Logger information
		bkg = MeanBackground(SigmaClip(sigma=3.))
		bkg_value = bkg.calc_background(Data)												# Calculates background form image
		Data = Data.subtract(bkg_value*np.ones(np.shape(Data))*u.adu,propagate_uncertainties=True,handle_meta='first_found') # Removes background from image
		logging.debug("Background removed from " + ficheiro)								# Logger information
		Data.header['BG'] = 'True'															# Writes in header
	else:
		Data.header['BG'] = 'False'															# Writes in header
	
	C = time.time() - C											
	
	D = time.time()					
	if config["CrossTalk"]:																	# Stars crostalk correction if enabled
		logging.debug("Starting CrossTalk correction of " + ficheiro)						# Logger information
		alfa1 = config['alfa1']																# Parameters for correction
		alfa2 = config['alfa2']
		convolution_kernel = np.array([[0, alfa1, 0], 										# Deconvolution kernel for cross-talk correction
									   [alfa1, 1-4*alfa1, alfa1], 
									   [0, alfa1, 0]])

		Data.data = cv2.filter2D(np.array(Data.data), -1, convolution_kernel)				# Applies kernel to image
		logging.debug("Crosstalk corrected from " + ficheiro)
		Data.header['CT'] = 'True'															# Writes in header
	else:
		Data.header['CT'] = 'False'															# Writes in header
	D = time.time()-D
	
	E = time.time()
	if config["CosmicRays"]:
		logging.debug("Starting CosmicRays Correction of " + ficheiro)						# Logger information
		crmask, Data.data = astroscrappy.detect_cosmics(Data.data, inmask = None, sepmed = False, cleantype='medmask', fsmode = 'median', verbose=False)  #deteta os raios cósmicos numa matriz numpy
		logging.debug("CosmicRays removed from " + ficheiro)								# Logger information
		Data.header['CR'] = 'True'															# Writes in header
	else:
		Data.header['CR'] = 'False'														# Writes in header
	E = time.time() - E
	
	DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		
	Data.write(DirImage)
	
	if config["Astrometry"]:
		#Runs sextractor on reduced image and produces sextractor catalogue
		logging.debug("Starting Sextractor on " + ficheiro)
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'
		SexTime = time.time()
		resultado = sextractor(DirImage,directory,SexFile)     
		SexTime = time.time() - SexTime
		
		#Runs astrometry with sextractor catalogue and outputs in astrometry folder (corr,wcs,xyls) files
		logging.debug("Starting Astrometry on " + ficheiro)
		AstTime = time.time()
		result,status = astrometry(diretoria + '/SextractorCatalogue/' + ficheiro[:-5] + '.cat',Data.header,diretoria + '/AstrometryFiles',config) 
		AstTime = time.time() - AstTime
		
		#Transforms from XY to RaDec with WCS information from astrometry
		nomeficheiro = ficheiro.split('/')[-1].split('.')[0]							# Saves name of file (removes rest of extension)
		FileWc = diretoria + '/AstrometryDefesa/' + nomeficheiro + '.wcs'				# Gets directory : WCS file created from astrometry 
		FileXY = diretoria + '/SextractorCatalogue/' + nomeficheiro + '.cat'													# Gets directory : catalogue created from sextractor
		FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro														# Output file for calculation of RaDec coordinate of all stars
		TransfTime = time.time()
		result, status	= XYToRaDec(FileWc,FileXY,FileOutput)			# Calculates RaDec positions for all stars obtained from sextractor
		TransfTime = time.time() - TransfTime

		#Gets the information from the sextractor file and the calculated positions
		Posicoes = fits.open(FileOutput)												# Gets positions calculated
		Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']							# Saves RaDec positions
		SexData = FuncoesReducao.ReadSextractorFileFit(SexFile)			# Reads sextractor file obtained in the beggining and saves desired values
		X,Y,Mags,MagsErro,BackG,fluxo,FWHM,thresh= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE'],SexData['THRESHOLD']
		
		#Matches the objects found and the catalogue
		stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
		logging.debug("Finding match between catalog stars and detected objects in image " + ficheiro)
		MatchTime = time.time()
		idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=0.06 * u.deg) 						# Starts matching between these stars and catalogue stars
		MatchTime = time.time() - MatchTime
		EncontradosX,EncontradosY,EncontradosMags = X[found],Y[found],Mags[found] 									#Saves X,Y and Mags positions found
		logging.debug("found " + len(EncontradosX) + " catalog objects in image " + ficheiro)
		



	if config['GRAPH']:
		fig,axes = plt.subplots(1,2)
		ax = axes.ravel()
		vmin,vmax = None,None
		ax[0].set_title('Imagem Inicial')
		ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].set_title('Imagem Final')
		ax[1].imshow(Data.data,cmap='gray',vmin=vmin or np.nanpercentile(Data.data, 0.1),vmax=vmax or np.nanpercentile(Data.data, 99))
		plt.show() 
	
	logging.info("Tempos de :" + ficheiro)
	logging.info('TempoDark: {:.2f}'.format(A))
	logging.info('TempoFlat: {:.2f}'.format(B))
	logging.info('TempoBackGround: {:.2f}'.format(C))
	logging.info('TempoCrossTalk: {:.2f}'.format(D))
	logging.info('TempoCosmicRays: {:.2f}'.format(E))
	logging.info('TempoSextractor: {:.2f}'.format(SexTime))
	# ~ logging.info('TempoAstrometry: {:.2f}'.format(AstTime))
	# ~ logging.info('TempoTransformar: {:.2f}'.format(TransfTime))
	# ~ logging.info('TempoMatch: {:.2f}'.format(MatchTime))

	
def SeeFitsTables(filedir):
	hdul = fits.open(filedir)
	return hdul[1].data

def ReadCatalogue(catalog_path,max_magnitude=None,min_magnitude=None,max_variability=None):
	''' Function to read hipparcos catalogue (old)
	catalog_path -> Directory of catalog file
	max_mangitude -> maximum magnitude of object to save
	min_magnitude -> minimum magnitude of object to save
	max_variability -> maximum variability of object to save
	
	returns a SkyCoord object with Ra and Dec catalogue and an array with corresponding magnitude
	'''
	
	
	catalog = Table.read(catalog_path)
	mask = np.isfinite(catalog['ra']) & np.isfinite(catalog['dec'])
	catalog = catalog[mask]

	if max_magnitude is not None:
		catalog = catalog[catalog['v_mag'] <= max_magnitude]
        
	if min_magnitude is not None:
		catalog = catalog[catalog['v_mag'] >= min_magnitude]

	if max_variability is not None:
		catalog = catalog[catalog['variability'] == max_variability]
	
	print(type(catalog['ra']))
	# ~ return SkyCoord(ra=Angle(catalog['ra']), dec=Angle(catalog['dec']))
	return SkyCoord(ra = Longitude(catalog['ra'], unit = u.deg), dec = np.array(catalog['dec'])*u.deg,frame = 'cirs'), catalog['v_mag']

def ReadDucati(catalog_path,max_magnitude=None,min_magnitude=None):
	''' Function to read Duccati catalog from Vizier: https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/237/colors&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
	catalog_path -> Directory of catalog file
	max_mangitude -> maximum magnitude of object to save
	min_magnitude -> minimum magnitude of object to save
	
	
	
	returns a SkyCoord object with Ra and Dec catalogue and an array with corresponding magnitude
	'''
	
	catalog = {}
	HDU = fits.open(catalog_path)
	Data = HDU[1].data
	Names = Data.names
	for i in Names:
		if i == 'RAJ1950':
			catalog[i] = HMS2deg(Data[i])
		elif i == 'DEJ1950':
			catalog[i] = DMS2deg(Data[i])
		else:
			catalog[i] = np.array(Data[i])

	return SkyCoord(ra = Longitude(catalog['RAJ1950'], unit = u.deg), dec = np.array(catalog['DEJ1950'])*u.deg,frame = 'cirs'), catalog

def ReadTychoCatalogue(catalog_path,max_magnitude=None,min_magnitude=None):
	''' Function to read Tycho catalogue in txt form taken from http://archive.eso.org/ASTROM/TYC-2/data/
	catalog_path -> Directory of catalog file
	max_mangitude -> maximum magnitude of object to save
	min_magnitude -> minimum magnitude of object to save
	
	
	
	returns a SkyCoord object with Ra and Dec catalogue and an array with corresponding magnitude
	'''
	
	catalog = {'RA':[],'DEC':[],'BTMag':[],'VTMag':[]}
	with open(catalog_path,'r') as text:
		while True:
			line = text.readline()
			if line == '':
				break
			else:
				BTMag = line[110:116].strip()
				VTMag = line[123:129].strip()
				
				if BTMag != '' and VTMag != '':
					catalog['RA'].append(float(line[152:164].strip()))
					catalog['DEC'].append(float(line[165:177].strip()))
					catalog['BTMag'].append(float(BTMag))
					catalog['VTMag'].append(float(VTMag))
	catalog['BTMag'] = np.array(catalog['BTMag'])
	catalog['VTMag'] = np.array(catalog['VTMag'])
	return SkyCoord(ra = Longitude(catalog['RA'], unit = u.deg), dec = np.array(catalog['DEC'])*u.deg,frame = 'cirs'), catalog['VTMag'], catalog['BTMag']
	
def HMS2deg(RA):
	'''
	Input data: Ra coordinate array with values in strings of +HH:MM:SS or -HH MM SS form
	Returns decimal degrees correspondent
	'''
	
	degrees = []
	for value in RA:
		H = float(value[1:3])
		M = float(value[4:6])
		S = float(value[7:10])
		degree = 360/24 * H + 360/(24*60) * M + 360/(24*60*60) * S
		degrees.append(degree)
		
	return degrees

def DMS2deg(DEC):
	'''
	Input data: Dec coordinate array with values in strings of +DD:MM:SS or -DD MM SS form
	Returns decimal degrees correspondent
	'''
	
	degrees = []
	for value in DEC:
		D = float(value[1:3])
		M = float(value[4:6])
		S = float(value[7:10])
		degree = D + M/60 + S/3600
		degrees.append(degree)
		
	return degrees

def find_matching_stars(catalog_stars, image_stars, max_sep=0.5 * u.deg):
	''' Function to find the closest match between two sets of positions. 
	This is used to find the closest match between found objects and stars in catalogue.
	catalog_stars -> SkyCoord object with catalogue stars coordinates
	image_stars -> SkyCoord object with coordinates of objects detected in image 
	max_sep -> maximum separation to be accepted for closest match (otherwise it is not found)
	
	Returns 
	idx: index of image stars that were found
	mask: True/False matrix matching the image_stars array
	'''
	
	
	if len(image_stars) == 0:
		return np.array([], dtype=int), np.zeros(len(catalog_stars), dtype=bool)

	idx, d2d, d3d = match_coordinates_sky(catalog_stars, image_stars,storekdtree = True)
	mask = d2d < max_sep

	return idx[mask], mask
