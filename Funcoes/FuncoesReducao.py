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



def sextractor(image,directory,catin,catout,cattype,VIGNET=None):
	''' Function to run sextractor on an image and output catalogue file
	INPUTS:
	image -> image directory to give to sextractor 														(str)
	directory -> directory of the Pipeline folder														(str)
	catin -> Configuration catalog used for Sextractor													(str)
	catout -> directory and name of the output catalog													(str)
	cattype -> Type of output catalog, FITS_1.0 for normal, and FITS_LDAC for PSFEX						(str)
	VINGET -> Size of array to save from each image, for PSFEX											(int)
	
	OUTPUTS:
	Writes a sextractor catalog at the day Sextractor folder with the same name as the original image
	
	'''
	
	
	
	SexDir = directory + '/Files/Sextractor/' 						# The sextractor directory is inside the pipeline folder
	cmd = {}		
	cmd['source-extractor'] = image									# Chooses the image directory as the image to analyse
	cmd['-c'] = SexDir + catin										# Chooses the sextractor catalogue as the default in the pipeline folder
	cmd['-CATALOG_NAME'] =  catout									# Chooses the catalogue output as the catout given in the function
	# ~ cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'				# Chooses the needed files from the sextractor folder
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	cmd['-CATALOG_TYPE'] = cattype
	
	if cattype == 'FITS_LDAC':
		cmd['-PARAMETERS_NAME'] = SexDir + 'default1.param'				# Chooses the needed files from the sextractor folder
	else:
		cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'				# Chooses the needed files from the sextractor folder
		
		
		
	# ~ if cattype == 'FITS_LDAC':
		# ~ with open(SexDir + 'default.param', 'a') as myfile:
			# ~ myfile.write('VIGNET{}\n'.format((VIGNET,VIGNET)))
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items()))) # Compacts the cmd dictionary created so it can be run through terminal
	result = subprocess.run(list_cmd, capture_output=True) 			  # Runs the program in terminal and outputs the result and status (0 = good, 255 = bad)
	status = result.returncode
	
	# ~ if cattype == 'FITS_LDAC':
		# ~ with open(SexDir + 'default.param', '+r') as myfile:
			# ~ lines = myfile.readlines()
			# ~ myfile.seek(0)
			# ~ myfile.truncate()
			# ~ myfile.writelines(lines[:-1])
			
	return result
	    		
def astrometry(DirCatalogo,header,DirAstro,config):
	''' Function to run astrometry.net on an image and output created files
	INPUTS:
	DirCatalogo -> directory of the sextractor catalogue with X,Y positions from the desired image								(str)
	Header -> Image header, to know information such as width, height and center-pixel Ra and Dec								(str)
	DirAstro -> Directory where we are going to save the output files															(str)
	config -> Config file for the chosen telescope, with informations regarding the parameters to be used to run astrometry		(str)
	
	OUTPUTS:
	Writes a .wcs and a .xyls file on astrometry, with the same name as the original image
	'''
	
	
	cmd = ['solve-field', '--no-plots',
	       '--x-column', 'X_IMAGE',					# Columns in sextractor 
           '--y-column', 'Y_IMAGE',					# Columns in sextractor
           '--no-remove-lines',
           '--uniformize', '0',
           '--width', str(header['NAXIS1']),    	# Width of image
           '--height', str(header['NAXIS2']),       # Height of image
           '--crpix-center',						# The reference pixel is in the center
           '--no-verify',
           '--cpulimit', config['cpulimit'],		# Max time for computation, stops after this time
           '--scale-low', config['scale-low'],		# Minimum scale to be used by astrometry
           '--scale-high', config['scale-high'],	# Maximum scale to be used by astrometry
           '--scale-units', config['scale-units'], 	# In what units are given each scale
           '--ra', header['RA'].replace(' ',':'),					# Center-pixel RA coordinate in HH:MM:SS or degrees
           '--dec', header['DEC'].replace(' ',':'),					# Center-pixel Dec coordinate in DD:MM:SS or degrees
           '--pixel-error','2',
           '--tweak-order','2',
           # ~ '--rdls', 'none',						# none in all files which we do not want, if you want
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
	cmd_str = ' '.join(cmd) 						  # Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True) # Runs astrometry.net with given parameters
	status = result.returncode						  # Returns result, status = 0 if good and 255 if bad
	return result, status	      		
	
def psfex(catin,file_config):
	cmd = ['psfex', catin, '-c', file_config] 
	cmd_str = ' '.join(cmd) 						  # Compacts the cmd dictionary created so it can be run through terminal
	print(cmd_str)
	
	result = subprocess.run(cmd, capture_output=True) # Runs psfex with given parameters
	print(result)
	status = result.returncode						  # Returns result, status = 0 if good and 255 if bad
	return result, status	
      			      		      		
def XYToRaDec(FileWc,FileXY,FileOutput):
	''' Function created to transform XY positions of an image into RaDec coordinates
	using WCS created by previous astrometry run. IT uses an astrometry function wcs-xy2rd.
	
	FileWc -> directory of WCS file of astrometry for corresponding image										(str)
	FileXY -> directory of XY coordinates file given by sextractor												(str)
	FileOutput -> directory of output file to write with Ra Dec coordinates of corresponding XY coordinates		(str)
	'''
	
	
	cmd = ['wcs-xy2rd', 
		'-w',FileWc, 									# WCS file directory
		'-i', FileXY,									# Sextractor catalogue with XY positions file directory
		'-o', FileOutput,								# Output file directory
		'-X','X_IMAGE',
		'-Y','Y_IMAGE']
		
	cmd_str = ' '.join(cmd)								# Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode  						# Returns result, status = 0 if good and 255 if bad
	return result, status	           		
	      
def RaDecToXY(FileWc,FileRaDec,FileOutput):
	cmd = ['wcs-rd2xy', 
		'-w',FileWc, 									# WCS file directory
		'-i', FileRaDec,									# Sextractor catalogue with XY positions file directory
		'-o', FileOutput]								# Output file directory
		
	cmd_str = ' '.join(cmd)								# Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode  						# Returns result, status = 0 if good and 255 if bad
	return result, status

def CreateRDLSFile(Ra,Dec,FileName):
	''' Function that changes the data from the RLDS file to a set of (Ra,Dec) coordinates.
	INPUTES:
	Ra -> array with Ra coordinates
	Dec -> array with Dec coordinates
	FileName -> Name of rdls file
	
	'''
	with fits.open(FileName, mode = 'update') as A:
		Dados = A[1].data
		
		x = np.array(Ra,Dec,dtype = [('RA', '<f8'),('DEC','<f8')])
		Dados = Dados[Dados['RA'] > 180]
		Dados = np.append(Dados,x)

		A.close()
		      		
def ReadSextractorFileTxt(SexFile):
	''' Function created to read catalogue created from sextractor in txt type. It outputs 
	a dictionary with the information related to each column as values of keys
	SexFile -> Directory of sextractor fits catalogue							(str)
	
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
	SexFile -> Directory of sextractor fits catalogue							(str)
	
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
	valor -> string to be converted 		(str)
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
	ficheiro -> Config file directory							(str)
	
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

def BestDark(Directory,time):
	''' Checks to see whats the best dark to use for a given image, checks filter
	and closest exposure time
	Directory -> MasterDarks directory 					(str)
	time -> Exposure time for image 					(float)
	'''
	files = os.listdir(Directory)
	files = [B for B in files if B[6:10] == 'Dark' and int(B[10:13]) <= time]
	# ~ files = [B for B in files if B[6:10] == 'Dark' and B[10] == Filter and int(B[11:14]) <= time]
	times = [B[11:14] for B in files]
	if times == []:
		return []
	else:
		return files[times.index(max(times))]

def AnalisarImagem(ficheiro,Masters,logging,config,diretoria,directory,catalogo,mag):
	''' Main function to analise an image. 
	Image reduction -> Sextractor twice -> PSFEX? -> Astrometry -> Comparison with catalog -> End
	ficheiro -> Directory of image to analise									(str)
	Masters -> Dictionary if available masters									(dict)
	logging -> Logger															(log)
	config -> Telescope config dictionary										(str)
	diretoria -> Images directory for that day									(str)
	directory -> Pipeline directory												(str)
	catalogo -> Catalogue to be used for astrometry in SkyCoord type			(numpy.list)
	mag -> Corrrespoinding magnitudes for catalogo positions					(numpy.list)
	
	Needs in header:
	EXPTIME/EXPOSURE
	FILTER
	ǴAIN
	
	
	
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
	# ~ exposure = float(Data.header['EXPTIME'].replace(',','.'))							# Bruno
	exposure = float(Data.header['EXPTIME'])												# Defesa
	
	A = time.time()
	if config["DARK"]:																		# If DARK is True does dark correction
		logging.debug("Starting Dark correction of file " + ficheiro)						# Logger information
		Master = BestDark(diretoria + '/Masters',exposure)									# Sees what is the best masterdark to apply
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
		
	if config['GAIN']:
		Data = ccdproc.gain_correct(Data, float(Data.header['GAIN']))#*u.electron/u.adu
		Data.header['GAIN'] = '1'
		
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
		alfa1 = float(config['alfa1'])																# Parameters for correction
		alfa2 = float(config['alfa2'])
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
		Data.header['CR'] = 'False'															# Writes in header
	E = time.time() - E
	
	DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro	
	Data.write(DirImage, overwrite = True)

	if config["Astrometry"]:
		#Runs sextractor on reduced image and produces sextractor catalogue
		logging.debug("Starting Sextractor on " + ficheiro)
		SexTime = time.time()
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
		resultado = sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_LDAC',int(config['VIGNET']))  
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'
		resultado = sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_1.0')  
		RemoveBorders(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat',2,int(config['VIGNET']))
		RemoveBorders(diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat',1,int(config['VIGNET']))	
		SexTime = time.time() - SexTime
		logging.debug("Finished Sextractor on " + ficheiro)

		
		
		# ~ #Runs astrometry with sextractor catalogue and outputs in astrometry folder (corr,wcs,xyls) files
		logging.debug("Starting Astrometry on " + ficheiro)
		AstTime = time.time()
		result,status = astrometry(diretoria + '/SextractorCatalogue/' + ficheiro[:-5] + '.cat',Data.header,diretoria + '/AstrometryFiles',config) 
		AstTime = time.time() - AstTime
		logging.debug("Finished Astrometry on " + ficheiro)
		
		# ~ #Transforms from XY to RaDec with WCS information from astrometry
		nomeficheiro = ficheiro.split('/')[-1].split('.')[0]														# Saves name of file (removes rest of extension)
		FileWc = diretoria + '/AstrometryFiles/' + nomeficheiro + '.wcs'											# Gets directory : WCS file created from astrometry 
		FileXY = diretoria + '/SextractorCatalogue/' + nomeficheiro + '.cat'										# Gets directory : catalogue created from sextractor
		FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro													# Output file for calculation of RaDec coordinate of all stars
		TransfTime = time.time()
		result, status	= XYToRaDec(FileWc,FileXY,FileOutput)														# Calculates RaDec positions for all stars obtained from sextractor
		TransfTime = time.time() - TransfTime

		#Gets the information from the sextractor file and the calculated positions
		Posicoes = fits.open(FileOutput)												# Gets positions calculated
		Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'							# Saves RaDec positions
		SexData = ReadSextractorFileFit(SexFile)			# Reads sextractor file obtained in the beggining and saves desired values
		X,Y,Mags,MagsErro,BackG,fluxo,FWHM= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE']
		
		#Matches the objects found and the catalogue
		stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
		MatchTime = time.time()
		WCS = fits.open(FileWc)
		scale = float(WCS[0].header['COMMENT'][-8].split()[1])
		
		max_sep = float(config['max_sep'])*np.median(FWHM)*scale/3600
		idx,found = find_matching_stars(stars, catalogo, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars
		MatchTime = time.time() - MatchTime
		EncontradosX = X[found]
		logging.debug("found " + len(EncontradosX) + " catalog objects in image " + ficheiro)
		



	if config['GRAPH']:
		fig,axes = plt.subplots(2,2)
		ax = axes.ravel()
		vmin,vmax = None,None
		ax[0].set_title('Imagem Inicial')
		ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].set_title('Imagem Final')
		ax[1].imshow(Data.data,cmap='gray',vmin=vmin or np.nanpercentile(Data.data, 0.1),vmax=vmax or np.nanpercentile(Data.data, 99))
		ax[2].set_title('Sextractor Detections')
		ax[2].imshow(Data.data,cmap='gray',vmin=vmin or np.nanpercentile(Data.data, 0.1),vmax=vmax or np.nanpercentile(Data.data, 99))
		
		ax[3].set_title('Sextractor')
		ax[3].imshow(Data.data,cmap='gray',vmin=vmin or np.nanpercentile(Data.data, 0.1),vmax=vmax or np.nanpercentile(Data.data, 99))
		plt.show() 
	
	logging.info("Tempos de :" + ficheiro)
	logging.info('TempoDark: {:.2f}'.format(A))
	logging.info('TempoFlat: {:.2f}'.format(B))
	logging.info('TempoBackGround: {:.2f}'.format(C))
	logging.info('TempoCrossTalk: {:.2f}'.format(D))
	logging.info('TempoCosmicRays: {:.2f}'.format(E))
	logging.info('TempoSextractor: {:.2f}'.format(SexTime))
	logging.info('TempoAstrometry: {:.2f}'.format(AstTime))
	logging.info('TempoTransformar: {:.2f}'.format(TransfTime))
	logging.info('TempoMatch: {:.2f}'.format(MatchTime))
	
	print(status)
	print(ficheiro)
	return None
	
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

	return idx, mask, d2d

def CreateMasters(DayDirectory,config):
	''' Function to create master files from available Flats and Darks. It saves with exposure
	and type of filter in name.
	DayDirectory -> Directory for folder with darks and flats files
	config -> Dictionary with configurations from Config File
	
	'''
	Darks = OrganizeDarks(DayDirectory) # Sees available Darks
	Flats = os.listdir(DayDirectory + '/Flats') # Sees available Flats
	Masters = SeeMasters(DayDirectory)	

	for Time in Darks:
		if Time not in Masters['Dark']:
			DFiles = Darks[Time]
			DarkFiles = []
			for i in DFiles:
				DarkFiles.append(DayDirectory + '/Darks/' + i)
			CreateMaster(DarkFiles,DayDirectory,Time,'Dark',config)	

	for Filter in Flats:
		if Filter not in Masters['Flat']:
			FFiles = os.listdir(DayDirectory + '/Flats/' + Filter)
			FlatFiles =[]
			for i in FFiles:
				FlatFiles.append(DayDirectory + '/Flats/' + Filter + '/' + i)
			CreateMaster(FlatFiles,DayDirectory,Filter,'Flat',config)
							
def CreateMaster(Files,DayDirectory,Filter,Type,config):
	''' Function to create single master file from set of files.
	Files -> Array with directory for files to be used 
	DayDirectory -> Directory for folder with darks and flats files
	Filter -> Filter used in files, to be used to save in name (B,V,etc... or Time in dark case)
	Type -> Type of master file (Dark,Flat)
	config -> Dictionary with configurations from Config File
	
	'''
	FilesRead = [CCDData.read(X, unit='adu', Format = 'fits') for X in Files]
	Master = ccdproc.combine(FilesRead,method='median')
	Master.header['FILTER'] = Filter
	if (Type == 'Flat') & config['FLATdark']:
		Dark = BestDark(DayDirectory + '/Masters',float(Master.header['EXPTIME'].replace(',','.')))							# Sees what is the best masterdark to apply
		if Dark == []:
			Master.header['DARK'] = 'False'												# Writes in header
		else:
			Dark = DayDirectory + '/Masters/' + Dark										# Gets Dark directory
			MasterDark = CCDData.read(Dark) 												# Reads and then applies master
			Master = ccdproc.ccd_process(Master, dark_frame=MasterDark,dark_exposure = float(MasterDark.header["EXPOSURE"].replace(',','.')) * u.s,data_exposure = float(MasterDark.header["EXPOSURE"].replace(',','.')) *u.s)
			Master.header['DARK'] = 'True'													# Writes in header
		
	MasterName = "Master%s%s.fits" % (Type, Filter)
	
	Master.write(DayDirectory + '/Masters/' + MasterName)
	return Master
		
def OrganizeDarks(DayDirectory):
	''' Returns DarkFiles available in a dictionary of the type:
	Dic = {100:[Dark100_001.fits,Dark100_002.fits],200:[Dark200_001.fits,Dark200_002.fits]}
	
	DayDirectory -> Directory for folder with darks and flats files
	'''
	Darks = os.listdir(DayDirectory + '/Darks')
	times = {}
	for i in Darks:
		time = i[4:7]
		if time in times:
			times[time].append(i)
		else:
			times[time] = [i]
	
	return times
			
def SeeMasters(DayDirectory):
	''' Returns a dictionary of the type:
	Dic = {Flat: ['A','B','C']; Dark:['50','70','100','200']}
	
	DayDirectory -> Directory for folder with darks and flats files
	
	'''
	MastersFiles = os.listdir(DayDirectory + '/Masters')
	Masters = {'Flat': [],'Dark':[]}
	for i in MastersFiles:
		Type = i[6:10]
		if Type =='Flat':
			Masters[Type].append(i[10])
		elif Type =='Dark':
			Masters[Type].append(i[10:13])
						
	return Masters
	
def RemoveBorders(File,n,P):
	with fits.open(File, mode = 'update') as A:
		Dados = A[n].data
		A[n].data = Dados[((Dados['Y_IMAGE']>P) & (Dados['Y_IMAGE']<2048-P)) & ((Dados['X_IMAGE']>P) & (Dados['X_IMAGE']<2048-P))]
		A.close()
