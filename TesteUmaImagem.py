from Funcoes import FuncoesReducao, FuncoesMaster
import ccdproc 
from datetime import date
import os
from astropy.nddata import CCDData
import astropy.units as u
from photutils import make_source_mask, Background2D, MeanBackground
from astropy.stats import SigmaClip
import numpy as np
import logging
import cv2
import time
import matplotlib.pyplot as plt
from pathlib import Path
import astroscrappy
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky


#diretoria -> Ficheiro do dia e directory -> Pasta inicial Pipeline


directory = os.getcwd()
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationBruno.txt")
config['LOGGER'] = int(config['LOGGER'])
# ~ data = date.today()
data = "2022-03-20"

catalog_path = directory + '/Files/Catalogos/tycho2.dat'
diretoria = "%s/Files/Imagens/%s" %(directory,str(data)) 
if not os.path.isdir(diretoria):
	print('Não existem dados ainda deste dia')
	os.mkdir(diretoria)
	os.mkdir(diretoria + "/Darks")
	os.mkdir(diretoria + "/Flats")
	os.mkdir(diretoria + "/Imagens") 
	os.mkdir(diretoria + "/ImagensTratadas")
	os.mkdir(diretoria + "/ImagensReduzidas")
	os.mkdir(diretoria + "/Masters")
	os.mkdir(diretoria + "/SextractorCatalogue")
	os.mkdir(diretoria + "/AstrometryFiles")

if type(config['LOGGER']) != int or 0>config['LOGGER'] or 50<config['LOGGER']:
	print('Logger fora de range')
LOG_FORMAT = "%(levelname)s %(asctime)s : %(message)s"
log_file_name = "%s/%s.log" %(diretoria,str(data))
logging.basicConfig(filename = log_file_name, level = config['LOGGER'],format = LOG_FORMAT, filemode ='w')
logger = logging.getLogger()
logging.debug("Starting create Master Files")
A = time.time()
FuncoesMaster.CreateMasters(diretoria,config)
logging.info("Tempos para criarMasters : {:.2f}".format(time.time()-A) )
logging.debug("Master Files Created")
Masters = FuncoesMaster.SeeMasters(diretoria)

logging.debug("Starting to read catalogue")
A = time.time()
catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)
logging.info("Tempos para ler catalogo : {:.2f}".format(time.time()-A) )
logging.debug("Catalogue Read")


ficheiro = '1037-293B_001.fits'
# ~ ficheiro = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fits'

logging.debug("Startins to analise " + ficheiro)										# Logger information
Data = CCDData.read(diretoria + "/ImagensTratadas/" + ficheiro, unit = 'adu') 			# Reads Image
raw = Data.data.copy()																	# Copys data so it conserves initial image
Filter = Data.header['FILTER']															# Gets filter from header
exposure = int(float(Data.header['EXPTIME'].replace(',','.')))							# Imagens Bruno
# ~ exposure = int(Data.header['EXPTIME'])												# Imagens Defesa

A = time.time()

if config["DARK"]:																		# If DARK is True does dark correction
	logging.debug("Starting Dark correction of file " + ficheiro)						# Logger information
	Master = FuncoesMaster.BestDark(diretoria + '/Masters',exposure)							# Sees what is the best masterdark to apply
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
	Data.header['CR'] = 'False'														# Writes in header
E = time.time() - E

DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		
Data.write(DirImage)

if config["Astrometry"]:
		#Runs sextractor on reduced image and produces sextractor catalogue
		logging.debug("Starting Sextractor on " + ficheiro)
		SexTime = time.time()
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
		resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_LDAC')  
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'
		resultado = FuncoesReducao.sextractor(DirImage,directory,config['SexConfig'],SexFile,'FITS_1.0',int(config['VIGNET']))  
		SexTime = time.time() - SexTime
		
		# ~ #Runs psfex
		logging.debug("Starting PSFex on " + ficheiro)
		PsfTime = time.time()
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'LDAC.cat'
		PsfFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'psfex.cat'
		PsfConfig = directory + '/Files/Sextractor/' + config['PSFexConfig']
		resultado = FuncoesReducao.psfex(SexFile,PsfConfig)  
		PsfTime = time.time() - PsfTime
		
		#Runs astrometry with sextractor catalogue and outputs in astrometry folder (corr,wcs,xyls) files
		logging.debug("Starting Astrometry on " + ficheiro)
		AstTime = time.time()
		result,status = FuncoesReducao.astrometry(diretoria + '/SextractorCatalogue/' + ficheiro[:-5] + '.cat',Data.header,diretoria + '/AstrometryFiles',config) 
		AstTime = time.time() - AstTime
		
		#Transforms from XY to RaDec with WCS information from astrometry
		nomeficheiro = ficheiro.split('/')[-1].split('.')[0]							# Saves name of file (removes rest of extension)
		FileWc = diretoria + '/AstrometryFiles/' + nomeficheiro + '.wcs'				# Gets directory : WCS file created from astrometry 
		FileXY = diretoria + '/SextractorCatalogue/' + nomeficheiro + '.cat'													# Gets directory : catalogue created from sextractor
		FileOutput = diretoria + '/AstrometryFiles/' + nomeficheiro				# Output file for calculation of RaDec coordinate of all stars														
		TransfTime = time.time()
		result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)			# Calculates RaDec positions for all stars obtained from sextractor
		TransfTime = time.time() - TransfTime
		
		#Gets the information from the sextractor file and the calculated positions
		Posicoes = fits.open(FileOutput)												# Gets positions calculated
		Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
		SexFile = diretoria + '/SextractorCatalogue/'+ ficheiro[:-5] +'.cat'							# Saves RaDec positions
		SexData = FuncoesReducao.ReadSextractorFileFit(SexFile)			# Reads sextractor file obtained in the beggining and saves desired values
		X,Y,Mags,MagsErro,BackG,fluxo,FWHM,thresh= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE'],SexData['THRESHOLD']
		
		#Matches stars found with catalogue
		stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
		MatchTime = time.time()
		WCS = fits.open(FileWc)
		scale = float(WCS[0].header['COMMENT'][-8].split()[1])
		
		max_sep = float(config['max_sep'])*np.mean(FWHM)*scale/60
		idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars
		MatchTime = time.time() - MatchTime
		# ~ EncontradosX,EncontradosY,EncontradosMags = X[found],Y[found],Mags[found] 									#Saves X,Y and Mags positions found

		vmin = None
		vmax = None
		fig,axes = plt.subplots(1,2)
		ax = axes.ravel()
		ax[0].imshow(Data,cmap='gray',vmin=vmin or np.nanpercentile(Data, 0.1),vmax=vmax or np.nanpercentile(Data, 99))
		ax[1].imshow(Data,cmap='gray',vmin=vmin or np.nanpercentile(Data, 0.1),vmax=vmax or np.nanpercentile(Data, 99))
		for contador in range(len(Y)):
			if found[contador]:
				c1 = plt.Circle((X[contador],Y[contador]),0.4, color = 'green', linewidth=2,fill=False)
				ax[1].add_patch(c1)
			else:
				c1 = plt.Circle((X[contador],Y[contador]),0.4, color = 'red', linewidth=2,fill=False)
				ax[1].add_patch(c1)
				
		plt.show()
