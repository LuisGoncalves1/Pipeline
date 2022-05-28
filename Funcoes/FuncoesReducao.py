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
	cmd['-CATALOG_NAME'] =  catout
	cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items())))
	result = subprocess.run(list_cmd, capture_output=True)
	return result
	
def astrometry(DirImagem,header,DirAstro):
	cmd = ['solve-field', '--no-plots',
           '--no-remove-lines',
           '--uniformize', '0',
           '--width', str(header['NAXIS1']),
           '--height', str(header['NAXIS2']),           
           '--crpix-center',
           '--cpulimit', '120',
           '--scale-low', '0.1',
           '--scale-high', '1',
           '--scale-units', 'app',
           '--ra', header['RA'],
           '--dec', header['DEC'],
           '--tweak-order','2',
           '--radius', '1.6',
           '--overwrite',
           '--dir', DirAstro,
            DirImagem]
            
    # log cmd executed
	cmd_str = ' '.join(cmd)
	print(cmd_str)
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode
	return result, status
	      		
def astrometrytesteComSex(DirImagem,header,DirAstro):
	cmd = ['solve-field', '--no-plots',
           '--no-remove-lines',
           '--uniformize', '0',
           '--width', str(header['NAXIS1']),
           '--height', str(header['NAXIS2']),           
           '--crpix-center',
           '--cpulimit', '120',
           '--scale-low', '1',
           '--scale-high', '5',
           '--scale-units', 'app',
           '--ra', '11:40:04',
           '--dec', '+06:00:23',
           '--tweak-order','2',
           '--use-sextractor',
           '--sextractor-config', '/home/luis/Desktop/Codigos/Pipeline/Files/Sextractor/default.sex',
           '--radius', '3',
           '--overwrite',
           '--dir', DirAstro,
           '-s', 'MAG_ISOCOR',
           '-a', '-X', 'X_IMAGE', '-Y', 'Y_IMAGE',
            DirImagem
    ]
    # log cmd executed
	cmd_str = ' '.join(cmd)
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode
	return result, status	      		

def astrometrytesteComFicheiroSex(DirCatalogo,header,DirAstro,NameImage):
	cmd = ['solve-field', '--no-plots',
	       '--x-column', 'X_IMAGE',
           '--y-column', 'Y_IMAGE',
           '--no-remove-lines',
           '--uniformize', '0',
           '--width', str(header['NAXIS1']),
           '--height', str(header['NAXIS2']),           
           '--crpix-center',
           '--cpulimit', '120',
           '--scale-low', '1',
           '--scale-high', '5',
           # ~ '--scale-low', '0.1',
           # ~ '--scale-high', '1',
           '--scale-units', 'app',
           # ~ '--ra', header['RA'],
           # ~ '--dec', header['DEC'],
           '--ra', '11:40:04',
           '--dec', '+06:00:23',
           '--tweak-order','2',
           '--rdls', 'none',
           # ~ '--wcs', 'none',
           '--solved','none',
           '--axy','none',
           '--match','none',
           # ~ '--new-fits', DirAstro + NameImage,
           '--radius', '3',
           # ~ '--radius', '1.6',
           '--overwrite',
           '--dir', DirAstro,
            DirCatalogo
    ]
    # log cmd executed
	cmd_str = ' '.join(cmd)
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode
	return result, status	      		
	      		
def XYToRaDec(FileWc,FileXY,FileOutput):
	cmd = ['wcs-xy2rd', 
		'-w',FileWc,
		'-i', FileXY,
		'-o', FileOutput,
		'-X','X_IMAGE',
		'-Y','Y_IMAGE']
		
	cmd_str = ' '.join(cmd)
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode
	return result, status	           		
	      		
def ReadSextractorFileTxt(SexFile):
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
	config = {}
	Ficheiro = fits.open(SexFile)
	Dados = Ficheiro[1]
	Colunas = Dados.columns
	for i in Colunas:
		Coluna = i.name
		config[Coluna] = Dados.data[Coluna]
		
		
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

def AnalisarImagem(ficheiro,Masters,logging,config,diretoria,directory):
	logging.debug("Startins to analise " + ficheiro)
	Data = CCDData.read(diretoria + "/ImagensTratadas/" + ficheiro, unit = 'adu') 
	raw = Data.data.copy()
	Filter = Data.header['FILTER']
	# ~ exposure = Data.header['EXPTIME']
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
		
	
	DirImage = diretoria  + '/ImagensReduzidas/' + ficheiro		
	Data.write(DirImage)
	
	
	if config["Astrometry"]:
		resultado = FuncoesReducao.sextractor(DirImage,directory,diretoria + '/SextractorCatalogue/%s.cat'(ficheiro[:-5]))
		result,status = FuncoesReducao.astrometrytesteComFicheiroSex(diretoria + '/SextractorCatalogue/%s.cat'(ficheiro[:-5]),raw.header,directory + '/AstrometryFiles')
	
	
	



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
	
	print(type(catalog['ra']))
	# ~ return SkyCoord(ra=Angle(catalog['ra']), dec=Angle(catalog['dec']))
	return SkyCoord(ra = Longitude(catalog['ra'], unit = u.deg), dec = np.array(catalog['dec'])*u.deg,frame = 'cirs'), catalog['v_mag']

def ReadDucati(catalog_path,max_magnitude=None,min_magnitude=None):
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
	Input data: +HH:MM:SS or -HH MM SS, returns decimal degrees correspondent
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
	Input data: +DD:MM:SS or -DD MM SS, returns decimal degrees correspondent
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
    if len(image_stars) == 0:
        return np.array([], dtype=int), np.zeros(len(catalog_stars), dtype=bool)

    idx, d2d, d3d = match_coordinates_sky(catalog_stars, image_stars,storekdtree = True)
    mask = d2d < max_sep

    return idx[mask], mask
