from Funcoes import FuncoesReducao
import os
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import numpy as np
import time
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky
import astropy.units as u5

catalog_path = 'Files/catalog.dat'
A = time.time()
directory = os.getcwd()
# ~ catalogo, mag = FuncoesReducao.ReadCatalogue(directory + '/Files/hipparcos.fits')
# ~ catalog,VTMag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)
catalog = FuncoesReducao.ReadTychoCatalogue(catalog_path)

print(time.time()-A)
print('----------------------------')
print(max(catalog['RA']))
print(min(catalog['RA']))
print('----------------------------')
print(max(catalog['DEC']))
print(min(catalog['DEC']))
print('----------------------------')
print(max(catalog['BTMag']))
print(min(catalog['BTMag']))
print('----------------------------')
print(max(catalog['VTMag']))
print(min(catalog['VTMag']))
print('----------------------------')



# ~ Teste = []
# ~ for i in range(len(catalog['BTMag'])):
	# ~ if catalog['BTMag'][i] == '':
		# ~ Teste.append(i)

# ~ print(Teste)
		



# ~ FileWc = '/home/luis/Desktop/Codigos/Pipeline/AstrometryDefesa/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.wcs'
# ~ FileXY = 'MAPA_SURV_4_085283_220417_114004+060023_scan047_02.cat'
# ~ FileOutput = 'Teste'
# ~ result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)
