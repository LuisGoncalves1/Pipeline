from Funcoes import FuncoesReducao
import os
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import numpy as np
import time
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky
import astropy.units as u


'''
Ficheiro para testar nas nossas imagens

'''

directory = os.getcwd()
imagem = '/Files/Imagens/2022-03-20/ImagensReduzidas/1037-293R_002.fits'

raw = CCDData.read(directory + imagem, unit = 'adu') 
header = raw.header
print(header['RA'])
resultado = FuncoesReducao.sextractor(directory+imagem,directory,'test.cat')
result,status = FuncoesReducao.astrometry(directory + imagem,header,directory + '/Files/Astrometry')
print(result)


catalogo, mag = FuncoesReducao.ReadCatalogue(directory + '/Files/Astrometry/hipparcos.fits')
Corr = fits.open(directory + '/Files/Astrometry/1037-293R_002.corr')
Ra,Dec,X,Y = Corr[1].data['field_ra'],Corr[1].data['field_dec'],Corr[1].data['field_x'],Corr[1].data['field_y']

stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')
idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=0.5 * u.deg)




# ~ print('''
# ~ Ficheiro para testar nas imagens do departamento da defesa com Sex

# ~ ''')


# ~ directory = os.getcwd()
# ~ imagem = '/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit'

# ~ raw = CCDData.read(imagem, unit = 'adu') 

# ~ result,status = FuncoesReducao.astrometrytesteComSex(imagem,raw.header,directory + '/AstrometryDefesa')

# ~ print(result)

# ~ print('-------------')
# ~ print(status)











# ~ SexData = FuncoesReducao.ReadSextractorFile('test.cat')
# ~ SebData = FuncoesReducao.ReadSextractorFile('valores.txt')
# ~ fig, axes = plt.subplots(1, 3)
# ~ vmin,vmax = None,None
# ~ axes[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
# ~ axes[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
# ~ axes[2].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
# ~ axes[0].set_title('Sextractor - %s contagens' %(len(SexData['X_IMAGE'])))
# ~ axes[1].set_title('Sep - %s contagens' %(len(SebData['X'])))
# ~ axes[2].set_title('Deteções utilizadas para o astrometry')


# ~ for i in range(len(SexData['X_IMAGE'])):
	# ~ x,y,r = float(SexData['X_IMAGE'][i]),float(SexData['Y_IMAGE'][i]),float(SexData['FLUX_ISOCOR'][i])/2000
	# ~ axes[0].scatter(x, y, r, marker='o', facecolors='none',color='Blue')


# ~ for i in range(len(SebData['X'])):
	# ~ x,y,r = float(SebData['X'][i]),float(SebData['Y'][i]),float(SebData['FLUX'][i])/1000
	# ~ axes[1].scatter(x, y, r, marker='o', facecolors='none',color='Blue')


# ~ for i in range(len(X)):
	# ~ x,y,r = X[i],Y[i],2
	# ~ axes[2].scatter(x, y, r, marker='o', facecolors='none',color='Blue')

# ~ plt.show()
