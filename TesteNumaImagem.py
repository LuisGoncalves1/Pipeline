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

print('''
Ficheiro para testar nas imagens do departamento da defesa com FicheiroSex

''')

directory = os.getcwd()
imagem = '/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit'
# ~ imagem = '/home/luis/Desktop/Codigos/Pipeline/Files/Imagens/2022-03-20/ImagensReduzidas/1037-293B_001.fits'

raw = CCDData.read(imagem, unit = 'adu') 
nomeficheiro = imagem.split('/')[-1].split('.')[0]


resultado = FuncoesReducao.sextractor(imagem,directory,nomeficheiro + '.cat')
result,status = FuncoesReducao.astrometrytesteComFicheiroSex(nomeficheiro+ '.cat',raw.header,directory + '/AstrometryDefesa','MAPA_SURV_4_085283_220417_114004+060023_scan047_02')


catalogo, mag = FuncoesReducao.ReadCatalogue(directory + '/Files/hipparcos.fits')

# ~ catalog_path = 'Files/catalog.dat'
# ~ catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)

FileWc = directory + '/AstrometryDefesa/' + nomeficheiro + '.wcs'
FileXY = nomeficheiro + '.cat'
FileOutput = nomeficheiro
result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)

Posicoes = fits.open(nomeficheiro)
Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
SexData = FuncoesReducao.ReadSextractorFileFit(nomeficheiro + '.cat')
X,Y,Mags = SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR']

print(SexData['FWHM_IMAGE'].data)


Corr = fits.open(directory + '/AstrometryDefesa/'+ nomeficheiro +'.corr' )
Ra2,Dec2,X2,Y2,Mags2 = Corr[1].data['field_ra'],Corr[1].data['field_dec'],Corr[1].data['field_x'],Corr[1].data['field_y'],Corr[1].data['MAG_ISOCOR']




stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')


idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=0.02 * u.deg)
EncontradosX,EncontradosY,EncontradosMags = X[found],Y[found],Mags[found]


print(idx)
print('-------------------------------------')
print('Diferença entre mags, minimizado')
print(mag[idx]-Mags[found] - min(mag[idx]-Mags[found] ))
# ~ print(mag[idx].data-Mags[found] - min(mag[idx].data-Mags[found] ))
print('-------------------------------------')
print('Mags astronomy, minimizado')
print(Mags[found] - min(Mags[found]))
print('-------------------------------------')
print('Mags catalogo, minimizado')
print(mag[idx] - min(mag[idx]))
# ~ print(mag[idx].data - min(mag[idx].data))



fig, axes = plt.subplots(1, 3)
vmin,vmax = None,None
axes[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
axes[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
axes[2].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
axes[0].set_title('Sextractor - %s contagens' %(len(SexData['X_IMAGE'])))
axes[1].set_title('Astrometry - %s contagens' %(len(X2)))
axes[2].set_title('Match')

for i in range(len(SexData['X_IMAGE'])):
	x,y,r = float(SexData['X_IMAGE'][i]),float(SexData['Y_IMAGE'][i]),2*float(SexData['FWHM_IMAGE'][i])
	axes[0].scatter(x, y, r, marker='o', facecolors='none',color='Blue')


for i in range(len(X2)):
	x,y,r = X2[i],Y2[i],2
	axes[1].scatter(x, y, r, marker='o', facecolors='none',color='Blue')
	
for i in range(len(EncontradosX)):
	x,y,r = EncontradosX[i],EncontradosY[i],2
	axes[2].scatter(x, y, r, marker='o', facecolors='none',color='Blue')


plt.show()
