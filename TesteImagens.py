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

A = time.time()
directory = os.getcwd()
catalog_path = 'Files/Catalogos/tycho2.dat'
catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)
print('lido em %s' %(time.time()-A))


imagens = ['/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit']

fig, axes = plt.subplots(1, 4)
vmin,vmax = None,None
for contador in range(len(imagens)):
	imagem = imagens[contador]
	raw = CCDData.read(imagem, unit = 'adu') 

	nomeficheiro = imagem.split('/')[-1].split('.')[0]
	resultado = FuncoesReducao.sextractor(imagem,directory,nomeficheiro + '.cat')
	result,status = FuncoesReducao.astrometrytesteComFicheiroSex(nomeficheiro+ '.cat',raw.header,directory + '/AstrometryDefesa',nomeficheiro)

	FileWc = directory + '/AstrometryDefesa/' + nomeficheiro + '.wcs'
	FileXY = nomeficheiro + '.cat'
	FileOutput = nomeficheiro
	result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)

	Posicoes = fits.open(nomeficheiro)
	Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
	SexData = FuncoesReducao.ReadSextractorFileFit(nomeficheiro + '.cat')
	X,Y,Mags,MagsErro,BackG,fluxo,FWHM,thresh= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE'],SexData['THRESHOLD']


	Corr = fits.open(directory + '/AstrometryDefesa/'+ nomeficheiro +'.corr' )
	Ra2,Dec2,X2,Y2,Mags2 = Corr[1].data['field_ra'],Corr[1].data['field_dec'],Corr[1].data['field_x'],Corr[1].data['field_y'],Corr[1].data['MAG_ISOCOR']
	stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')


	idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=0.1 * u.deg)
	EncontradosX,EncontradosY,EncontradosMags = X[found],Y[found],Mags[found]
	axes[contador].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))
	print(len(found))
	print(len(SexData['Y_IMAGE']))
	
	for i in range(len(EncontradosX)):
		x,y,r = EncontradosX[i],EncontradosY[i],2
		axes[contador].scatter(x, y, r, marker='o', facecolors='none',color='Blue')

	
	
plt.show()

