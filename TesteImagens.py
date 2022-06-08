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
config = FuncoesReducao.ReadPipelineConfig(directory+"/Files/Config/ConfigurationDefesa.txt")
catalog_path = 'Files/Catalogos/tycho2.dat'
# ~ catalogo,mag,BTMag = FuncoesReducao.ReadTychoCatalogue(catalog_path)
print('lido em %s' %(time.time()-A))


imagens = ['/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fit',
'/home/luis/Desktop/imagensMadeiraSV/Reduzidas/MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fit']

# ~ imagens = ['/home/luis/Desktop/Codigos/Pipeline/Files/Imagens/2022-03-20/Imagens/MAPA_SURV_4_085283_220417_114004+060023_scan047_01.fit',
# ~ '/home/luis/Desktop/Codigos/Pipeline/Files/Imagens/2022-03-20/Imagens/MAPA_SURV_4_085283_220417_114004+060023_scan047_02.fit',
# ~ '/home/luis/Desktop/Codigos/Pipeline/Files/Imagens/2022-03-20/Imagens/MAPA_SURV_4_085283_220417_114004+060023_scan047_03.fit',
# ~ '/home/luis/Desktop/Codigos/Pipeline/Files/Imagens/2022-03-20/Imagens/MAPA_SURV_4_085283_220417_114004+060023_scan047_04.fit']

fig, axes = plt.subplots(3, 4)
vmin,vmax = None,None
for contador in range(len(imagens)):
	imagem = imagens[contador]
	raw = CCDData.read(imagem, unit = 'adu') 

	if contador == 0:
		imagem2 = raw
		
	nomeficheiro = imagem.split('/')[-1].split('.')[0]
	resultado = FuncoesReducao.sextractor(imagem,directory,nomeficheiro + '.cat')
	# ~ result,status = FuncoesReducao.astrometry(nomeficheiro+ '.cat',raw.header,directory + '/AstrometryDefesa',config)

	FileWc = directory + '/AstrometryDefesa/' + nomeficheiro + '.wcs'
	FileXY = nomeficheiro + '.cat'
	FileOutput = nomeficheiro
	# ~ result, status	= FuncoesReducao.XYToRaDec(FileWc,FileXY,FileOutput)

	# ~ Posicoes = fits.open(nomeficheiro)
	# ~ Ra,Dec = Posicoes[1].data['RA'],Posicoes[1].data['DEC']
	SexData = FuncoesReducao.ReadSextractorFileFit(nomeficheiro + '.cat')
	X,Y,Mags,MagsErro,BackG,fluxo,FWHM,thresh= SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'],SexData['MAGERR_ISOCOR'],SexData['BACKGROUND'],SexData['FLUX_ISOCOR'],SexData['FWHM_IMAGE'],SexData['THRESHOLD']


	# ~ Corr = fits.open(directory + '/AstrometryDefesa/'+ nomeficheiro +'.corr' )
	# ~ Ra2,Dec2,X2,Y2,Mags2 = Corr[1].data['field_ra'],Corr[1].data['field_dec'],Corr[1].data['field_x'],Corr[1].data['field_y'],Corr[1].data['MAG_ISOCOR']
	# ~ stars = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')


	# ~ idx,found = FuncoesReducao.find_matching_stars(stars, catalogo, max_sep=0.1 * u.deg)
	# ~ EncontradosX,EncontradosY,EncontradosMags = X[found],Y[found],Mags[found]
	
	ImagemDiferença = np.array(raw)-np.array(imagem2)
	Binarizacao = ImagemDiferença > 1000
	axes[0,contador].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))	
	axes[1,contador].imshow(ImagemDiferença,cmap='gray',vmin=vmin or np.nanpercentile(raw, 0.1),vmax=vmax or np.nanpercentile(raw, 99))	
	axes[2,contador].imshow(Binarizacao,cmap = 'gray',vmin = 0,vmax=1)	
	for i in range(len(X)):
		x,y,r = X[i],Y[i],FWHM[i]
		axes[0,contador].scatter(x, y, r, marker='o', facecolors='none',color='Blue')
	
	
	
	
plt.show()

# Ver: Numero de fontes em cada imagem, diferença entre fontes e magnitude de diferenças

