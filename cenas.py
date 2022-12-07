# ~ from astroquery.vizier import Vizier
# ~ from astropy.io import fits
# ~ Vizier.ROW_LIMIT = 99999999
# ~ Catalogo = Vizier.get_catalogs('UCAC5')
# ~ Catalogo[0].write('teste.fits')


from Funcoes import FuncoesReducao

# ~ FileIn ='/home/luis/Desktop/Catalogs/UCAC5/z001'
# ~ FileOutput ='/home/luis/Desktop/Codigos/Pipeline/Files/Catalogos/UCAC5/z001.txt'
# ~ FuncoesReducao.UCAC5ToTxt(FileIn,FileOutput)





# ~ FileIn ='/home/luis/Desktop/Catalogos/CDS-322A/UCAC4/u4b/z'
# ~ FileOutput ='/home/luis/Desktop/Codigos/Pipeline/Files/Catalogos/UCAC4/'
# ~ for i in range(1,901):
	# ~ numero = f'{i:03d}'
	# ~ print(numero)
	# ~ FuncoesReducao.UCAC4ToTxt(FileIn+str(numero),FileOutput+str(numero)+'.txt')




FileIn ='/home/luis/Desktop/Catalogs/UCAC5/z'
FileOutput ='/home/luis/Desktop/Codigos/Pipeline/Files/Catalogos/UCAC5/'
for i in range(1,901):
	numero = f'{i:03d}'
	print(numero)
	FuncoesReducao.UCAC5ToTxt(FileIn+str(numero),FileOutput+str(numero)+'.txt')
