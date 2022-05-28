from Funcoes import FuncoesReducao
from astropy.io import fits

catalog_path = 'Files/Catalogos/ducati2002.fit'
catalogo,dados = FuncoesReducao.ReadDucati(catalog_path)
print(dados)








# ~ A = fits.open(catalog_path)
# ~ Data = A[1].data
# ~ Ra = Data['RAJ1950']
# ~ Dec = Data['DEJ1950']
# ~ Ragraus = FuncoesReducao.HMS2deg(Ra)
# ~ Decgraus = FuncoesReducao.DMS2deg(Dec)
# ~ print(Ragraus)
# ~ print(Decgraus)
