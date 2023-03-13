import ccdproc 
from datetime import date
import datetime
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
from astropy.coordinates import SkyCoord, Longitude,  Angle, match_coordinates_sky, AltAz, EarthLocation
import struct
import requests
import json
# CATALOG FUNCTIONS

def UCAC5ToTxtTotal(FileIn,FileOutput):
	'''
	Function that reads the UCAC5 original file and returns a Txt file with all information
	
	INPUTS:
	FileIn 				-> Path to UCAC5 zone file									(str)
	FileOutput		 	-> Path to output file name									(str)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	FileOutput			-> File with all UCAC5 information							(txt)
	'''
	with open(FileIn,'rb') as ficheiro:
		A = ficheiro.read()
	
	with open(FileOutput,'a') as ficheiro:	
		for i in range(int(len(A)/52)):
			DadosEstrela = A[i*52:(i+1)*52]
			
			#Set de 15 bytes (15)
			srcid = struct.unpack('ii',DadosEstrela[0:8])[0]
			Ra = struct.unpack('i',DadosEstrela[8:12])[0]
			Dec = struct.unpack('i',DadosEstrela[12:16])[0]
			Erg = int.from_bytes(DadosEstrela[16:18],'little')
			Edg = int.from_bytes(DadosEstrela[18:20],'little')
			flg = int(DadosEstrela[20])
			nu = int(DadosEstrela[21])
			epu = int.from_bytes(DadosEstrela[22:24],'little')
			ira = struct.unpack('i',DadosEstrela[24:28])[0]
			idc = struct.unpack('i',DadosEstrela[28:32])[0]
			pmur = int.from_bytes(DadosEstrela[32:34],'little')
			pmud = int.from_bytes(DadosEstrela[34:36],'little')
			pmer = int.from_bytes(DadosEstrela[36:38],'little')
			pmed = int.from_bytes(DadosEstrela[38:40],'little')
			
			gmag = int.from_bytes(DadosEstrela[40:42],'little')
			umag = int.from_bytes(DadosEstrela[42:44],'little')
			rmag = int.from_bytes(DadosEstrela[44:46],'little')
			jmag = int.from_bytes(DadosEstrela[46:48],'little')
			hmag = int.from_bytes(DadosEstrela[48:50],'little')
			kmag = int.from_bytes(DadosEstrela[50:52],'little')
		
			linha2 = f'{srcid:10},{Ra:11},{Dec:11},{Erg:5},{Edg:5},{flg:2},{nu:2},{epu:4},{ira:11},{idc:11},{pmur:5},{pmud:5},{pmer:5},{pmed:5},{gmag:5},{umag:5},{rmag:5},{jmag:5},{hmag:5},{kmag:5}\n'
			ficheiro.write(linha2)

def UCAC5ToTxt(FileIn,FileOutput):
	'''
	Function that reads the UCAC5 original file and returns a Txt file with Ra, Dec and g/u/r/j/h/k magnitudes
	
	INPUTS:
	FileIn 				-> Path to UCAC5 zone file									(str)
	FileOutput		 	-> Path to output file name									(str)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	FileOutput			-> File with all UCAC5 information							(txt)
	'''
	with open(FileIn,'rb') as ficheiro:
		A = ficheiro.read()
	
	with open(FileOutput,'a') as ficheiro:	
		for i in range(int(len(A)/52)):
			DadosEstrela = A[i*52:(i+1)*52]
			
			#Set de 15 bytes (15)
			Ra = struct.unpack('i',DadosEstrela[8:12])[0]
			Dec = struct.unpack('i',DadosEstrela[12:16])[0]	
			gmag = int.from_bytes(DadosEstrela[40:42],'little')
			umag = int.from_bytes(DadosEstrela[42:44],'little')
			rmag = int.from_bytes(DadosEstrela[44:46],'little')
			jmag = int.from_bytes(DadosEstrela[46:48],'little')
			hmag = int.from_bytes(DadosEstrela[48:50],'little')
			kmag = int.from_bytes(DadosEstrela[50:52],'little')
		
			linha2 = f'{Ra:11},{Dec:11},{gmag:5},{umag:5},{rmag:5},{jmag:5},{hmag:5},{kmag:5}\n'
			ficheiro.write(linha2)

def UCAC4ToTxtTotal(FileIn,FileOutput):
	'''
	Function that reads the UCAC4 original file and returns a Txt file with all information
	
	INPUTS:
	FileIn 				-> Path to UCAC4 zone file									(str)
	FileOutput		 	-> Path to output file name									(str)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	FileOutput			-> File with all UCAC4 information							(txt)
	'''
	with open(FileIn,'rb') as ficheiro:
		A = ficheiro.read()
	
	with open(FileOutput,'a') as ficheiro:	
		for i in range(int(len(A)/78)):
			DadosEstrela = A[i*78:(i+1)*78]
			
			#Set de 15 bytes (15)
			Ra = struct.unpack('i',DadosEstrela[0:4])[0]
			Dec = struct.unpack('i',DadosEstrela[4:8])[0]
			Magm = int.from_bytes(DadosEstrela[8:10],'little')
			Maga = int.from_bytes(DadosEstrela[10:12],'little')
			SigmaG = int(DadosEstrela[12])
			ObjT = int(DadosEstrela[13])
			Cdf = int(DadosEstrela[14])
			
			
			#Set de 5 bytes (20)
			sigra = int(DadosEstrela[15])
			if sigra > 128/2:
				sigra -= 128
			sigdc = int(DadosEstrela[16])
			if sigdc > 128/2:
				sigdc -= 128
			na1 = int(DadosEstrela[17])
			nu1 = int(DadosEstrela[18])
			cu1 = int(DadosEstrela[19])

			
			#Set de 10 bytes (30)
			cepra = int.from_bytes(DadosEstrela[20:22],'little')
			cepdc = int.from_bytes(DadosEstrela[22:24],'little')
			pmrac = int.from_bytes(DadosEstrela[24:26],'little')
			if pmrac > 65536/2:
				pmrac -= 65536
			pmdc = int.from_bytes(DadosEstrela[26:28],'little')
			if pmdc > 65536/2:
				pmdc -= 65536
			sigpmr = int(DadosEstrela[28])
			if sigpmr > 128/2:
				sigpmr -= 128
			sigpmd = int(DadosEstrela[29])
			if sigpmd > 128/2:
				sigpmd -= 128

			
			#Set de 16 bytes (46)
			pts_key = struct.unpack('i',DadosEstrela[30:34])[0]
			j_m = int.from_bytes(DadosEstrela[34:36],'little')
			h_m = int.from_bytes(DadosEstrela[36:38],'little')
			k_m = int.from_bytes(DadosEstrela[38:40],'little')
			icqflg1 = int(DadosEstrela[40])
			icqflg2 = int(DadosEstrela[41])
			icqflg3 = int(DadosEstrela[42])
			e2mpho1 = int(DadosEstrela[43])
			e2mpho2 = int(DadosEstrela[44])
			e2mpho3 = int(DadosEstrela[45])

			
			
			#Set de 16 bytes (62)
			apasm1 = int.from_bytes(DadosEstrela[46:48],'little')
			apasm2 = int.from_bytes(DadosEstrela[48:50],'little')
			apasm3 = int.from_bytes(DadosEstrela[50:52],'little')
			apasm4 = int.from_bytes(DadosEstrela[52:54],'little')
			apasm5 = int.from_bytes(DadosEstrela[54:56],'little')
			apase1 = int(DadosEstrela[56])
			apase2 = int(DadosEstrela[57])
			apase3 = int(DadosEstrela[58])
			apase4 = int(DadosEstrela[59])
			apase5 = int(DadosEstrela[60])
			gcflg = int(DadosEstrela[61])

			
		
			#Set de 4 bytes (66)
			icf = struct.unpack('i',DadosEstrela[62:66])[0]
			icf = f'{icf:09d}'
			icf1 = int(icf[0])
			icf2 = int(icf[1])
			icf3 = int(icf[2])
			icf4 = int(icf[3])
			icf5 = int(icf[4])
			icf6 = int(icf[5])
			icf7 = int(icf[6])
			icf8 = int(icf[7])
			icf9 = int(icf[8])

			
			#Set de 12 bytes (78)
			leda = int(DadosEstrela[66])
			x2m = int(DadosEstrela[67])
			rnm = struct.unpack('i',DadosEstrela[68:72])[0]
			zn2 = int.from_bytes(DadosEstrela[72:74],'little')
			rn2 = struct.unpack('i',DadosEstrela[74:78])[0]
	
			
			
			# ~ linha = '{10:},{10:},{5:},{5:},{2:},{1:},{2:},{3:},{3:},{2:},{2:},{2:},{5:},{5:},{6:},{6:},{3:},{3:},{10:},{5:},{5:},{5:},{2:},{2:},{2:},{2:},{2:},{2:},{5:},{5:},{5:},{5:},{5:},{3:},{3:},{3:},{3:},{3:},{1:},{1:},{1:},{1:},{1:},{1:},{1:},{1:},{1:},{1:},{2:},{2:},{9:},{3:},{6:}\n'.format(Ra,Dec,Magm,Maga,SigmaG,ObjT,Cdf,sigra,sigdc,na1,nu1,cu1,cepra,cepdc,pmrac,pmdc,sigpmr,sigpmd,pts_key,j_m,h_m,k_m,icqflg1,icqflg2,icqflg3,e2mpho1,e2mpho2,e2mpho3,apasm1,apasm2,apasm3,apasm4,apasm5,apase1,apase2,apase3,apase4,apase5,gcflg,icf1,icf2,icf3,icf4,icf5,icf6,icf7,icf8,icf9,leda,x2m,rnm,zn2,rn2)
			linha2 = f'{Ra:10},{Dec:9},{Magm:5},{Maga:5},{SigmaG:2},{ObjT:1},{Cdf:2},{sigra:3},{sigdc:3},{na1:2},{nu1:2},{cu1:2},{cepra:5},{cepdc:5},{pmrac:6},{pmdc:6},{sigpmr:3},{sigpmd:3},{pts_key:10},{j_m:5},{h_m:5},{k_m:5},{icqflg1:2},{icqflg2:2},{icqflg3:2},{e2mpho1:2},{e2mpho2:2},{e2mpho3:2},{apasm1:5},{apasm2:5},{apasm3:5},{apasm4:5},{apasm5:5},{apase1:3},{apase2:3},{apase3:3},{apase4:3},{apase5:3},{gcflg:1},{icf1:1},{icf2:1},{icf3:1},{icf4:1},{icf5:1},{icf6:1},{icf7:1},{icf8:1},{icf9:1},{leda:2},{x2m:2},{rnm:9},{zn2:3},{rn2:6}\n'
			ficheiro.write(linha2)

def UCAC4ToTxt(FileIn,FileOutput):
	'''
	Function that reads the UCAC4 original file and returns a Txt file with Ra,Dec and UCAC4m/a B/V/g/r/i magnitudes
	
	INPUTS:
	FileIn 				-> Path to UCAC5 zone file									(str)
	FileOutput		 	-> Path to output file name									(str)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	FileOutput			-> File with all UCAC5 information							(txt)
	'''
	with open(FileIn,'rb') as ficheiro:
		A = ficheiro.read()
	
	with open(FileOutput,'a') as ficheiro:
		# ~ linha2 = f'{"RA":10} {"DEC":9} {"Magm":5} {"Maga":5} {"B":5} {"V":5} {"g":5} {"r":5} {"i":5}\n'
		# ~ ficheiro.write(linha2)	
		for i in range(int(len(A)/78)):
			DadosEstrela = A[i*78:(i+1)*78]
			
			#Set de 15 bytes (15)
			Ra = struct.unpack('i',DadosEstrela[0:4])[0]#/(1000*60*60)
			Dec = struct.unpack('i',DadosEstrela[4:8])[0]#/(1000*60*60)
			Magm = int.from_bytes(DadosEstrela[8:10],'little')
			Maga = int.from_bytes(DadosEstrela[10:12],'little')
			
			#Set de 16 bytes (62)
			apasm1 = int.from_bytes(DadosEstrela[46:48],'little')
			apasm2 = int.from_bytes(DadosEstrela[48:50],'little')
			apasm3 = int.from_bytes(DadosEstrela[50:52],'little')
			apasm4 = int.from_bytes(DadosEstrela[52:54],'little')
			apasm5 = int.from_bytes(DadosEstrela[54:56],'little')
	
			
			
			linha2 = f'{Ra:10},{Dec:9},{Magm:5},{Maga:5},{apasm1:5},{apasm2:5},{apasm3:5},{apasm4:5},{apasm5:5}\n'
			ficheiro.write(linha2)
	
def ZonesUCAC(MinDec,MaxDec):
	'''
	Function that reads calculates which zone files to read from UCAC4/UCAC5 catalogues
	
	INPUTS:
	MinDec 				-> Minimum Dec value in field								(float)
	MaxDec		 		-> Maximum Dec value in field								(float)
	
	RETURNS:
	Zones				-> list with all zone files to read							(list)
	
	FILES CREATED:
	NONE
	'''
	Decs = np.array([-90.0, -89.8, -89.6, -89.4, -89.2, -89.0, -88.8, -88.6, -88.4, -88.2, -88.0, -87.8, -87.6, -87.4, -87.2, -87.0, -86.8, -86.6, -86.4, -86.2, -86.0, -85.8, -85.6, -85.4, -85.2, -85.0, -84.8, -84.6, -84.4, -84.2, -84.0, -83.8, -83.6, -83.4, -83.2, -83.0, -82.8, -82.6, -82.4, -82.2, -82.0, -81.8, -81.6, -81.4, -81.2, -81.0, -80.8, -80.6, -80.4, -80.2, -80.0, -79.8, -79.6, -79.4, -79.2, -79.0, -78.8, -78.6, -78.4, -78.2, -78.0, -77.8, -77.6, -77.4, -77.2, -77.0, -76.8, -76.6, -76.4, -76.2, -76.0, -75.8, -75.6, -75.4, -75.2, -75.0, -74.8, -74.6, -74.4, -74.2, -74.0, -73.8, -73.6, -73.4, -73.2, -73.0, -72.8, -72.6, -72.4, -72.2, -72.0, -71.8, -71.6, -71.4, -71.2, -71.0, -70.8, -70.6, -70.4, -70.2, -70.0, -69.8, -69.6, -69.4, -69.2, -69.0, -68.8, -68.6, -68.4, -68.2, -68.0, -67.8, -67.6, -67.4, -67.2, -67.0, -66.8, -66.6, -66.4, -66.2, -66.0, -65.8, -65.6, -65.4, -65.2, -65.0, -64.8, -64.6, -64.4, -64.2, -64.0, -63.8, -63.6, -63.4, -63.2, -63.0, -62.8, -62.6, -62.4, -62.2, -62.0, -61.8, -61.6, -61.4, -61.2, -61.0, -60.8, -60.6, -60.4, -60.2, -60.0, -59.8, -59.6, -59.4, -59.2, -59.0, -58.8, -58.6, -58.4, -58.2, -58.0, -57.8, -57.6, -57.4, -57.2, -57.0, -56.8, -56.6, -56.4, -56.2, -56.0, -55.8, -55.6, -55.4, -55.2, -55.0, -54.8, -54.6, -54.4, -54.2, -54.0, -53.8, -53.6, -53.4, -53.2, -53.0, -52.8, -52.6, -52.4, -52.2, -52.0, -51.8, -51.6, -51.4, -51.2, -51.0, -50.8, -50.6, -50.4, -50.2, -50.0, -49.8, -49.6, -49.4, -49.2, -49.0, -48.8, -48.6, -48.4, -48.2, -48.0, -47.8, -47.6, -47.4, -47.2, -47.0, -46.8, -46.6, -46.4, -46.2, -46.0, -45.8, -45.6, -45.4, -45.2, -45.0, -44.8, -44.6, -44.4, -44.2, -44.0, -43.8, -43.6, -43.4, -43.2, -43.0, -42.8, -42.6, -42.4, -42.2, -42.0, -41.8, -41.6, -41.4, -41.2, -41.0, -40.8, -40.6, -40.4, -40.2, -40.0, -39.8, -39.6, -39.4, -39.2, -39.0, -38.8, -38.6, -38.4, -38.2, -38.0, -37.8, -37.6, -37.4, -37.2, -37.0, -36.8, -36.6, -36.4, -36.2, -36.0, -35.8, -35.6, -35.4, -35.2, -35.0, -34.8, -34.6, -34.4, -34.2, -34.0, -33.8, -33.6, -33.4, -33.2, -33.0, -32.8, -32.6, -32.4, -32.2, -32.0, -31.8, -31.6, -31.4, -31.2, -31.0, -30.8, -30.6, -30.4, -30.2, -30.0, -29.8, -29.6, -29.4, -29.2, -29.0, -28.8, -28.6, -28.4, -28.2, -28.0, -27.8, -27.6, -27.4, -27.2, -27.0, -26.8, -26.6, -26.4, -26.2, -26.0, -25.8, -25.6, -25.4, -25.2, -25.0, -24.8, -24.6, -24.4, -24.2, -24.0, -23.8, -23.6, -23.4, -23.2, -23.0, -22.8, -22.6, -22.4, -22.2, -22.0, -21.8, -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2, -20.0, -19.8, -19.6, -19.4, -19.2, -19.0, -18.8, -18.6, -18.4, -18.2, -18.0, -17.8, -17.6, -17.4, -17.2, -17.0, -16.8, -16.6, -16.4, -16.2, -16.0, -15.8, -15.6, -15.4, -15.2, -15.0, -14.8, -14.6, -14.4, -14.2, -14.0, -13.8, -13.6, -13.4, -13.2, -13.0, -12.8, -12.6, -12.4, -12.2, -12.0, -11.8, -11.6, -11.4, -11.2, -11.0, -10.8, -10.6, -10.4, -10.2, -10.0, -9.8, -9.6, -9.4, -9.2, -9.0, -8.8, -8.6, -8.4, -8.2, -8.0, -7.8, -7.6, -7.4, -7.2, -7.0, -6.8, -6.6, -6.4, -6.2, -6.0, -5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2, -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2, 16.4, 16.6, 16.8, 17.0, 17.2, 17.4, 17.6, 17.8, 18.0, 18.2, 18.4, 18.6, 18.8, 19.0, 19.2, 19.4, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2, 21.4, 21.6, 21.8, 22.0, 22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6, 23.8, 24.0, 24.2, 24.4, 24.6, 24.8, 25.0, 25.2, 25.4, 25.6, 25.8, 26.0, 26.2, 26.4, 26.6, 26.8, 27.0, 27.2, 27.4, 27.6, 27.8, 28.0, 28.2, 28.4, 28.6, 28.8, 29.0, 29.2, 29.4, 29.6, 29.8, 30.0, 30.2, 30.4, 30.6, 30.8, 31.0, 31.2, 31.4, 31.6, 31.8, 32.0, 32.2, 32.4, 32.6, 32.8, 33.0, 33.2, 33.4, 33.6, 33.8, 34.0, 34.2, 34.4, 34.6, 34.8, 35.0, 35.2, 35.4, 35.6, 35.8, 36.0, 36.2, 36.4, 36.6, 36.8, 37.0, 37.2, 37.4, 37.6, 37.8, 38.0, 38.2, 38.4, 38.6, 38.8, 39.0, 39.2, 39.4, 39.6, 39.8, 40.0, 40.2, 40.4, 40.6, 40.8, 41.0, 41.2, 41.4, 41.6, 41.8, 42.0, 42.2, 42.4, 42.6, 42.8, 43.0, 43.2, 43.4, 43.6, 43.8, 44.0, 44.2, 44.4, 44.6, 44.8, 45.0, 45.2, 45.4, 45.6, 45.8, 46.0, 46.2, 46.4, 46.6, 46.8, 47.0, 47.2, 47.4, 47.6, 47.8, 48.0, 48.2, 48.4, 48.6, 48.8, 49.0, 49.2, 49.4, 49.6, 49.8, 50.0, 50.2, 50.4, 50.6, 50.8, 51.0, 51.2, 51.4, 51.6, 51.8, 52.0, 52.2, 52.4, 52.6, 52.8, 53.0, 53.2, 53.4, 53.6, 53.8, 54.0, 54.2, 54.4, 54.6, 54.8, 55.0, 55.2, 55.4, 55.6, 55.8, 56.0, 56.2, 56.4, 56.6, 56.8, 57.0, 57.2, 57.4, 57.6, 57.8, 58.0, 58.2, 58.4, 58.6, 58.8, 59.0, 59.2, 59.4, 59.6, 59.8, 60.0, 60.2, 60.4, 60.6, 60.8, 61.0, 61.2, 61.4, 61.6, 61.8, 62.0, 62.2, 62.4, 62.6, 62.8, 63.0, 63.2, 63.4, 63.6, 63.8, 64.0, 64.2, 64.4, 64.6, 64.8, 65.0, 65.2, 65.4, 65.6, 65.8, 66.0, 66.2, 66.4, 66.6, 66.8, 67.0, 67.2, 67.4, 67.6, 67.8, 68.0, 68.2, 68.4, 68.6, 68.8, 69.0, 69.2, 69.4, 69.6, 69.8, 70.0, 70.2, 70.4, 70.6, 70.8, 71.0, 71.2, 71.4, 71.6, 71.8, 72.0, 72.2, 72.4, 72.6, 72.8, 73.0, 73.2, 73.4, 73.6, 73.8, 74.0, 74.2, 74.4, 74.6, 74.8, 75.0, 75.2, 75.4, 75.6, 75.8, 76.0, 76.2, 76.4, 76.6, 76.8, 77.0, 77.2, 77.4, 77.6, 77.8, 78.0, 78.2, 78.4, 78.6, 78.8, 79.0, 79.2, 79.4, 79.6, 79.8, 80.0, 80.2, 80.4, 80.6, 80.8, 81.0, 81.2, 81.4, 81.6, 81.8, 82.0, 82.2, 82.4, 82.6, 82.8, 83.0, 83.2, 83.4, 83.6, 83.8, 84.0, 84.2, 84.4, 84.6, 84.8, 85.0, 85.2, 85.4, 85.6, 85.8, 86.0, 86.2, 86.4, 86.6, 86.8, 87.0, 87.2, 87.4, 87.6, 87.8, 88.0, 88.2, 88.4, 88.6, 88.8, 89.0, 89.2, 89.4, 89.6, 89.8, 90.0])
	mask = (Decs>MinDec) & (Decs<MaxDec)
	Zones =  [i for i, x in enumerate(mask) if x]
	Zones.append(Zones[-1]+1)
	return Zones

def ReadUCAC4(Zones,Folder,MagLim):
	'''
	Function that reads UCAC4 txt zone files and returns SkyCoord Object with Ra/Dec and Dictionary with all magnitudes
	
	INPUTS:
	Zones 				-> list with all zone files to read	 						(float)
	Folder		 		-> Folder with UCAC4 zone files								(float)
	MagLim		 		-> Maximum magnitude to read								(float)
	
	RETURNS:
	Catalog				-> SkyCoord object with Ra/Dec information					(SkyCoord)
	Mags				-> Dictionary with magnitudes corresponding with Catalog	(dict)
	
	FILES CREATED:
	NONE
	'''
	Folder += 'UCAC4/'
	Ra = []
	Dec = []
	Mags = {'UCAC4':[],'B':[],'V':[],'g':[],'r':[],'i':[]}
	for i in Zones:
		zone = f'{i:03d}'
		with open(Folder + zone + '.txt','r') as text:
			while True:
				line = text.readline()
				if line == '':
					break
				
				line = line.strip().split(',')
				Mag = float(line[2])/1000
				if Mag < MagLim:
					Ra.append(int(line[0])/(60*60*1000))
					Dec.append(int(line[1])/(60*60*1000)-90)
					Mags['UCAC4'].append(Mag)
					Mags['B'].append(float(line[4])/1000)
					Mags['V'].append(float(line[5])/1000)
					Mags['g'].append(float(line[6])/1000)
					Mags['r'].append(float(line[7])/1000)
					Mags['i'].append(float(line[8])/1000)
			
	return SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg,frame = 'cirs'), Mags

def ReadUCAC5(Zones,Folder,MagLim):
	'''
	Function that reads UCAC5 txt zone files and returns SkyCoord Object with Ra/Dec and Dictionary with all magnitudes
	
	INPUTS:
	Zones 				-> list with all zone files to read	 						(float)
	Folder		 		-> Folder with UCAC4 zone files								(float)
	MagLim		 		-> Maximum magnitude to read								(float)
	
	RETURNS:
	Catalog				-> SkyCoord object with Ra/Dec information					(SkyCoord)
	Mags				-> Dictionary with magnitudes corresponding with Catalog	(dict)
	
	FILES CREATED:
	NONE
	'''
	Ra = []
	Dec = []
	Folder += 'UCAC5/'
	Mags = {'G':[],'U':[],'R':[],'J':[],'H':[],'K':[]}
	for i in Zones:
		zone = f'{i:03d}'
		with open(Folder + zone + '.txt','r') as text:
			while True:
				line = text.readline()
				if line == '':
					break
				
				line = line.strip().split(',')
				Mag = float(line[2])/1000
				if Mag < MagLim:
					Ra.append(int(line[0])/(60*60*1000))
					Dec.append(int(line[1])/(60*60*1000))
					Mags['G'].append(Mag)
					Mags['U'].append(float(line[3])/1000)
					Mags['R'].append(float(line[4])/1000)
					Mags['J'].append(float(line[5])/1000)
					Mags['H'].append(float(line[6])/1000)
					Mags['K'].append(float(line[7])/1000)
			
	return SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg,frame = 'cirs'), Mags
		
def ReadHipparcos(Zones,Folder,MagLim):
	''' Function to read hipparcos catalogue (old)
	
	INPUTS:
	catalog_path 		-> Directory of catalog file								(str)
	max_mangitude 		-> maximum magnitude of object to save						(float)
	min_magnitude 		-> minimum magnitude of object to save						(float)
	max_variability 	-> maximum variability of object to save					(float)
	
	RETURNS:
	
	Catalogue			-> SkyCoord object with Ra and Dec catalogue				(SkyCoord)
	mags				-> array with corresponding magnitude						(np.array)
	
	FILES CREATED:
	NONE
	'''
	
	Folder += 'hipparcos.fits'
	catalog = Table.read(Folder)
	mask = np.isfinite(catalog['ra']) & np.isfinite(catalog['dec'])
	catalog = catalog[mask]

	if MagLim is not None:
		catalog = catalog[catalog['v_mag'] <= MagLim]

	return SkyCoord(ra = Longitude(catalog['ra'], unit = u.deg), dec = np.array(catalog['dec'])*u.deg,frame = 'cirs'), catalog['v_mag']

def ReadDucati(Zones,Folder,MagLim):
	''' Function to read Duccati catalog from Vizier: https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=II/237/colors&-out.max=50&-out.form=HTML%20Table&-out.add=_r&-out.add=_RAJ,_DEJ&-sort=_r&-oc.form=sexa
	
	INPUTS:
	catalog_path 		-> Directory of catalog file								(str)
	max_mangitude 		-> maximum magnitude of object to save						(float)
	min_magnitude 		-> minimum magnitude of object to save						(float)
	
	RETURNS:
	Catalogues			-> SkyCoord object with Ra and Dec catalogue				(SkyCoord)
	mags				-> array with corresponding magnitude						(np.array)
	
	FILES CREATED:
	NONE
	'''
	Folder += 'ducati2002.fit'
	catalog = {}
	HDU = fits.open(Folder)
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

def ReadTychoCatalogue(Zones,Folder,MagLim):
	''' 
	Function to read Tycho catalogue in txt form taken from http://archive.eso.org/ASTROM/TYC-2/data/
	
	INPUTS:
	catalog_path 		-> Directory of catalog file								(str)
	max_mangitude		-> maximum magnitude of object to save						(float)
	min_magnitude 		-> minimum magnitude of object to save						(float)
	
	RETURNS:
	Catalogue			-> SkyCoord object with Ra and Dec catalogue				(SkyCoord)
	mags				-> array with corresponding magnitude						(np.array)
	
	FILES CREATED:
	NONE
	'''
	
	Folder += 'tycho2.dat'
	Ra = []
	Dec = []
	Mag = {'BTMag':[],'VTMag':[]}
	with open(Folder,'r') as text:
		while True:
			line = text.readline()
			if line == '':
				break
			else:
				BTMag = line[110:116].strip()
				VTMag = line[123:129].strip()
				
				if BTMag != '' and VTMag != '':
					Ra.append(float(line[152:164].strip()))
					Dec.append(float(line[165:177].strip()))
					Mag['BTMag'].append(float(BTMag))
					Mag['VTMag'].append(float(VTMag))
	# ~ catalog['BTMag'] = np.array(catalog['BTMag'])
	# ~ catalog['VTMag'] = np.array(catalog['VTMag'])
	return SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg,frame = 'cirs'), Mag
	
ReadCatalogue = {
    'Tycho2': ReadTychoCatalogue,
    'UCAC4': ReadUCAC4,
    'UCAC5': ReadUCAC5,
    'Hipparchos': ReadHipparcos,
    'Ducati': ReadDucati, 
}
	
##################################################################################################################################################################

# Coordinate transformation
			      		      		
def XYToRaDec(FileWc,FileXY,FileOutput):
	''' 
	Function created to transform XY positions of an image into RaDec coordinates
	using WCS created by previous astrometry run. IT uses an astrometry function wcs-xy2rd.
	
	INPUTS:
	FileWc 				-> WCS file directory										(wcs)
	FileXY 				-> XY coordinates file given by sextractor					(fits)
	FileOutput 			-> output file to write with Ra Dec of corresponding XY		(str)
	
	OUTPUTS:
	result				-> result of wcs-xy2rd function								(str)
	status				-> status of wcs-xy2rd function								(str)
	
	FILES CREATED:
	NONE
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
	''' 
	Function created to transform RaDec positions of an image into XY coordinates
	using WCS created by previous astrometry run. IT uses an astrometry function wcs-xy2rd.
	
	INPUTS:
	FileWc 				-> WCS file directory										(wcs)
	FileRaDec 			-> RaDec coordinates file 									(fits)
	FileOutput 			-> output file to write with Ra Dec of corresponding XY		(str)
	
	OUTPUTS:
	result				-> result of wcs-rd2xy function								(str)
	status				-> status of wcs-rd2xy function								(str)
	
	FILES CREATED:
	NONE
	'''
	cmd = ['wcs-rd2xy', 
		'-w',FileWc, 									# WCS file directory
		'-i', FileRaDec,									# Sextractor catalogue with XY positions file directory
		'-o', FileOutput]								# Output file directory
		
	cmd_str = ' '.join(cmd)								# Compacts the cmd dictionary created so it can be run through terminal
	
	result = subprocess.run(cmd, capture_output=True)
	status = result.returncode  						# Returns result, status = 0 if good and 255 if bad
	return result, status

def HMS2deg(RA):
	''' 
	Function created to transform HMS to degrees
	
	INPUTS:
	Ra 					-> Ra value in HMS												(str)
	
	OUTPUTS:
	degrees				-> degrees value corresponding to input							(float)
	
	FILES CREATED:
	NONE
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
	Function created to transform DMS to degrees
	
	INPUTS:
	Dec 				-> Dec value in DMS												(str)
	
	OUTPUTS:
	degrees				-> degrees value corresponding to input							(float)
	
	FILES CREATED:
	NONE
	'''
	
	degrees = []
	for value in DEC:
		D = float(value[1:3])
		M = float(value[4:6])
		S = float(value[7:10])
		degree = D + M/60 + S/3600
		degrees.append(degree)
		
	return degrees

def str2bool(BoolStr):
	''' 
	Function created to transform string to boolean
	
	INPUTS:
	BoolStr 			-> Boolean in string											(str)
	
	OUTPUTS:
	BoolStr				-> Boolean in bool												(bool)
	
	FILES CREATED:
	NONE
	'''
	if BoolStr == "True":
		return True
	if BoolStr == "False":
		return False
	else:
		return BoolStr

def AngleDifferencePoints(Coord1,Coord2):
	Ra1  = np.radians(Coord1[0])
	Dec1 = np.radians(Coord1[1])
	Ra2  = np.radians(Coord2[0])
	Dec2 = np.radians(Coord2[1])
	
	difference = np.arccos(np.sin(Dec1)*np.sin(Dec2)+np.cos(Dec1)*np.cos(Dec2)*np.cos(Ra1-Ra2))
	return np.degrees(difference)
	
def AngleDifferenceLists(L1,L2):
	''' Calculates Angular difference between two sets of data Keys: 'RA' and 'DEC'
	L -> {RA:[Ra1,Ra2],DEC:[Dec1,Dec2]}
	
	
	'''
	if len(L1) != len(L2):
		print('Different size lists! Please use two lists of same size')
		return None
	
	Dists = []
	for i in range(len(L1['RA'])):
		Coord1 = [L1['RA'][i],L1['DEC'][i]]
		Coord2 = [L2['RA'][i],L2['DEC'][i]]
		Dists.append(AngleDifferencePoints(Coord1,Coord2))
	
	return np.array(Dists)

def decdeg2dms(dd):
    mult = -1 if dd < 0 else 1
    mnt,sec = divmod(abs(dd)*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return round(mult*deg), round(mult*mnt), float("{:.1f}".format(mult*sec))
    
##################################################################################################################################################################

# Masters Functions

def BestDark(Directory,time):
	''' 
	Function that checks to see whats the best dark to use for a given image in a directory of masters, checks for closest exposure time
	
	INPUTS:
	Directory 			-> MasterDarks directory 										(str)
	time 				-> Exposure time for image 										(float)
	
	RETURNS:
	Masters				-> List with Darks name											(str)
	
	FILES CREATED:
	NONE
	'''
	files = os.listdir(Directory)
	files = [B for B in files if B[6:10] == 'Dark' and int(B[10:13]) <= time]
	times = [B[11:14] for B in files]
	if times == []:
		return []
	else:
		return files[times.index(max(times))]

def CreateMasters(DayDirectory,config):
	''' 
	Function to create master files from available Flats and Darks. It saves with exposure and type of filter in name.
	
	INPUTS:
	DayDirectory 		-> Directory for folder with darks and flats files				(str)
	config 				-> Dictionary with configurations from Config File				(dict)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	MasterFiles			-> Creates MasterFlat and MasterDark for available files		(fits)
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
	''' 
	Function to create single master file from set of files.
	
	INPUTS:
	Files 				-> Array with directory for files to be used 					(list)								
	DayDirectory 		-> Directory for folder with darks and flats files				(str)
	Filter 				-> Filter/time used in files, to be used to save in name		(str)
	Type 				-> Type of master file (Dark,Flat)								(str)
	config 				-> Dictionary with configurations from Config File				(dict)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	Master				-> Creates Master file											(fits)
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
		Master.data = Master.data/np.max(Master.data)
		
	MasterName = "Master%s%s.fits" % (Type, Filter)
	
	Master.write(DayDirectory + '/Masters/' + MasterName)
		
def OrganizeDarks(DayDirectory):
	''' 
	Function that checks available DarkFiles:
	
	INPUTS:
	DayDirectory 		-> Directory for folder with darks and flats files				(str)
	
	RETURNS
	Dic 				-> Dictionary with exposure:files								(dict)
	
	FILES CREATED:
	NONE
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
	''' 
	Function that checks which masters are available
	
	INPUTS:
	DayDirectory 		-> Directory for folder with darks and flats files				(str)
	
	RETURNS:
	Masters 			-> {Flat: ['A','B','C']; Dark:['50','70','100','200']}			(dict)
	
	FILES CREATED:
	None
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
	

##################################################################################################################################################################

# Config Files reader
   
def ReadPipelineConfig(ficheiro):
	''' 
	Function to read pipeline config files. It skips over blank lines and lines which start with an #. 
	Returns in the form of a dictionary with properties as keys and True/False or numbers as values
	
	INPUTS:
	ficheiro 			-> Config file directory										(str)
	
	RETURNS:
	dic					-> Dictionary with config keys									(dict)
	
	FILES CREATED:
	NONE
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
       		
def ReadSextractorFileTxt(SexFile):
	''' 
	Function created to read catalogue created from sextractor in txt type. It outputs 
	a dictionary with the information related to each column as values of keys
	
	INPUTS:
	SexFile 			-> Directory of sextractor fits catalogue						(txt)
	
	RETURNS:
	config				-> Dictionary with information from sextractor catalog (txt)	(dict)
	
	FILES CREATED:
	NONE
	
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
	''' 
	Function created to read catalogue created from sextractor in fits type. It outputs 
	a dictionary with the information related to each column as values of keys
	
	INPUTS:
	SexFile 			-> Directory of sextractor fits catalogue						(fits)
	
	RETURNS:
	config				-> Dictionary with information from sextractor catalog (txt)	(dict)
	
	FILES CREATED:
	NONE
	
	'''
	
	
	config = {}
	Ficheiro = fits.open(SexFile)
	Dados = Ficheiro[-1]
	Colunas = Dados.columns
	for i in Colunas:
		Coluna = i.name
		config[Coluna] = Dados.data[Coluna]
		
		
	return config

def SeeFitsTables(filedir):
	hdul = fits.open(filedir)
	return hdul[1].data

def SeePSF(File):
	''' 
	Function that plots PSF from PSFex output file.
	
	INPUTS:
	File 				-> Directory of psf file										(fits)
	
	RETURNS:
	NONE
	
	FILES CREATED:
	NONE
	
	'''
	A = fits.open(File)
	PSF = A[1].data['PSF_MASK'][0,0,:,:]
	
	size = np.shape(PSF)[-1]
	x = np.linspace(0,size-1,size)
	y = np.linspace(0,size-1,size)
	X,Y = np.meshgrid(x,y)
	
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	ax.plot_wireframe(X,Y,PSF,color = 'black')
	plt.show()
	
			   
##################################################################################################################################################################

# File Changer

def Rem0FWHM(FileName):
	''' 
	Removes detections with FWHM = 0, before it is passed to psfex.
	
	INPUTS:
	FileName -> Name of file
	
	RETURNS:
	
	FILES CREATED:
	
	
	'''
	A = fits.open(FileName, mode = 'update')
	A[2].data = A[2].data[A[2].data['FWHM_IMAGE'] > 0]
	A.close()
  
def CreateRDLSFile(Ra,Dec,FileName):
	''' 
	Function that changes the data from the RLDS file to a set of (Ra,Dec) coordinates.
	
	INPUTS:
	Ra 					-> array with Ra coordinates									(list)
	Dec 				-> array with Dec coordinates									(list)
	FileName 			-> Name of rdls file											(rdls)
	
	RETURNS:
	NONE
	
	FILES CHANGED:
	FileName			-> Changes Ra,Dec information to match with Ra,Dec given		(rdls)
	'''
	with fits.open(FileName, mode = 'update') as A:
		Dados = A[1].data
		
		x = np.array(Ra,Dec,dtype = [('RA', '<f8'),('DEC','<f8')])
		Dados = Dados[Dados['RA'] > 180]
		Dados = np.append(Dados,x)

		A.close()
	
def RemoveBorders(File,P):
	''' 
	Function that removes detections too close to the border (P distance to border)
	
	INPUTS:
	
	File 				-> sextractor file with detections								(list)
	P 					-> Distance threshold											(rdls)
		
	RETURNS:
	NONE
	
	FILES CHANGED:
	File				-> Removes information about detections too close to border		(rdls)
	'''
	with fits.open(File, mode = 'update') as A:
		Dados = A[-1].data
		A[-1].data = Dados[((Dados['Y_IMAGE']>P) & (Dados['Y_IMAGE']<2048-P)) & ((Dados['X_IMAGE']>P) & (Dados['X_IMAGE']<2048-P))]
		A.close()
	
def WriteOutput_Txt(File,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM):
	'''
	Creates txt Output File with Information about detections:
	X,Y,(Ra,Dec) calculated, (Ra,Dec) catalogue, Mag calculated, Mag catalogue
	
	'''
	logging.info("Started Writing Output of image " + File)
	Base = File.split('.')[0]
	with open(DayDirectory + '/OutputFiles/' + Base,'a') as ficheiro:
		linha = f'{"X":20},{"Y":20},{"SNR":5},{"FWHM":10},{"Match":6},{"RaCalculated":20},{"DecCalculated":20},{"RaCatalogue":20},{"DecCatalogue":20},{"MagInst":19},{"MagReal":8}\n'
		ficheiro.write(linha)
		for i in range(len(Detections)):
			Xdetect = X[i]
			Ydetect = Y[i]
			RaCalculated = Detections.ra[i].value
			DecCalculated = Detections.dec[i].value
			MagInst = DetectionsMags[i]
			Match = found[i]
			if Match:
				RaCatalogue = Catalogue.ra[i].value
				DecCatalogue = Catalogue.dec[i].value
				MagReal = CatalogueMags[i]
			else:
				RaCatalogue = '-'
				DecCatalogue = '-'
				MagReal = 0
			
			
			
			linha = f'{Xdetect:20},{Ydetect:20},{SNR[i]:5},{FWHM[i]:10},{Match:6},{RaCalculated:20},{DecCalculated:20},{RaCatalogue:20},{DecCatalogue:20},{MagInst:19},{MagReal:8}\n'
			ficheiro.write(linha)
	
	logging.info("Finished Writing Output of image " + File)

def ReadOutput_Txt(OutputFile):
	'''
	Reads an output file created by WriteOutput and returns a dictionary with:
	{Variable1:[A1,A2,A3],Variable2:[B1,B2,B3]} (strings)
 	
	'''
	Dict = {}
	with open(OutputFile) as File:
		FirstLine = File.readline()
		FirstLine = FirstLine.replace(' ','')
		Variables = FirstLine.strip().split(',')
		for i in Variables:
			Dict[i] = np.array([])
		
		while True:
			Line = File.readline()
			if Line == '':
				break
			Line = Line.replace(' ','')
			Data = Line.strip().split(',')
			for i in range(len(Data)):
				Dict[Variables[i]] = np.append(Dict[Variables[i]],Data[i])
			
	
	return Dict
		
def ReadTLES_txt(Links,Files):
	if len(Links) != len(Files):
		print('The length of Links array is different from Files array!')
		return None
	
	Sucess = [False]*len(Links)
	
	for i in range(len(Links)):
		response = requests.get(Links[i])
		if response.status_code == 200:
			with open(Files[i],'w') as Teste:
				Teste.write(response.text)
				Sucess[i] = True
		
	return Sucess
	
def ReadTLES_json(Links):
	
	'''
	ReadTLES(Links) function uploads tle information from a url and saves it in a .json file (a list of dicts).
	
	INPUTS:
	Links		->   Celestrak links of Satellites TLEs                               (lst)
	
	FILES CREATED:
	GEO_tles.json
	'''
	
	print('[+] Downloading TLE data...')
	tle_jsonfile=[]
	
	for i in range(len(Links)):
		response = requests.get(Links[i])

	if response.status_code == 200: #The request was successful and the server responded with the data you were requesting
		tmp_dict = {}
		tle_data=response.text.split('\n')
		count=0
		
		for i in range(0,len(tle_data)-1):
			
			if tle_data[i][0] == '1':
				tmp_dict['tle_1'] = tle_data[i].strip()
			if tle_data[i][0] == '2':
				tmp_dict['tle_2'] = tle_data[i].strip()
			if tle_data[i][0] != '1' and tle_data[i][0] != '2':
				tmp_dict['satellite_name'] = tle_data[i].strip()

			if "tle_1" in tmp_dict and "tle_2" in tmp_dict and "satellite_name" in tmp_dict:
				tle_jsonfile.append(tmp_dict)
				tmp_dict = {}
			else:
				pass
				
		with open('GEO_tles.json', 'w') as f:
			json.dump(tle_jsonfile, f, indent=3)
			print('[+] TLE data downloaded in GEO_tles.json')
	else:
		print('[+] ERROR in downloading TLE data')

def WriteSextractorFile_Json(DayDirectory):
	'''
	Function created to read catalogue created from sextractor in json type. It outputs 
	a list of dictionaries with the information related to each column as values of keys
	
	INPUTS:
	DayDirectory		->   directory of the present day									(str)
	
	FILES CREATED:
	SExtractorOutput.json
	'''
	
	all_files=files(DayDirectory)
	
	SEx_directory = DayDirectory + '/SextractorCatalogue/'
	
	lst_jsonfile=[]
	
	for ficheiro in all_files:
		
		filename= SEx_directory + f"{ficheiro}.cat"

		A=fits.open(filename)
	
		for i in A[1].data:
			tmp_dict = {}
			tmp_dict['image']= str(ficheiro[-2])+str(ficheiro[-1])
			tmp_dict['number'] = str(i[0])
			tmp_dict['Flux_isocor'] = str(i[1])
			tmp_dict['Mag_isocor'] = str(i[2])
			tmp_dict['MagErr_isocor'] = str(i[3])
			tmp_dict['Flux_aper'] = str(i[4])
			tmp_dict['FluxErr_aper'] = str(i[5])
			tmp_dict['x_image'] = str(i[6])
			tmp_dict['y_image'] = str(i[7])
			tmp_dict['Flags'] = str(i[8])
			tmp_dict['FWHM_Image'] = str(i[9])
			tmp_dict['Elongation'] = str(i[10])
			tmp_dict['Flux_Radius'] = str(i[11])
			tmp_dict['SNR_WIN'] = str(i[12])
			lst_jsonfile.append(tmp_dict)
	
	with open('SExtractorOutput.json', 'w') as f:
		json.dump(lst_jsonfile, f, indent=3)

def ReadFilesHeader(DayDirectory,header,File):
	
	'''
	Function created to return a list of dicts with the header information (time,date,scale,mean fwhm,dimensions) of all processed 
	files. The size of the list is associated with the number os processed files
	
	INPUTS:
	DayDirectory	->   directory of the present day									(str)
	
	RETURNS:
	header	->   header information of all processed images                  		(list of dicts)
	
	FILES CREATED:
	None
	'''
	

	Header_information={}

	Header_information['Image'] =File[-6]+File[-5] #Stores the last two strings of file name, that is associated to the image number, so you know to which image the header information is associated to.
	Header_information['Dimensions'] =header['NAXIS1']
	Header_information['Scale']=float(header['SCALE'])
	Header_information['Alt']=float(header['ALTUR'])
	Header_information['Az']=float(header['AZIMUTH'])
	
	#In current file headers, date and time information is stored in the format: 21/04/23 and 13:09:12.2 (an example)
	#The code bellow splits the caracters from : and / and returns the strings that begin with 0 with only the second 
	#element (04 --> 4, 09 --> 9).This is helpful when using datetime module in further analysis, that can't interpret 
	#month 04 or hour 09. It also adds the string 20 to the caracter asscociated to the year for the same reason (23 --> 2023).
		
	Day=header['Date'].split("/") 
	for i in Day: 
		if i[0]=='0':
			idx=Day.index(i)
			Day[idx]=i[1]

	Time_UT=header['FULLTIME'].split(":")

	for i in Time_UT:
		if i[0]=='0':
			idx=Time_UT.index(i)
			Time_UT[idx]=i[1]
	
	Header_information['Year']=int('20'+Day[2])
	Header_information['Month']=int(Day[1])
	Header_information['Day']=int(Day[0])
	Header_information['Hour']=int(Time_UT[0])
		
		
	# This part of the algoritm adds half of the exposure time to the time when the image was taken.
	exp_time=float(header['EXPTIME'])
	minute=int(Time_UT[1])
	second=float(Time_UT[2])+exp_time/2
		
	if second >= 60:
		minute=minute+1
		if minute >= 60:
			Header_information['Minute']=int(minute-60) 
		else:
			Header_information['Minute']=int(minute)
			
		Header_information['Second']=int(second-60) 
		Header_information['Microsecond']=int(abs(int(second)-second)*10**6)
	
	else:
		Header_information['Minute']=minute
		Header_information['Second']=int(second)
		Header_information['Microsecond']=int(abs(int(second)-second)*10**6)
		

	return data
	
def WriteOutput_Json(File,DayDirectory,header,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM):
	 
	'''
	Function that creates a json file with the information about detections from Astrometry:
	X,Y,(Ra,Dec) calculated, (Ra,Dec) catalogue, Mag calculated, Mag catalogue, Altitude and Azimuthe. 
	
	INPUTS:
	File 				-> File name without path and ending 						(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	Detections			-> SkyCoord Object with image detections (Ra/Dec)			(SkyCoord)
	DetectionsMags		-> Instrumental magnitudes corresponding to Detections		(np.array)
	Catalogue			-> SkyCoord Object with catalogue matches (Ra/Dec)			(SkyCoord)
	CatalogueMags		-> G magnitudes corresponding to Detections					(np.array)
	found				-> Boolean array with matches								(np.array)
	
	OUTPUTS:
	NONE
	
	FILES CREATED:
	{File}.json --> the number of files created is related to the number of files
	
	'''
	logging.info("Started Writing Output of image " + File)

	items_2_store=[]

	#The location where the images were taken is necessary to convert RA,DEC coordinates to Alt,Az coordinates
	location=EarthLocation(
			lat=float(config['Latitude'])*u.deg,
			lon=float(config['Longitude'])*u.deg,
			height=float(config['Height'])*u.m)
	
	#Header Information of proessed files

	header_data=ReadFilesHeader(DayDirectory,header,File)

	#The algoritm only continues when the image of which the header information is associated to is the same as the image of current Base
	time = datetime.datetime(header_data['Year'],header_data['Month'],header_data['Day'],header_data['Hour'],header_data['Minute'],header_data['Second'],header_data['Microsecond'])

	for i in range(len(Detections)):
			
			Match = found[i]
			items={}
			Base = File.split('.')[0]
			items['image']=File.split('.')[0][-2:]
			items['Xdetect']=float(X[i])
			items['Ydetect']= float(Y[i])
			items['FWHM']= float(FWHM[i])
			items['SNR']= float(SNR[i])
			items['RaCalculated'] = float(Detections.ra[i].value)
			items['DecCalculated'] = float( Detections.dec[i].value)
			items['MagInst '] = str(DetectionsMags[i])
			items['Match']= int(Match)
			
			if Match:
				items['RaCatalogue'] = str(Catalogue.ra[i].value)
				items['DecCatalogue'] = str(Catalogue.dec[i].value)
				items['MagReal'] =str(CatalogueMags[i])
				
			if not Match:
				items['MagReal'] ='0'

			

			dec=(Detections.dec[i].value)*u.deg
			ra=(Detections.ra[i].value)*u.deg
			stars = SkyCoord(ra, dec, frame='icrs')

			#Converts Ra,Dec to Alt,Az
			stars_altaz = stars.transform_to(AltAz(obstime=time, location=location))
			sky_coords = [stars_altaz.alt.degree,stars_altaz.az.degree]

			items['Alt']=float(sky_coords[0])
			items['Az']=float(sky_coords[1])

			items_2_store.append(items)
			
	with open(f"{DayDirectory}/OutputFiles/{File.split('.')[0]}.json", 'w') as f: 
		json.dump(items_2_store, f, indent=3)
						
	logging.info("Finished Writing Output of image " + File)

##################################################################################################################################################################

# Specific Softwares

def sextractorOld(File,PipelineDirectory,DayDirectory,SexConfig,cattype,VIGNET=None,fit_psf=False):
	''' 
	Function to run sextractor on an image and output catalogue file
	
	INPUTS:
	Base 				-> image directory to give to sextractor						(str)
	PipelineDirectory 	-> directory of the Pipeline folder								(str)
	DayDirectory 		-> directory of the present day									(str)
	SexConfig 			-> Configuration catalog used for Sextractor					(str)
	cattype 			-> Type of output catalog										(str)
	VINGET 				-> Size of array to save from each image, for PSFEX				(int)
	fit_psf 			-> Boolean to check wheter to use psf for photometry			(int)
	
	RETURNS:
	result				-> result sextractor software									(str)
	
	FILES CREATED:
	Writes a sextractor catalog at the day Sextractor folder with the same name as the original image
	
	
	'''
	Base = File.split('.')[0]
	FileBase = f"{DayDirectory}/SextractorCatalogue/{Base}"
	SexFile = f"{DayDirectory}/ReducedImages/{File}"
	SexDir = f"{PipelineDirectory}/Files/Sextractor/" 				# The sextractor directory is inside the pipeline folder
	
	if fit_psf:
		FileBase = FileBase + '_psf'
	
	if cattype == 'FITS_1.0':
		catout = FileBase +'.cat'
	elif cattype == 'FITS_LDAC':
		catout = FileBase +'LDAC.cat'
	else:
		raise Exception('Not an acceptable catalog type')
	
	
	cmd = {}		
	cmd['source-extractor'] = SexFile								# Chooses the image directory as the image to analyse			
	cmd['-c'] = SexDir + SexConfig									# Chooses the sextractor catalogue as the default in the pipeline folder
	cmd['-CATALOG_NAME'] =  catout									# Chooses the catalogue output as the catout given in the function
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	cmd['-CATALOG_TYPE'] = cattype
	cmd['-PARAMETERS_NAME'] = SexDir + 'default.param'				# Chooses the needed files from the sextractor folder

	
	if fit_psf:
		cmd['-PSF_NAME'] = psf_out = FileBase.replace('_psf', 'LDAC_psf.fits')
		with open(SexDir + 'default.param', 'a') as myfile:
			myfile.write('MAG_PSF\n')
		
	if VIGNET != None:
		with open(SexDir + 'default.param', 'a') as myfile:
			myfile.write('VIGNET{}\n'.format((VIGNET,VIGNET)))
	
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items()))) # Compacts the cmd dictionary created so it can be run through terminal
	result = subprocess.run(list_cmd, capture_output=True) 			  # Runs the program in terminal and outputs the result and status (0 = good, 255 = bad)
	status = result.returncode
	
	if VIGNET != None:
		with open(SexDir + 'default.param', '+r') as myfile:
			lines = myfile.readlines()
			myfile.seek(0)
			myfile.truncate()
			myfile.writelines(lines[:-1])
	
	if fit_psf:
		with open(SexDir + 'default.param', '+r') as myfile:
			lines = myfile.readlines()
			myfile.seek(0)
			myfile.truncate()
			myfile.writelines(lines[:-1])
	
			
	return result

def sextractor(File,PipelineDirectory,DayDirectory,SexConfig,n):
	''' 
	Function to run sextractor on an image and output catalogue file
	
	INPUTS:
	File 				-> image name to give to sextractor								(str)
	PipelineDirectory 	-> directory of the Pipeline folder								(str)
	DayDirectory 		-> directory of the present day									(str)
	SexConfig 			-> Configuration catalog used for Sextractor					(str)
	n					-> 1 = default.param, 2 = default1.param, 3 = default2.param	(int)
	
	RETURNS:
	result				-> result sextractor software									(str)
	
	FILES CREATED:
	Writes a sextractor catalog at the day Sextractor folder with the same name as the original image
	
	
	'''
	Base = File.split('.')[0]
	FileBase = f"{DayDirectory}/SextractorCatalogue/{Base}"
	SexFile = f"{DayDirectory}/ReducedImages/{File}"
	SexDir = f"{PipelineDirectory}/Files/Sextractor/" 				# The sextractor directory is inside the pipeline folder
	
	if n == 1:
		Params = 'default.param'
		catout = FileBase +'.cat'
		cattype = 'FITS_1.0'
		
	elif n == 2:
		Params = 'default1.param'
		catout = FileBase +'LDAC.cat'
		cattype = 'FITS_LDAC'
		
	elif n == 3:
		Params = 'default2.param'
		FileBase = FileBase + '_psf'
		catout = FileBase +'LDAC.cat'
		cattype = 'FITS_LDAC'
	
	cmd = {}		
	cmd['source-extractor'] = SexFile								# Chooses the image directory as the image to analyse			
	cmd['-c'] = SexDir + SexConfig									# Chooses the sextractor catalogue as the default in the pipeline folder
	cmd['-CATALOG_NAME'] =  catout									# Chooses the catalogue output as the catout given in the function
	cmd['-FILTER_NAME'] = SexDir + 'default.conv'
	cmd['-WEIGHT_IMAGE'] = SexDir + 'weight.fits'
	cmd['-CATALOG_TYPE'] = cattype
	cmd['-PARAMETERS_NAME'] = SexDir + Params						# Chooses the needed files from the sextractor folder
	if n == 3:
		cmd['-PSF_NAME'] = FileBase.replace('_psf', 'LDAC_psf.fits')
	
	list_cmd = list(itertools.chain.from_iterable(list(cmd.items()))) # Compacts the cmd dictionary created so it can be run through terminal
	result = subprocess.run(list_cmd, capture_output=True) 			  # Runs the program in terminal and outputs the result and status (0 = good, 255 = bad)
	
			
	return result
		    		
def astrometry(DirCatalogo,header,DirAstro,config):
	''' 
	Function to run astrometry.net on an image and output created files
	
	INPUTS:
	DirCatalogo 		-> directory of the sextractor catalogue						(str)
	Header 				-> Image header													(CCData.header)
	DirAstro 			-> Directory for output files									(str)
	config 				-> Config file for the chosen telescope							(dict)
	
	RETURNS:
	result				-> result astrometry software									(?)
	status				-> status of astrometry software								(?)
	
	FILES CREATED:
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
	
def psfex(catin,catout,file_config,VIGNET):
	''' 
	Function to run astrometry.net on an image and output created files
	
	INPUTS:
	Base 				-> image directory to give to sextractor						(str)
	catin 				-> FITS_LDAC sextractor catalog									(str)
	catout 				-> Path to output file											(str)
	file_config 		-> PSFEX configuration file Path								(str)
	VIGNET 				-> VIGNET value used for sextractor catalog						(str)
	
	RETURNS:
	result				-> result psfex software										(str)
	status				-> status of psfex software										(str)
	
	FILES CREATED:
	Writes a . psf file
	'''
	cmd = ['psfex', catin, '-c', file_config,'-OUTCAT_NAME',catout, '-PSF_SIZE','{},{}'.format(VIGNET,VIGNET)] 
	cmd_str = ' '.join(cmd) 						  # Compacts the cmd dictionary created so it can be run through terminal
	result = subprocess.run(cmd, capture_output=True) # Runs psfex with given parameters
	status = result.returncode						  # Returns result, status = 0 if good and 255 if bad
	return result, status
   
def find_matching_stars(catalog_stars, image_stars, max_sep=0.5 * u.deg):
	''' 
	Function to find the closest match between two sets of positions. 
	This is used to find the closest match between found objects and stars in catalogue.
	
	INPUTS:
	Base 				-> image directory to give to sextractor						(str)
	catalog_stars 		-> SkyCoord object with catalogue stars coordinates				(SkyCoord)
	image_stars 		-> SkyCoord object with coordinates of objects detected 		(SkuCoord)
	max_sep 			-> maximum separation to be accepted for closest match 			(u.deg)
	
	RETURNS: 
	idx					-> index of image stars that were found							(list)
	mask				-> True/False matrix matching the image_stars array				(list)
	
	FILES CREATED:
	NONE
	'''
	
	
	if len(image_stars) == 0:
		return np.array([], dtype=int), np.zeros(len(catalog_stars), dtype=bool)

	idx, d2d, d3d = match_coordinates_sky(catalog_stars, image_stars,storekdtree = True)
	mask = d2d < max_sep

	return idx, mask


##################################################################################################################################################################

# ReductionSteps

def DarkCorrection(DayDirectory,Data,File,exposure,logging):
	''' 
	Function that corrects image for dark current.
	
	INPUTS:
	DayDirectory		-> directory of the present day									(str)
	File		 		-> SkyCoord object with catalogue stars coordinates				(CCDData)
	exposure	 		-> Image exposure time									 		(float)
	logging 			-> logger used in pipeline							 			(logger)
	
	RETURNS: 
	Data				-> Data corrected with dark										(CCDData)
	
	FILES CREATED:
	NONE
	'''
	logging.debug("Starting Dark correction of file " + File)						# Logger information
	Master = BestDark(DayDirectory + '/Masters',exposure)									# Sees what is the best masterdark to apply
	if Master == []:
		logging.warning('There are no Darks available for exposure ' + exposure)				# IF there is none, it skips over dark correction
		Data.header['DARK'] = 'False'													# Writes in header
	else:
		Master = DayDirectory + '/Masters/' + Master										# Gets master directory
		MasterDark = CCDData.read(Master) 												# Reads and then applies master
		Data = ccdproc.ccd_process(Data, dark_frame=MasterDark,dark_exposure = float(MasterDark.header["EXPOSURE"].replace(',','.')) * u.s,data_exposure = float(Data.header["EXPOSURE"].replace(',','.')) *u.s)
		Data.header['DARK'] = 'True'													# Writes in header
		logging.debug("Finished Dark correction of file " + File)					# Logger information

	return Data

def FlatCorrection(DayDirectory,Data,File,Filter,Masters,logging):
	''' 
	Function that performs flat correction on image data.
	
	INPUTS:
	DayDirectory		-> directory of the present day									(str)
	File		 		-> SkyCoord object with catalogue stars coordinates				(CCDData)
	Filter		 		-> Filter used in image									 		(float)
	Masters		 		-> Image exposure time									 		(float)
	logging 			-> logger used in pipeline							 			(logger)
	
	RETURNS: 
	Data				-> Data corrected with dark										(CCDData)
	
	FILES CREATED:
	NONE
	'''
	
	logging.debug("Starting Fits correction of file " + File)						# Logger information
	if Filter not in Masters['Flat']:
		logging.warning('There are no Flats available of type ' + Filter)				# Logger information if there is no flat
		Data.header['FLAT'] = 'False'													# Writes in header
	else:
		Master = DayDirectory + '/Masters/MasterFlat' + Filter + '.fits'					# Gets master directory
		MasterFlat = CCDData.read(Master)												# Reads master
		Data = ccdproc.ccd_process(Data, master_flat = MasterFlat)						# Applies Master
		Data.header['FLAT'] = 'True'													# Writes in header
		logging.debug("Finished Fits correction of file " + File)					# Logger information
	
	return Data
	
def BackgroundCorrection(File,Data,logging):
	''' 
	Function that removes background from image
	
	INPUTS:
	File		 		-> File name													(str)
	Data	 			-> Image data											 		(CCDAata)
	logging 			-> logger used in pipeline							 			(logger)
	
	RETURNS: 
	Data				-> Data corrected												(CCDData)
	
	FILES CREATED:
	NONE
	'''
	logging.debug("Eliminating Background of " + File)								# Logger information
	bkg = MeanBackground(SigmaClip(sigma=3.))
	bkg_value = bkg.calc_background(Data)												# Calculates background form image
	Data = Data.subtract(bkg_value*np.ones(np.shape(Data))*u.adu,propagate_uncertainties=True,handle_meta='first_found') # Removes background from image
	logging.debug("Background removed from " + File)								# Logger information
	Data.header['BG'] = 'True'															# Writes in header	
	return Data

def CrosstalkCorrection(File,Data,logging,config):
	''' 
	Function that removes crosstalk effect from image
	
	INPUTS:
	File		 		-> File name													(str)
	Data	 			-> Image data											 		(CCDAata)
	logging 			-> logger used in pipeline							 			(logger)
	config	 			-> configuration file for telescope					 			(dict)
	
	RETURNS: 
	Data				-> Data corrected												(CCDData)
	
	FILES CREATED:
	NONE
	'''
	logging.debug("Starting CrossTalk correction of " + File)						# Logger information
	alfa1 = float(config['alfa1'])																# Parameters for correction
	alfa2 = float(config['alfa2'])
	convolution_kernel = np.array([[0, alfa1, 0], 										# Deconvolution kernel for cross-talk correction
								   [alfa1, 1-4*alfa1, alfa1], 
								   [0, alfa1, 0]])
	Data.data = cv2.filter2D(np.array(Data.data), -1, convolution_kernel)				# Applies kernel to image
	logging.debug("Crosstalk corrected from " + File)
	Data.header['CT'] = 'True'														# Writes in header
	return Data															

def CosmicRayCorrection(File,Data,logging):
	''' 
	Function that removes cosmicrays from image
	
	INPUTS:
	File		 		-> File name													(str)
	Data	 			-> Image data											 		(CCDAata)
	logging 			-> logger used in pipeline							 			(logger)
	
	RETURNS: 
	Data				-> Data corrected												(CCDData)
	
	FILES CREATED:
	NONE
	'''
	logging.debug("Starting CosmicRays Correction of " + File)						# Logger information
	crmask, Data.data = astroscrappy.detect_cosmics(Data.data, inmask = None, sepmed = False, cleantype='medmask', fsmode = 'median', verbose=False)  #deteta os raios cósmicos numa matriz numpy
	logging.debug("CosmicRays removed from " + File)								# Logger information
	Data.header['CR'] = 'True'															# Writes in header
	
	return Data

###########################################################
# Drawings

def ImagePlot(Images):		
	for i in range(len(Images)):
		Image = Images[i]
		Data = CCDData.read(f"{DayDirectory}/ReducedImages/{Image}.fit", unit = 'adu')
		raw = Data.data
		
		OutputFile = f"{DayDirectory}/OutputFiles/{Image}"
		Dados = SPADEFunctions.ReadOutput(OutputFile)
		X,Y,found = Dados['X'],Dados['Y'],Dados['Match']
		
		fig,axes = plt.subplots(2,2)
		ax = axes.ravel()
		vmin,vmax = None,None
		fig.suptitle(f'Image: {Imagens[i]} ', fontsize=16)
		ax[0].set_title('Imagem Inicial')
		ax[0].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,3),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].set_title('Detecoes')
		ax[1].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,3),vmax=vmax or np.nanpercentile(raw, 99))
		ax[1].set_xlabel(len(X))
		ax[2].set_title('Encontradas')
		ax[2].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,3),vmax=vmax or np.nanpercentile(raw, 99))
		ax[2].set_xlabel(sum(found))
		ax[3].set_title('Em Falta')
		ax[3].imshow(raw,cmap='gray',vmin=vmin or np.nanpercentile(raw,3),vmax=vmax or np.nanpercentile(raw, 99))
		ax[3].set_xlabel(len(X)-sum(found))
		size = 3
		for i in range(len(X)):
			c1 = plt.Circle((float(X[i]),float(Y[i])),size, color = 'red', linewidth=0.4,fill=True)
			ax[1].add_patch(c1)
			if int(found[i]):
				c1 = plt.Circle((float(X[i]),float(Y[i])),size, color = 'red', linewidth=0.4,fill=True)
				ax[2].add_patch(c1)
			else:
				c1 = plt.Circle((float(X[i]),float(Y[i])),size, color = 'red', linewidth=0.4,fill=True)
				ax[3].add_patch(c1)


	plt.show()		
		

##################################################################################################################################################################

# Main Code	 				

def ImageReduction(File,PipelineDirectory,DayDirectory,Masters,config,logging):
	'''
	Function that reduces each Image, correcting for Dark,Gain,Flat,Background,Crosstalk and CosmicRays according to the config file.
	
	INPUTS:
	File 				-> File name with path										(str)
	PipelineDirectory 	-> Directory of Pipeline folder								(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	logging				-> logger used in Pipeline									(logger)
	
	RETURNS:
	Image Header
	
	FILES CREATED (ReducedImages):
	File	 			-> Corrected image file										(fits)
	'''
	logging.info("Starting Reduction on image " + File)										# Logger information
	Data = CCDData.read(f"{DayDirectory}/ProcessedImages/{File}" , unit = 'adu') 			# Reads Image
	Filter = Data.header['FILTER']															# Gets filter from header
	
	# ~ exposure = float(Data.header['EXPTIME'].replace(',','.'))							# Bruno
	exposure = float(Data.header['EXPTIME'])												# Defesa
	
	if config["DARK"]:																		# If DARK is True does dark correction
		Data = DarkCorrection(DayDirectory,Data,File,exposure,logging)
	else:	
		Data.header['DARK'] = 'False'														# Writes in header
		logging.debug("Skipped Dark correction of file " + File)
	
	if config['GAIN']:
		Data.data = Data.data/float(Data.header['GAIN'])
		Data.header['GAIN'] = '1'
		logging.debug("Finished Gain correction of file " + File)
	else:
		Data.header['GAIN'] = 'False'														# Writes in header
		logging.debug("Skipped Gain correction of file " + File)
		
	if config["FLAT"]:
		Data = FlatCorrection(DayDirectory,Data,File,Filter,Masters,logging)
	else:
		Data.header['FLAT'] = 'False'														# Writes in header
		logging.debug("Skipped Flat correction of file " + File)
	
	if config["Background"]:
		Data = BackgroundCorrection(File,Data,logging)
	else:
		Data.header['BG'] = 'False'															# Writes in header
		logging.debug("Skipped Background correction of file " + File)															
	
	if config["CrossTalk"]:																	# Stars crostalk correction if enabled
		Data = CrosstalkCorrection(File,Data,logging,config)
	else:
		Data.header['CT'] = 'False'															# Writes in header
		logging.debug("Skipped Crosstalk correction of file " + File)
	
	if config["CosmicRays"]:
		CosmicRayCorrection(File,Data,logging)
	else:
		Data.header['CR'] = 'False'															# Writes in header
		logging.debug("Skipped CosmicRay correction of file " + File)
	
	DirImage = f"{DayDirectory}/ReducedImages/{File}"	
	Data.write(DirImage, overwrite = True)
	logging.info("Finished Reduction of image " + File)
	return Data.header

def SextractorPart(File,PipelineDirectory,DayDirectory,config):
	'''
	Function that processes all sextractor Steps. 
	Firstly it runs sextractor twice, in order to produce a FITS_LDAC (for psfex) and a FITS_1.0 (for astrometry.net). Detections too close to the borders are removed.
	Secondly it runs PSFEX to obtain the average PSF in the image, after removing all detections with FWHM=0.
	Lastly, it runs sextractor again with the PSF model obtained, in order to improve photometry on sources.
	
	INPUTS:
	Base 				-> File name without path and ending 						(str)
	PipelineDirectory 	-> Directory of Pipeline folder								(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	
	RETURNS:
	No returns
	
	FILES CREATED (SextractorCatalogue):
	Base.cat 			-> Sextractor catalog with FITS_1.0 types					(fits)
	BaseLDAC.cat 		-> Sextractor catalog with FITS_LDAC types					(FITS_LDAC)
	BaseLDAC_psf.fits	-> PSF model obtained with PSFEX							(Psf file)
	Base_out.cat		-> Ouput file with PSFEX									(txt)
	Base_psfLDAC.cat	-> Sextractor catalog obtained using PSF model				(Final File)
	'''
	Base = File.split('.')[0]
	FileBase = f"{DayDirectory}/SextractorCatalogue/{Base}"
	VIGNET = int(config['VIGNET'])
	logging.info("Started Sextractor Pass of image " + File)
	
	# First Sextractor pass: creates 2 catalogs: one in Fits_1.0 format for use in astrometry.net and 1 in Fits_LDAC for use in PSFEX, then sources too close to the
	# border are elliminated. The threshold for the distance is the VIGNET parameter given in the configuration file
	result = sextractor(File,PipelineDirectory,DayDirectory,config['SexConfig'],1)
	if config['VERBOSE']:
		print(result)
	result = sextractor(File,PipelineDirectory,DayDirectory,config['SexConfig'],2)
	if config['VERBOSE']:
		print(result)
	RemoveBorders(FileBase +'LDAC.cat',int(config['VIGNET']))
	RemoveBorders(FileBase +'.cat',int(config['VIGNET']))	
	
	
	#Secondly, PSFEX is run with the sextractor catalog written before.
	CatIn 	= f"{FileBase}LDAC.cat"
	CatOut 	= f"{FileBase}_out.cat"
	ConfigFile =  f"{PipelineDirectory}/Files/Sextractor/{config['PSFexConfig']}"
	Rem0FWHM(CatIn)
	result, status = psfex(CatIn,CatOut,ConfigFile,VIGNET)
	
	psf_in = CatIn.replace('.cat', '.psf')
	psf_out = CatIn.replace('.cat', '_psf.fits')
	os.rename (psf_in, psf_out)
	
	
	#Lastly, Sextractor is run again with the PSF model, in order to improve photometry
	logging.info("Finished Sextractor Pass of " + File)
	resultado = sextractor(Base,PipelineDirectory,DayDirectory,config['SexConfig'],3)

def AstrometryPart(File,PipelineDirectory,DayDirectory,config,header):
	'''
	Function that runs astrometry and compares image detections with catalogue positions.
	Firstly it runs astrometry with sextractor catalogue created previously. Then, it transforms the detections XY coordinates to Ra/Dec system. 
	Secondly it reads portion of the designed catalogue to compare with image detections.
	Lastly, it runs sextractor again with the PSF model obtained, in order to improve photometry on sources.
	
	INPUTS:
	Base 				-> File name without path and ending 						(str)
	PipelineDirectory 	-> Directory of Pipeline folder								(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	
	RETURNS:
	No returns
	
	FILES CREATED (AstrometryFiles):
	Base.wcs 			-> WCS file created by astrometry							(fits)
	Base.rdls 			-> RDLS file created by astrometry							(FITS_LDAC)
	Base-indx.xyls		-> XY coordinates of stars used for plate-solving			(Psf file)
	Base				-> Ra/Dec coordinates of detections							(Final File)
	'''
	logging.info("Started Astrometry Pass of image " + File)
	Base = File.split('.')[0]
	result,status = astrometry(DayDirectory + '/SextractorCatalogue/' + Base + '.cat',header,DayDirectory + '/AstrometryFiles',config) 

	# ~ #Transforms from XY to RaDec with WCS information from astrometry
	FileWc = f"{DayDirectory}/AstrometryFiles/{Base}.wcs"															# Gets directory : WCS file created from astrometry 
	FileXY = f"{DayDirectory}/SextractorCatalogue/{Base}.cat"														# Gets directory : catalogue created from sextractor
	FileOutput = f"{DayDirectory}/AstrometryFiles/{Base}"																# Output file for calculation of RaDec coordinate of all stars
	result, status	= XYToRaDec(FileWc,FileXY,FileOutput)																# Calculates RaDec positions for all stars obtained from sextractor
	logging.info("Finished Astrometry Pass of image " + File)

def CatalogComparison(File,PipelineDirectory,DayDirectory,config,Filter):
	'''
	Function that compares image detections with the catalog.
	
	INPUTS:
	Base 				-> File name without path and ending 						(str)
	PipelineDirectory 	-> Directory of Pipeline folder								(str)
	DayDirectory 		-> Directory of present day									(str)
	config				-> Configuration Dictionary									(dic)
	
	RETURNS:
	Detections			-> SkyCoord Object with image detections (Ra/Dec)			(SkyCoord)
	DetectionsMags		-> Instrumental magnitudes corresponding to Detections		(np.array)
	Catalogue			-> SkyCoord Object with catalogue matches (Ra/Dec)			(SkyCoord)
	CatalogueMags		-> G magnitudes corresponding to Detections					(np.array)
	found				-> Boolean array with matches								(np.array)
	
	FILES CREATED (SextractorCatalogue):
	NONE
	'''
	logging.info("Started Catalogue Comparison of image " + File)
	Base = File.split('.')[0]
	FileOutput = f"{DayDirectory}/AstrometryFiles/{Base}"																# Output file for calculation of RaDec coordinate of all stars
	FileXY = f"{DayDirectory}/SextractorCatalogue/{Base}.cat"														# Gets directory : catalogue created from sextractor
	FileWc = f"{DayDirectory}/AstrometryFiles/{Base}.wcs"															# Gets directory : WCS file created from astrometry
	
	#Gets the information from the sextractor file and the calculated positions
	Posicoes = fits.open(FileOutput)												# Gets positions calculated
	Ra,Dec = Posicoes[-1].data['RA'],Posicoes[-1].data['DEC']
	SexData = ReadSextractorFileFit(FileXY)			# Reads sextractor file obtained in the beggining and saves desired values
	X,Y, DetectionsMags, FWHM, SNR = SexData['X_IMAGE'], SexData['Y_IMAGE'], SexData['MAG_ISOCOR'], SexData['FWHM_IMAGE'],  SexData['SNR_WIN']
	
	#Matches the objects found and the catalogue
	Detections = SkyCoord(ra = Longitude(Ra, unit = u.deg), dec = np.array(Dec)*u.deg, frame = 'cirs')  				# Creates SkyCoord object with RaDec positions obtained (needed for matching)
	WCS = fits.open(FileWc)
	scale = float(WCS[0].header['COMMENT'][-8].split()[1])
	max_sep = 0.8*float(config['max_sep'])*np.median(FWHM[FWHM!=0])*scale/3600
	
	catalog_path = PipelineDirectory + '/Files/Catalogues/'
	Catalogue,CatalogueMags = ReadCatalogue[config['Catalog']](ZonesUCAC(min(Dec),max(Dec)),catalog_path,float(config['MagLim']))
	idx,found = find_matching_stars(Detections, Catalogue, max_sep=max_sep * u.deg) 						# Starts matching between these stars and catalogue stars
	
	if Filter == 'none':
		Filter = 'G'
	CatalogueMags = np.array(CatalogueMags[Filter])[idx]
	Catalogue = Catalogue[idx]
	logging.info("Finished Catalogue Comparison of image " + File)
	return X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM
				
def ImageProcess(File,PipelineDirectory,DayDirectory,Masters,config,logging):
	header = ImageReduction(File,PipelineDirectory,DayDirectory,Masters,config,logging)
	SextractorPart(File,PipelineDirectory,DayDirectory,config)
	AstrometryPart(File,PipelineDirectory,DayDirectory,config,header)
	X, Y, Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM = CatalogComparison(File,PipelineDirectory,DayDirectory,config,header['FILTER'])
	# ~ WriteOutput_Txt(File,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)
	WriteOutput_Json(File,DayDirectory,header,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)

	#WriteOutput(File,DayDirectory,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)
	WriteAstrometryOutput_Json(File,DayDirectory,header,config,X,Y,Detections, DetectionsMags, Catalogue, CatalogueMags, found, SNR, FWHM)
	
def FindMatchJson(obj,ImageData,scale):  
	L1=[sum(obj['RaCalculated'])/len(obj['RaCalculated']),sum(obj['DecCalculated'])/len(obj['DecCalculated'])]
	for point in ImageData:  #fixates a point an image to compares to all other points from different images
		FWHM = max(float(point['FWHM']),float(obj['FWHM'])) 
		L2=[float(point['RaCalculated']),float(point['DecCalculated'])]
		distance=SPADEFunctions.AngleDifferencePoints(L1,L2)
	

		if float(distance) <= scale*FWHM/3600:
			return point
	
	return False
			
def CompareImagesJson(DayDirectory,Images,scale):
	if len(Images) == 1:
		raise Exception('Not enough images to compare!')
	dic = {}
	Image = Images[0]
	with open(f"{DayDirectory}/OutputFiles/{Image.split('.')[0]}.json", 'r') as stream:
		Data = json.load(stream)
	for count,obj in enumerate(Data):
		dic[count] = {}
		for i,key in enumerate(obj):
			dic[count][key] = []
			dic[count][key].append(float(obj[key]))
			
			
	for count, Image in enumerate(Images[1:-1]):
		with open(f"{DayDirectory}/OutputFiles/{Image.split('.')[0]}.json", 'r') as stream:
			Data = json.load(stream)
		for i in range(len(dic)):
			obj = dic[i]
			Point = FindMatchJson(obj,Data,scale) # Devolve None/False, ou  informacao Objeto 
			if not Point:
				dic[len(dic)] = obj
			else:
				for count, key in enumerate(list(dic[i].keys())):
					dic[i][key].append(float(Point[key]))
	
	return dic


##################################################################################################################################################################

# In Development

def CompareImagesTXT(Images):
	if len(Images) == 1:
		raise Exception('Not enough images to compare!')
	dic = {}
	Image = Images[0]
	OutputFile = f"{DayDirectory}/OutputFiles/{Image.split('.')[0]}"
	Data = SPADEFunctions.ReadOutput(OutputFile)
	for obj in range(len(Data['X'])):
		dic[obj] = {}
		dic['Image'] = [1]
		keys = Data.keys()
		for count, key in enumerate(keys):
			dic[obj][key] = [Data[key][obj]]
	
	
	for count, Image in enumerate(Images[1:-1]):
		OutputFile = f"{DayDirectory}/OutputFiles/{Image.split('.')[0]}"
		Data = SPADEFunctions.ReadOutput(OutputFile)
		time.sleep(500)
		
		for i in range(len(dic)):
			obj = dic[i]
			Valor = FindMatch(obj,Data) # Devolve None/False, ou  informacao Objeto 
			if Valor:
				for count, key in enumerate(dic[obj]):
					dic[obj][key] = dic[obj][key].append(Valor[key])
			else:
				dic[len(dic)] = Valor
