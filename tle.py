import json
import requests
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84

'''
ReadTLES(Links) function uploads tle information from a given url and saves it in a .json file (a list of dicts).
Currently, I am using celestrak to extract the tles, and foccusing on GeoProtectedZone and GeoProtectedZonePlus satellites.
'''

def ReadTLES(Links):
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
		
					
#a=ReadTLES(['https://celestrak.org/NORAD/elements/gp.php?SPECIAL=gpz&FORMAT=tle',
#'https://celestrak.org/NORAD/elements/gp.php?SPECIAL=gpz-plus&FORMAT=tle'])

'''
DetermineSatellitePosition() is still in construction
'''

def DetermineSatellitePosition():
	with open("GEO_tles.json", 'r') as stream:
		lista_ficheiro=json.load(stream)
		
	ts = load.timescale()

	for i in lista_ficheiro:
	
		line1=i['tle_1']
		line2=i['tle_2']
		name=i['satellite_name']
	
		# An Earth satellite loaded from a TLE file and propagated with SGP4. 
		satellite = EarthSatellite(line1, line2, name, ts)
		#print(satellite)
	
		#The “epoch” date of a satellite element set is the all-important date and time for which the set of elements is most accurate,
		#and before or after which they go rapidly out of date.
		print(f"Epoch data:{satellite.epoch.utc_jpl()}")
	
