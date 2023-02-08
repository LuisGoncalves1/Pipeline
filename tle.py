import json
import requests
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84

'''
ReadTLES(Links) function uploads tle information from a url and saves it in a .json file (a list of dicts).
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
SatellitePosition() loads the satellites information from the TLEs saved in the .json file from the ReadTLES(Links) function, 
and calculates their position in a given time and location.
'''

def SatellitePosition():
	with open("GEO_tles.json", 'r') as stream:
		lista_ficheiro=json.load(stream)
	
	ts = load.timescale() #The Timescale object returned by load.timescale() manages the conversions between different 
						  #time scales and is also how the programmer builds Time objects for specific dates. 
	for i in lista_ficheiro:

		line1=i['tle_1']
		line2=i['tle_2']
		nome=i['satellite_name']
	
		satellite = EarthSatellite(line1, line2, nome, ts) # An Earth satellite loaded from a TLE file and propagated with SGP4. 

		bluffton = wgs84.latlon(+40.182524, -7.873811) #coordinates of latitude and longitude of pampilhosa

		t0 = ts.utc(2023, 2, 8, 23, 47, 53) #today's date, time when the first image was takes
		t1 = ts.utc(2023, 2, 8, 23, 48, 58) #today's date, time when the last image was takes
		
		t, events = satellite.find_events(bluffton, t0, t1, altitude_degrees=30.0) 
		
		# find_events(topos, t0, t1, altitude_degrees=0.0) - Searches between t0 and t1, which should each be a Skyfield Time 
		# object, for passes of this satellite above the location topos that reach at least altitude_degrees above the horizon
		for ti, event in zip(t, events):
			
			name = ('rise above 30°', 'culminate', 'set below 30°')[event] #culminate-the moment when it is nearest to the Zenith
			print(f"Satellite {nome}: {ti.utc_strftime('%Y %b %d %H:%M:%S'),name}")
			
			difference = satellite - bluffton 
			topocentric = difference.at(ti)
			
			alt, az, distance = topocentric.altaz()
			print('Altitude:', alt)
			print('Azimuth:', az)

			ra, dec, distance = topocentric.radec()
			print(f"Right Ascention: {ra}")
			print(f"Declination: {dec}")

b=SatellitePosition()
