from skyfield.api import load, wgs84

stations_url = 'http://celestrak.com/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'
satellites = load.tle_file(stations_url)
print(satellites)
print('Loaded', len(satellites), 'satellites')
