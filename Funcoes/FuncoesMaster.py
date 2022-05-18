import os
from astropy.nddata import CCDData
import ccdproc



def CreateMasters(DayDirectory):
	Darks = OrganizeDarks(DayDirectory) # Sees available Darks
	Flats = os.listdir(DayDirectory + '/Flats') # Sees available Flats
	Masters = SeeMasters(DayDirectory)	

	for Filter in Flats:
		if Filter not in Masters['Flat']:
			FFiles = os.listdir(DayDirectory + '/Flats/' + Filter)
			FlatFiles =[]
			for i in FFiles:
				FlatFiles.append(DayDirectory + '/Flats/' + Filter + '/' + i)
			CreateMaster(FlatFiles,DayDirectory,Filter,'Flat')
			
	for Time in Darks:
		if Time not in Masters['Dark']:
			DFiles = Darks[Time]
			DarkFiles = []
			for i in DFiles:
				DarkFiles.append(DayDirectory + '/Darks/' + i)
			CreateMaster(DarkFiles,DayDirectory,Time,'Dark')
			
def CreateMaster(Files,DayDirectory,Filter,Type):
	FilesRead = [CCDData.read(X, unit='adu', Format = 'fits') for X in Files]
	Master = ccdproc.combine(FilesRead,method='median')
	Master.header['FILTER'] = Filter
	
	MasterName = "Master%s%s.fits" % (Type, Filter)
	
	Master.write(DayDirectory + '/Masters/' + MasterName)
	return Master
		
def OrganizeDarks(DayDirectory):
	''' Returns DarkFiles available in a dictionary of the type:
	Dic = {100:[Dark100_001.fits,Dark100_002.fits],200:[Dark200_001.fits,Dark200_002.fits]}
	
	
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
	''' Returns a dictionary of the type:
	Dic = {Flat: ['A','B','C']; Dark:['50','70','100','200']}
	
	
	
	
	
	
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
		
def BestDark(Directory,time,Filter):
	files = os.listdir(Directory)
	files = [B for B in files if B[6:10] == 'Dark' and int(B[10:13]) <= time]
	# ~ files = [B for B in files if B[6:10] == 'Dark' and B[10] == Filter and int(B[11:14]) <= time]
	times = [B[11:14] for B in files]
	if times == []:
		return []
	else:
		return files[times.index(max(times))]
