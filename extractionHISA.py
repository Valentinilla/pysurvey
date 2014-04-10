#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class extractionHISA(object):
	
	def __init__(self,obs,spatialConf,spectralConf,analysis):
		
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type
		self.analysis = analysis
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_HISA_'+self.analysis+'Search')
		self.logger.info("Open file and get data...")


		# Get HI emission data
		if self.survey == 'LAB':
			Tb = obs.observation[:,:,:]
		else:
			Tb = obs.observation[0,:,:,:]

		# In case of multiprocessing analysis split the array along the maximum axis
		maxis = 1+argmax(Tb.shape[-2:])
			
		lon = obs.xarray
		lat = obs.yarray
		
		vel = obs.zarray/1000.
		dv = fabs(obs.dz/1000.)		# [velocity] = km s-1
		
		nlon = obs.nx
		nvel = obs.nz

		# free memory
		del obs.observation
		del obs.xarray
		del obs.yarray
		del obs.zarray

		# Array to store results
		cubemap = zeros(Tb.shape,dtype=float32)		
		
		list = []
		if analysis == 'spatial':		
			list = [obs.zmin,obs.zmax,obs.dx,spatialConf]
		if analysis == 'spectral':		
			list = [obs.zmin,obs.zmax,obs.dx,dv,vel,spectralConf]
		
		# Using Multiprocessing if enough cpus are available
		import multiprocessing
		
		ncpu = glob_ncpu
		# Maximum number of cpus
		if ncpu > 16: ncpu = 16
		# Minimum number of cpus
		if Tb.shape[maxis] < ncpu:
			ncpu = Tb.shape[maxis]
		
		self.logger.info("Running on %i cpu(s)"%(ncpu))
		if ncpu > 1:
			import itertools
			
			arrays = array_split(Tb, ncpu, axis=maxis)
			pool = multiprocessing.Pool(processes = ncpu)
			if analysis == 'spatial':
				results = pool.map(spatialSearch, itertools.izip(arrays, itertools.repeat(list)))
			if analysis == 'spectral':
				results = pool.map(spectralSearch, itertools.izip(arrays, itertools.repeat(list)))
			pool.close()
			pool.join()
			cubemap = concatenate(results, axis = maxis)
			del arrays
			del results
		else:
			if analysis == 'spatial':
				cubemap = spatialSearch( (Tb,list) )
			if analysis == 'spectral':
				cubemap = spectralSearch( (Tb,list) )

		results = zeros(cubemap.shape,dtype='b')
		results = where(cubemap<spatialConf('min_hisa'),0,1)
		
		#newheader = pyfits.Header()
		#newheader["ctype1"] = ("GLON-CAR","Coordinate type")
		#newheader["crval1"] = (obs.keyword["crval1"],"Galactic longitude at reference pixel")
		#newheader["crpix1"] = (obs.keyword["crpix1"],"Reference pixel")
		#newheader["cdelt1"] = (obs.keyword["cdelt1"],"Longitude increment")
		#newheader["crota1"] = (obs.keyword["crota1"],"Longitude rotation")
		#newheader["cunit1"] = ("deg","Unit type")
		
		#newheader["ctype2"] = ("GLAT-CAR","Coordinate type")
		#newheader["crval2"] = (obs.keyword["crval2"],"Galactic latitude at reference pixel")
		#newheader["crpix2"] = (obs.keyword["crpix2"],"Reference pixel")
		#newheader["cdelt2"] = (obs.keyword["cdelt2"],"Latitude increment")
		#newheader["crota2"] = (obs.keyword["crota2"],"Latitude rotation")
		#newheader["cunit2"] = ("deg","Unit type")

		#newheader["ctype3"] = (obs.keyword["ctype3"],"Coordinate type")
		#newheader["crval3"] = (obs.keyword["crval3"],"Velocity at reference pixel")
		#newheader["crpix3"] = (obs.keyword["crpix3"],"Reference pixel")
		#newheader["cdelt3"] = (obs.keyword["cdelt3"],"Velocity increment")
		#newheader["crota3"] = (obs.keyword["crota3"],"Velocity rotation")
		#newheader["cunit3"] = ("m/s","Unit type")

		#newheader['bunit'] = (obs.keyword["bunit"],"Map units")
		#newheader['datamin'] = (amin(cubemap),"Min value")
		#newheader['datamax'] = (amax(cubemap),"Max value")
		
		#newheader['minfil'] = unravel_index(argmin(cubemap),cubemap.shape)[0]
		#newheader['mincol'] = unravel_index(argmin(cubemap),cubemap.shape)[1]
		#newheader['minrow'] = unravel_index(argmin(cubemap),cubemap.shape)[2]
		#newheader['maxfil'] = unravel_index(argmax(cubemap),cubemap.shape)[0]
		#newheader['maxcol'] = unravel_index(argmax(cubemap),cubemap.shape)[1]
		#newheader['maxrow'] = unravel_index(argmax(cubemap),cubemap.shape)[2]

		#newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		

		path = getPath(self.logger,'lustre_'+self.survey.lower()+'_hisa')
		# Output file
		#results = pyfits.PrimaryHDU(cubemap,obs.keyword)
		self.logger.info("Writing data to a fits file...")
		#results.writeto(path+self.survey+'_'+self.mosaic+'_'+analysis+'_search.fits', output_verify='fix')
		open(path+self.survey+'_'+self.mosaic+'_'+analysis+'_search.dat', 'wb').write(results)
		self.logger.info("Done")


