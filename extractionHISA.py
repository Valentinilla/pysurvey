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
		
		list = []
		if analysis == 'spatial':		
			list = [obs.zmin,obs.zmax,obs.dx,obs.nx,obs.ny,obs.nz,spatialConf]
		if analysis == 'spectral':		
			list = [obs.zmin,obs.zmax,obs.dx,obs.dz,obs.nx,obs.ny,obs.nz,obs.zarray,spectralConf]

		# Using Multiprocessing if enough cpus are available
		import multiprocessing

		ncpu = 1# glob_ncpu
		if ncpu > 16: ncpu = 16
		dim = min(obs.observation.shape[1:])
		if dim < ncpu:
			ncpu = dim
		
		self.logger.info("Running on %i cpu(s)"%(ncpu))
		if ncpu > 1:
			import itertools
			if analysis == 'spatial':	
				ax = argmin(obs.observation.shape[2:]) + 1
				samples_list = array_split(obs.observation[0,:,:,:], ncpu, axis = ax)
				pool = multiprocessing.Pool(processes = ncpu)
				results = pool.map(spatialSearch, itertools.izip(samples_list, itertools.repeat(list)))
				pool.close()
				pool.join()
				obs.observation[0,:,:,:] = concatenate(results, axis = ax)
				del samples_list
				del list_spat
				del results
			if analysis == 'spectral':
				ax = argmin(obs.observation.shape[2:]) + 1
				samples_list = array_split(obs.observation[0,:,:,:], ncpu, axis = ax)
				pool = multiprocessing.Pool(processes = ncpu)
				results = pool.map(spectralSearch, itertools.izip(samples_list, itertools.repeat(list)))
				pool.close()
				pool.join()
				obs.observation[0,:,:,:] = concatenate(results, axis = ax)
				del samples_list
				del list_spat
				del results
		else:
			if analysis == 'spatial':
				obs.observation[0,:,:,:] = spatialSearch( (obs.observation[0,:,:,:],list) )
			if analysis == 'spectral':
				obs.observation[0,:,:,:] = spectralSearch( (obs.observation[0,:,:,:],list) )

		path = getPath(self.logger,'lustre_'+self.survey.lower()+'_hisa')
		# Output file
		results = pyfits.PrimaryHDU(obs.observation,obs.keyword)
		self.logger.info("Writing data to a fits file...")
		results.writeto(path+self.survey+'_'+self.mosaic+'_'+analysis+'_search.fits', output_verify='fix')
		self.logger.info("Done")


