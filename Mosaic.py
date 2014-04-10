#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import random


try:
	from SimpleDialog import SimpleDialog, map, Param
except ImportError:
	pass

class Mosaic(object):
	
	def __init__(self,surveyConf,mosaicConf,type,species='',load=False):
		"""
		Read the fits file: HI, HISA, CO
		"""
		self.survey = surveyConf['survey']
		self.mosaic = mosaicConf['mosaic']
		self.species = species
		if self.species == '':
			self.species = surveyConf['species']
		self.type = type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_Mosaic')

		# Select the file according to Survey and Mosaic
		path = ''
		flag = ''
		sur = (self.survey).lower()

		if self.survey == 'Galprop':
			# If Galprop
			if self.species == 'HI':
				if not load:
					path = getPath(self.logger, key=sur+'_hi')
					self.mosaic = 'TOT'
					flag = self.species+'_column_density_rbands_image'
				else:
					if self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = self.species+'_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
			elif self.species == 'WCO':
					if not load:
						if self.type == glob_ITb:
							path = getPath(self.logger, key=sur+'_co')
							self.mosaic = 'TOT'
							flag = self.species+'_rbands_image'
						else:
							typeErrorMsg(self.logger,self.type,self.species)
					else:
						if self.type == glob_ITb:
							path = getPath(self.logger, key='lustre_'+sur+'_co')
							flag = self.species+'_line'
						elif self.type == glob_N:
							path = getPath(self.logger, key='lustre_'+sur+'_co_column_density')
							flag = self.species+'_column_density'
						else:
							typeErrorMsg(self.logger,self.type)
			else:
				self.logger.critical("Only HI and WCO species available for "+self.survey+" survey.")
		
		elif self.survey == 'LAB':
			# If LAB
			if self.species == 'HI':
				if not load:
					path = getPath(self.logger, key='lab_hi')
					self.mosaic = 'TOT'
					flag = self.species+'_line_image'
				else:
					if self.type == glob_Tb:
						path = getPath(self.logger, key='lustre_lab')
						flag = self.species+'_line'
					elif self.type == glob_N:
						path = getPath(self.logger, key='lustre_lab_hi_column_density')
						flag = self.species+'_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
			else:
				self.logger.critical("Only HI species available for "+self.survey+" survey.")
		
		elif self.survey == 'Dame':
			# If Dame
			if self.species == 'WCO':
				if not load:
					path = getPath(self.logger, key='dame_co')
					if self.type == glob_ITb:
						self.mosaic = 'TOT'
						flag = self.species+'_line_image'
				else:
					if self.type == glob_ITb:
						path = getPath(self.logger, key='lustre_dame')
						flag = self.species+'_line'
					elif self.type == glob_N:
						path = getPath(self.logger, key='lustre_dame_co_column_density')
						flag = self.species+'_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
			else:
				self.logger.critical("Only WCO species available for "+self.survey+" survey.")
		
		else:
			if not load:
				if self.species == 'HI':
					path = getPath(self.logger, key=sur+'_hi')
					if self.type == glob_Tb:
						flag = self.species+'_line_image'
				elif self.species == 'CO':
					path = getPath(self.logger, key=sur+'_co')
					if self.type == glob_Tb:
						flag = self.species+'_line_image'
			else:
				if self.species == 'HI':
					if self.type == glob_Tb:
						path = getPath(self.logger, key='lustre_'+sur+'_hi')
						flag = self.species+'_unabsorbed_line'
					elif self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = self.species+'_unabsorbed_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
				
				elif self.species == 'HISA':				
					if self.type == glob_Tb:
						path = getPath(self.logger, key='lustre_'+sur+'_hisa')
						flag = self.species+'_line'
					elif self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hisa_column_density')
						flag = self.species+'_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
				
				elif self.species == 'HI+HISA':
					if self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = self.species+'_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'CO':				
					if self.type == glob_Tb:
						path = getPath(self.logger, key='lustre_'+sur+'_co')
						flag = 'CO_line'
					elif self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_co_column_density')
						flag = 'H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'HI+CO':				
					if self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = 'HI+H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'HI+HISA+CO':				
					if self.type == glob_N:
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = 'HI+HISA+H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
	
		filename = path+self.survey+'_'+self.mosaic+'_'+flag+'.fits'
		checkForFiles(self.logger,[filename])
		
		# Open the file and set the variables
		self.filename = filename
		f = pyfits.open(filename)
		self.keyword = f[0].header
		
		# Bscale and bzero should be read as first keywords	
		if ('BSCALE' and 'BZERO') in self.keyword:	
			self.bscale = self.keyword['bscale']
			self.bzero = self.keyword['bzero']

		if not('CROTA1' or 'CROTA2') in self.keyword:
			self.keyword['CROTA1'] = 0.0
			self.keyword['CROTA2'] = 0.0
		
		# Build arrays
		try:
			self.x,self.y = self.keyword['CRVAL1'],self.keyword['CRVAL2']
			self.dx,self.dy = self.keyword['CDELT1'],self.keyword['CDELT2']
			self.x_px,self.y_px = self.keyword['CRPIX1'],self.keyword['CRPIX2']
			self.nx,self.ny = self.keyword['NAXIS1'],self.keyword['NAXIS2']
			self.xarray = self.x + self.dx * (arange(self.nx) + 1. - self.x_px) #self.lon_array
			self.yarray = self.y + self.dy * (arange(self.ny) + 1. - self.y_px) #self.lat_array
		except:
			self.logger.critical("Some keyword missing. Coordinate arrays cannot be built.")
		
		if 'ADC_AREA' in self.keyword:
			self.object = self.keyword['ADC_AREA']
			del self.keyword['ADC_AREA']
			self.keyword['OBJECT'] = self.object	
		if 'FREQ0' in self.keyword:
			self.band = self.keyword['FREQ0']
			del self.keyword['FREQ0']
			self.keyword['BAND'] = self.band
		elif 'RESTFREQ' in self.keyword:
			self.band = self.keyword['RESTFREQ']
			del self.keyword['RESTFREQ']
			self.keyword['BAND'] = self.band
		elif 'ADC_BAND' in self.keyword:
			self.band = self.keyword['ADC_BAND']
			del self.keyword['ADC_BAND']
			self.keyword['BAND'] = self.band
		
		if self.keyword['NAXIS'] > 2:
			if not 'CROTA3' in self.keyword:
				self.keyword['CROTA3'] = 0.0
			# Build array
			try:
				self.z = self.keyword['CRVAL3']
				self.dz = self.keyword['CDELT3']
				self.z_px = self.keyword['CRPIX3']
				self.nz = self.keyword['NAXIS3']
				self.zarray = self.z + self.dz * (arange(self.nz) + 1. - self.z_px)
			except:
				self.logger.critical("Some keyword missing. 3rd axis array cannot be built.")
		
		self.observation = f[0].data
		# free memory
		del f[0].data
		
		if self.type == glob_Tb or self.type == glob_ITb:	
			filter1 = self.observation < -1e4
			filter2 = isnan(self.observation)
			self.observation[filter1] = 0.
			self.observation[filter2] = 0.
			if self.keyword['NAXIS'] > 3:
				if self.survey == 'CGPS':
					if 'BAND' in self.keyword:
						if self.keyword['BAND'] == 'HI':	
							self.observation[:,:18,:,:] = 0.
							self.observation[:,271,:,:] = 0.
						#if self.keyword['BAND'] == 'CO':	
							#self.observation[:,:23,:,:] = 0.
							#self.observation[:,256:,:,:] = 0.
		
		self.mosaic = mosaicConf['mosaic']
		if not load:
			self._inputs = 'Created '+self.survey+' Mosaic object '+self.species+' '+self.type
		else:
			self._inputs = 'Loaded '+self.survey+' Mosaic object '+self.species+' '+self.type
	
	#def __getattr__(self, attrname):
		#try;
		#return getattr(self.observation, attrname)
		#except AttributeError:

	def __repr__(self):
		return self._inputs

	def _obsDialog(self, filename):
		paramDict = MyOrderedDict()
		if filename is None:
			paramDict['filename'] = Param('file', '*.fits')
		else:
			paramDict['filename'] = Param('file', filename)

		root = SimpleDialog(paramDict, title="Mosaic Elements:")
		root.mainloop()
		output = (paramDict['filename'].value())
		return output

	def state(self, output=sys.stdout):
		close = False
		try:
			output = open(output, 'w')
			close = False
		except:
			pass
		output.write( "from Mosaic import *\n" )
		output.write( "obs = MosaicObs(filename=%s)\n" %_quotefn(self.filename) )
		if close:
			output.close()


