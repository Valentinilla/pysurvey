#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits
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
		if self.survey == 'LAB':
			sur = 'lab'
		if self.survey == 'CGPS':
			sur = 'cgps'
		if self.survey == 'SGPS':
			sur = 'sgps'
		
		flag = ''
		if self.survey == 'LAB':
			# If LAB
			if self.species == 'HI':
				if not load:
					path = getPath(self.logger, key='lab_hi')
					self.mosaic = 'TOT'
					flag = self.species+'_'+'line_image'
				else:
					if self.type == 'brightness temperature':
						path = getPath(self.logger, key='lustre_lab')
						flag = self.species+'_'+'line'
					elif self.type == 'column density':
						path = getPath(self.logger, key='lustre_lab_hi_column_density')
						flag = self.species+'_'+'column_density'
					else:
						typeErrorMsg(self.logger,self.type)
			else:
				self.logger.critical("Only HI species available for "+self.survey+" survey.")
		else:
			if not load:
				if self.species == 'HI':
					path = getPath(self.logger, key=sur+'_hi')
					if self.type == 'brightness temperature':
						flag = self.species+'_'+'line_image'
				elif self.species == 'CO':
					path = getPath(self.logger, key=sur+'_co')
					if self.type == 'brightness temperature':
						flag = self.species+'_'+'line_image'
			else:
				if self.species == 'HI':
					if self.type == 'brightness temperature':
						path = getPath(self.logger, key='lustre_'+sur+'_hi')
						flag = self.species+'_'+'unabsorbed_line'
					elif self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = self.species+'_'+'unabsorbed_column_density'
					else:
						typeErrorMsg(self.logger,self.type)
				
				elif self.species == 'HISA':				
					if self.type == 'brightness temperature':
						path = getPath(self.logger, key='lustre_'+sur+'_hisa')
						flag = self.species+'_'+'line'
					elif self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_hisa_column_density')
						flag = self.species+'_'+'column_density'
					else:
						typeErrorMsg(self.logger,self.type)
				
				elif self.species == 'HI+HISA':
					if self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = self.species+'_'+'column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'CO':				
					if self.type == 'brightness temperature':
						path = getPath(self.logger, key='lustre_'+sur+'_co')
						flag = 'CO_line'
					elif self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_co_column_density')
						flag = 'H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'HI+CO':				
					if self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = 'HI+H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
				
				elif self.species == 'HI+HISA+CO':				
					if self.type == 'column density':
						path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
						flag = 'HI+HISA+H2_column_density'
					else:
						typeErrorMsg(self.logger,self.type,typeEntry=self.species)
	
		filename = path+self.survey+'_'+self.mosaic+'_'+flag+'.fits'
		checkForFiles(self.logger,[filename])
		
		# Open the file and set the variables
		self.filename = filename
		f = pyfits.open(filename)
		hdu = f[0]
		self.msc_hdu = hdu

		if self.survey == 'LAB' or self.survey == 'CGPS' or self.survey == 'SGPS':
			if self.type == 'brightness temperature' and self.survey == 'CGPS':
				self.msc_bitpix = hdu.header['bitpix']
				self.msc_bzero = hdu.header['bzero']
				self.msc_bscale = hdu.header['bscale']
				
			self.msc_size = hdu.header['NAXIS']	# 4
			self.msc_lon = hdu.header['NAXIS1']	# 1024
			self.msc_lat = hdu.header['NAXIS2']	# 1024
			self.msc_ref_lon = float(hdu.header['CRVAL1']) # deg (GLON-CAR)
			self.msc_del_lon = float(hdu.header['CDELT1']) # -4.9e-3 deg
			self.msc_ind_lon = float(hdu.header['CRPIX1']) # 513 px
			self.msc_ref_lat = float(hdu.header['CRVAL2']) # deg
			self.msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
			self.msc_ind_lat = float(hdu.header['CRPIX2']) # 513 px

			if self.survey == 'LAB' or self.survey == 'CGPS':
				self.msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
				self.msc_rot_lat = float(hdu.header['CROTA2']) # 0.0
	
			self.lon_array=self.msc_ref_lon+self.msc_del_lon*(arange(self.msc_lon)+1.-self.msc_ind_lon)
			self.lat_array=self.msc_ref_lat+self.msc_del_lat*(arange(self.msc_lat)+1.-self.msc_ind_lat)
	
			self.observation = hdu.data
	
			if self.type == 'brightness temperature':
				if self.survey == 'LAB':
					self.msc_area = hdu.header['object']
					self.msc_band = hdu.header['freq0'] # 1.4204057E+09 / rest frequency in Hz
				if self.survey == 'CGPS':
					self.msc_area = hdu.header['adc_area']
					self.msc_band = hdu.header['adc_band'] # HI,CO
				if self.survey == 'SGPS':
					self.msc_area = hdu.header['object']
					self.msc_band = hdu.header['restfreq'] # 1.4204057E+09 / rest frequency in Hz

				self.msc_vel = hdu.header['NAXIS3'] # 272
				if self.survey == 'CGPS' or self.survey == 'SGPS':
					self.msc_pol = hdu.header['NAXIS4']	# 1, polarization

				self.msc_ref_vel = float(hdu.header['CRVAL3']/1000.) #convert to km/s
				self.msc_del_vel = float(hdu.header['CDELT3']/1000.) #convert to km/s
				self.msc_ind_vel = float(hdu.header['CRPIX3']) # 145 px

				if self.survey == 'LAB' or self.survey == 'CGPS':
					self.msc_rot_vel = float(hdu.header['CROTA3']) # 0.0
				if self.survey == 'CGPS':
					if self.msc_band == 'HI':	
						self.observation[:,:18,:,:] = 0.
						self.observation[:,271,:,:] = 0.
						filter1 = self.observation < -1e4
						filter2 = isnan(self.observation)
						self.observation[filter1] = 0.
						self.observation[filter2] = 0.
					if self.msc_band == 'CO':	
						#self.observation[:,:23,:,:] = 0.
						#self.observation[:,256:,:,:] = 0.
						filter1 = self.observation < -1e3
						filter2 = isnan(self.observation)
						self.observation[filter1] = 0.
						self.observation[filter2] = 0.
		
				if self.survey == 'LAB':
					filter1 = self.observation < -1e3
					filter2 = isnan(self.observation)
					self.observation[filter1] = 0.
					self.observation[filter2] = 0.
				
				self.vel_array=self.msc_ref_vel+self.msc_del_vel*(arange(self.msc_vel)+1.-self.msc_ind_vel)
			
			if self.type == 'integrated brightness temperature':	
				filter1 = self.observation < -1e4
				filter2 = isnan(self.observation)
				self.observation[filter1] = 0.
				self.observation[filter2] = 0.
					
		self.header = hdu.header
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


