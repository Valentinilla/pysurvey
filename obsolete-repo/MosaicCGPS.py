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

class MosaicObsCGPS(object):
	
	def __init__(self, mosaicFile,type):
		"""
		Read the fits file: HI, HISA, CO
		"""

		self._inputs = 'Mosaic map: '+str(mosaicFile)
		self.mosaicFile = mosaicFile
		self.mosaicType = type
		
		f = pyfits.open(mosaicFile)
		hdu = f[0]
		self.msc_hdu = hdu
		
		if self.mosaicType == 'brightness temperature':
			self.msc_bitpix = hdu.header['bitpix']
			self.msc_bzero = hdu.header['bzero']
			self.msc_bscale = hdu.header['bscale']

		self.msc_size = hdu.header['NAXIS']	# 4
		self.msc_lon = hdu.header['NAXIS1']	# 1024
		self.msc_lat = hdu.header['NAXIS2']	# 1024
		self.msc_ref_lon = float(hdu.header['CRVAL1']) # deg (GLON-CAR)
		self.msc_del_lon = float(hdu.header['CDELT1']) # -4.9e-3 deg
		self.msc_ind_lon = float(hdu.header['CRPIX1']) # 513 px
		#self.msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
		self.msc_ref_lat = float(hdu.header['CRVAL2']) # deg
		self.msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
		self.msc_ind_lat = float(hdu.header['CRPIX2']) # 513 px
		#self.msc_rot_lat = float(hdu.header['CROTA2']) # 0.0

		self.lon_array=self.msc_ref_lon+self.msc_del_lon*(arange(self.msc_lon)+1.-self.msc_ind_lon)
		self.lat_array=self.msc_ref_lat+self.msc_del_lat*(arange(self.msc_lat)+1.-self.msc_ind_lat)

		self.observation = hdu.data

		if self.mosaicType == 'brightness temperature':
			#self.msc_area = hdu.header['adc_area']
			self.msc_band = ''#hdu.header['adc_band'] # HI,CO
			self.msc_vel = hdu.header['NAXIS3']	# 272
			self.msc_pol = hdu.header['NAXIS4']	# 1, polarization
			self.msc_ref_vel = float(hdu.header['CRVAL3']/1000.) #convert to km/s
			self.msc_del_vel = float(hdu.header['CDELT3']/1000.) #convert to km/s
			self.msc_ind_vel = float(hdu.header['CRPIX3']) # 145 px
			#self.msc_rot_vel = float(hdu.header['CROTA3']) # 0.0
	
			self.vel_array=self.msc_ref_vel+self.msc_del_vel*(arange(self.msc_vel)+1.-self.msc_ind_vel)

			if self.msc_band == 'HI' or self.msc_band == 'CO':	
				self.observation[:,:18,:,:] = 0.
				self.observation[:,271,:,:] = 0.
				filter1 = self.observation < -1e4
				filter2 = isnan(self.observation)
				self.observation[filter1] = 0.
				self.observation[filter2] = 0.

		if self.mosaicType == 'integrated brightness temperature':	
			filter1 = self.observation < -1e4
			filter2 = isnan(self.observation)
			self.observation[filter1] = 0.
			self.observation[filter2] = 0.
			
				
		self.header = hdu.header

		if self.mosaicType == 'brightness temperature':
			self._inputs = 'Created CGPS Tb Mosaic object'
		if self.mosaicType == 'column density':
			self._inputs = 'Created CGPS NHI Mosaic object'
		if self.mosaicType == 'integrated brightness temperature':
			self._inputs = 'Created CGPS ITb Mosaic object'

	#def __getattr__(self, attrname):
		#try;
		#return getattr(self.observation, attrname)
		#except AttributeError:

	def __repr__(self):
		return self._inputs

	def _obsDialog(self, mosaicFile):
		paramDict = MyOrderedDict()
		if mosaicFile is None:
			paramDict['mosaicFile'] = Param('file', '*.fits')
		else:
			paramDict['mosaicFile'] = Param('file', mosaicFile)

		root = SimpleDialog(paramDict, title="Mosaic Elements:")
		root.mainloop()
		output = (paramDict['mosaicFile'].value())
		return output

	def state(self, output=sys.stdout):
		close = False
		try:
			output = open(output, 'w')
			close = False
		except:
			pass
		output.write( "from Mosaic import *\n" )
		output.write( "obs = MosaicObs(mosaicFile=%s)\n" %_quotefn(self.mosaicFile) )
		if close:
			output.close()


