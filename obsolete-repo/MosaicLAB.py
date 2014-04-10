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

class MosaicObsLAB(object):
	def __init__(self, fileName=None):
		if fileName is None:
			fileName = self._obsDialog(fileName)
		self._inputs = 'Mosaic map: '+str(fileName)
		self.fileName = fileName

		#data, header = pyfits.getdata(fileName, 0, header = True)
		f = pyfits.open(fileName)
		hdu = f[0]
		self.msc_hdu = hdu
		self.msc_bitpix = hdu.header['bitpix'] # 16
		#self.msc_bzero = hdu.header['bzero']   # 70
		#self.msc_bscale = hdu.header['bscale'] # 0.0025

		self.msc_area = hdu.header['object'].split('/')[0]
		self.msc_band = hdu.header['freq0']	# 1.420405752000E+09 / rest frequency in Hz
		self.msc_size = hdu.header['NAXIS']	# 3
		self.msc_lon = hdu.header['NAXIS1']	# 721
		self.msc_lat = hdu.header['NAXIS2']	# 361
		self.msc_vel = hdu.header['NAXIS3']	# 891
		self.msc_ref_lon = float(hdu.header['CRVAL1']) # 180 deg (GLON-CAR)
		self.msc_del_lon = float(hdu.header['CDELT1']) # -0.5 deg
		self.msc_ind_lon = float(hdu.header['CRPIX1']) # 1 px
		self.msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
		self.msc_ref_lat = float(hdu.header['CRVAL2']) # -90 deg
		self.msc_del_lat = float(hdu.header['CDELT2']) # 0.5 deg
		self.msc_ind_lat = float(hdu.header['CRPIX2']) # 1 px
		self.msc_rot_lat = float(hdu.header['CROTA2']) # 0.0
		self.msc_ref_vel = float(hdu.header['CRVAL3']/1000.) #convert to km/s
		self.msc_del_vel = float(hdu.header['CDELT3']/1000.) #convert to km/s
		self.msc_ind_vel = float(hdu.header['CRPIX3']) # 1 px
		self.msc_rot_vel = float(hdu.header['CROTA3']) # 0.0

		self.lon_array=self.msc_ref_lon+self.msc_del_lon*(arange(self.msc_lon)+1.-self.msc_ind_lon)
		self.lat_array=self.msc_ref_lat+self.msc_del_lat*(arange(self.msc_lat)+1.-self.msc_ind_lat)
		self.vel_array=self.msc_ref_vel+self.msc_del_vel*(arange(self.msc_vel)+1.-self.msc_ind_vel)

	        self.observation = hdu.data

		filter1 = self.observation < -1e3
		filter2 = isnan(self.observation)
		self.observation[filter1] = 0.
		self.observation[filter2] = 0.

		self.header = hdu.header
		self._inputs = 'Created LAB object'

	#def __getattr__(self, attrname):
		#try;
		#return getattr(self.observation, attrname)
		#except AttributeError:

	def __repr__(self):
		return self._inputs

	def _obsDialog(self, fileName):
		paramDict = MyOrderedDict()
		if fileName is None:
			paramDict['fileName'] = Param('file', '*.fits')
		else:
			paramDict['fileName'] = Param('file', fileName)

		root = SimpleDialog(paramDict, title="Mosaic Elements:")
		root.mainloop()
		output = (paramDict['fileName'].value())
		return output

	def state(self, output=sys.stdout):
		close = False
		try:
			output = open(output, 'w')
			close = False
		except:
			pass
		output.write( "from Mosaic import *\n" )
		output.write( "obs = MosaicObs(fileName=%s)\n" %_quotefn(self.fileName) )
		if close:
			output.close()


