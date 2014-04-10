#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits


class makeLABMosaic(object):
	
	def __init__(self,obs,mosaicConf):
		"""
		Allow to generate CGPS-like mosaics for LAB Survey.
		"""	
		
		self.logger = initLogger(mosaicConf['mosaic'], 'LAB_GenerateMosaic')
		self.logger.info("Open file and get data...")

		lon = float(mosaicConf['lon'])
		lat = float(mosaicConf['lat'])

		v1 = int(mosaicConf['vel1'])
		v2 = int(mosaicConf['vel2'])
	
		side = float(mosaicConf['side'])/2.

		l1,l2,b1,b2 = getMosaicCoordinate(obs,lon,lat,side)
		
		crpix1 = rint((b2-b1)/2.)
		crpix2 = rint((l2-l1)/2.)
		#print l1,l2, b1,b2

		self.logger.info("Mosaic properties:")
		self.logger.info("- (l0,b0) = (%s,%s) deg"%(lon,lat))
		self.logger.info("- (x0,y0) = (%s,%s) px"%(l1,b1))
		self.logger.info("- (v1,v2) = (%s,%s) px"%(v1,v2))
		self.logger.info("- (h,w) = (%s,%s) deg"%(side,side))
		self.logger.info("- (h,w) = (%s,%s) px"%(crpix1,crpix2))
		
		newmosaic = obs.observation[(v1-1):v2,b1:b2,l1:l2]
		#newmosaic = obs.observation[:,b1:b2,l1:l2]

		newheader = pyfits.Header()
		newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
		newheader.update(key="crval1", value=lon, comment="Galactic longitude of reference pixel")
		newheader.update(key="crpix1", value=crpix1, comment="Reference pixel in lon")
		newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
		newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
		newheader.update(key="cunit1", value="deg", comment="Unit type")

		newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
		newheader.update(key="crval2", value=lat, comment="Galactic latitude of reference pixel")
		newheader.update(key="crpix2", value=crpix2, comment="Reference pixel in lat")
		newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
		newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
		newheader.update(key="cunit2", value="deg", comment="Unit type")

		newheader.update(key="ctype3", value="VELO-LSR", comment="Coordinate type")
		#newheader.update(key="crval3", value=-164.892, comment="Velocity of reference pixel")
		newheader.update(key="crval3", value=obs.vel_array[min(v1,v2)-1]*1e3, comment="Velocity of reference pixel")
		newheader.update(key="crpix3", value=obs.msc_ind_vel, comment="Reference pixel in vel")
		newheader.update(key="cdelt3", value=obs.msc_del_vel*1e3, comment="Velocity increment")
		newheader.update(key="crota3", value=obs.msc_rot_vel, comment="Velocity rotation")
		newheader.update(key="cunit3", value="m/s", comment="Unit type")

		newheader.update(key="proj", value="FLAT", comment="Type of projection")
		newheader.update(key="system", value="GALACTIC", comment="Coordinate system")
		newheader.update(key="bunit", value="K", comment="Map units")
		newheader.update(key="equinox", value=2000., comment="Equinox of ref. coord.")
	
		newheader.update(key="datamin", value=amin(newmosaic))
		newheader.update(key="datamax", value=amax(newmosaic))

		newheader.update(key="object", value="LAB Mosaic %s"%mosaicConf['mosaic'], comment="Equivalent of CGPS Mosaic")
		newheader.update(key="freq0", value=obs.msc_band, comment="Rest frequency in Hz")

		results = pyfits.PrimaryHDU(newmosaic,newheader)
		#results.scale('int16', '', bscale=obs.msc_bscale, bzero=obs.msc_bzero)
		
		# Output file
		#self.logger.info("Write scaled data (int16) to a fits file...")
		self.logger.info("Write data to a fits file in...")
		path = getPath(self.logger, key="lustre_lab")
		results.writeto(path+'LAB_'+mosaicConf['mosaic']+'_HI_line.fits', output_verify='fix')
		self.logger.info("%s"%path)

		self.logger.info("Done")




