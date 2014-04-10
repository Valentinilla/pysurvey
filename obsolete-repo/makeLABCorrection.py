#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits


class makeLABCorrection(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf):

		self.logger = initLogger(mosaicConf['mosaic'], 'LAB_ColumnDensity')

		path = getPath(self.logger, key="lustre_lab_hi_column_density")
		file = path+'LAB_'+mosaicConf['mosaic']+'_HI_column_density.fits'
		checkForFiles(self.logger,[file],existence=True)

		self.logger.info("Open file and get data...")
		
		# Array to store results
		NHI = zeros((mosaic.msc_lat,mosaic.msc_lon),dtype=float)
		
		Ts = amax(mosaic.observation)	  # [Excitation (or Spin) Temperature] = K (150)
		Tbg = float(utilsConf['tcmb'])	  # [Cosmic Microwave Background (CMB)] = K
		dv = abs(mosaic.msc_del_vel)	  # [velocity] = km s-1
		C = float(utilsConf['c'])	  # [costant] = cm-2
		
		# Corrections for latitude (almost negligible)
		#pxsize = mosaic.msc_del_lat*(pi/180.)	 # [pixel size] = rad
		#Apx = power(pxsize, 2)	 		 # [Area pixel] = sr
		cosdec = cos( (pi/180.) * mosaic.lat_array) # [lat. correction] = rad

		# Get LAB emission data		
		Tb = mosaic.observation[:,:,:]

		# Setting the negative/0-values to the background temperature 
		#Tb = where( (Tb<0.) | (Tb==0.),Tbg,Tb)

		# Without the continuum component
		index = where(Tb>=Ts)
		Tb[index] = 0.999*Ts # Tb must be < Ts, otherwise problems with log

		self.logger.info("Calculating NHI...")
		# Optical depth correction
		# With the continuum component
		#cTb = log( (Ts-Tc)/(Ts-Tc-Tb) ) * Ts
		# Without the continuum component
		cTb = log( (Ts)/(Ts-Tb) ) * Ts
		
		# Integrated brightness temperature over velocity (axis = 0)
		ITb = sum(cTb,axis=0)
		
		# Column density
		NHI = C * ITb * dv	# [NHI] = cm-2
		# Corrected column density
		NHI = NHI*cosdec
		
		# python convention
		# a[y,x] i.e., indexing from right to left
		# y/j = latitude, x/i = longitude
		
		newheader = pyfits.Header()
		newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
		newheader.update(key="crval1", value=mosaic.msc_ref_lon, comment="Galactic longitude of reference pixel")
		newheader.update(key="crpix1", value=mosaic.msc_ind_lon, comment="Reference pixel in lon")
		newheader.update(key="cdelt1", value=mosaic.msc_del_lon, comment="Longitude increment")
		newheader.update(key="crota1", value=mosaic.msc_rot_lon, comment="Longitude rotation")
		newheader.update(key="cunit1", value="deg", comment="Unit type")

		newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
		newheader.update(key="crval2", value=mosaic.msc_ref_lat, comment="Galactic latitude of reference pixel")
		newheader.update(key="crpix2", value=mosaic.msc_ind_lat, comment="Reference pixel in lat")
		newheader.update(key="cdelt2", value=mosaic.msc_del_lat, comment="Latitude increment")
		newheader.update(key="crota2", value=mosaic.msc_rot_lat, comment="Latitude rotation")
		newheader.update(key="cunit2", value="deg", comment="Unit type")

		newheader.update(key="proj", value="FLAT", comment="Type of projection")
		newheader.update(key="system", value="GALACTIC", comment="Coordinate system")
		newheader.update(key="bunit", value="cm-2", comment="Map units")
		newheader.update(key="equinox", value=2000., comment="Equinox of ref. coord.")
		newheader.update(key="datamin", value="%e"%amin(NHI))
		newheader.update(key="datamax", value="%e"%amax(NHI))
		newheader.update(key="object", value="LAB Mosaic %s"%mosaicConf['mosaic'], comment="Equivalent of CGPS Mosaic")

		results = pyfits.PrimaryHDU(NHI, newheader)
		#results.scale('int16', '', bscale=mosaic.msc_bscale, bzero=mosaic.msc_bzero)

		# Output file
		self.logger.info("Write scaled data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")

