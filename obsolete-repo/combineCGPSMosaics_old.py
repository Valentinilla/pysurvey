#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits


class combineCGPSMosaics(object):
	
	def __init__(self,mosaicConf,species,type,flag):
		"""
		Allow to combine mosaics in 2D or 3D regions.
		"""
		self.logger = initLogger(mosaicConf['mosaic'], 'CGPS_CombineMosaics')

		path = '/afs/ifh.de/group/that/work-sf/survey/lab/src_old/info/'
		mosaiclist = path+'mosaic.list'		
		checkForFiles(self.logger,[mosaiclist])
		input = open(mosaiclist,"r")
		list = input.readlines()
		first_mosaic = list[0].split('\n')[0]
		mosaics_number = len(list)
		
		# if species== HI and type==column density				
		path = getPath(self.logger, key='lustre_cgps_hi_column_density')
		ref_mosaic = path+'CGPS_'+first_mosaic+'_HI_unabsorbed_column_density.fits'
		checkForFiles(self.logger,[ref_mosaic])
		
		f = pyfits.open(ref_mosaic)
		hdu = f[0]
		msc_size = hdu.header['NAXIS']	# 4
		msc_lon = hdu.header['NAXIS1']	# 1024
		msc_lat = hdu.header['NAXIS2']	# 1024
		msc_ref_lon = float(hdu.header['CRVAL1']) # deg (GLON-CAR)
		msc_del_lon = float(hdu.header['CDELT1']) # -4.9e-3 deg
		msc_ind_lon = float(hdu.header['CRPIX1']) # 513 px
		msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
		msc_ref_lat = float(hdu.header['CRVAL2']) # deg
		msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
		msc_ind_lat = float(hdu.header['CRPIX2']) # 513 px
		msc_rot_lat = float(hdu.header['CROTA2']) # 0.0

		lon_array = msc_ref_lon + msc_del_lon*(arange(msc_lon)+1.-msc_ind_lon)
		lat_array = msc_ref_lat + msc_del_lat*(arange(msc_lat)+1.-msc_ind_lat)

		glon = msc_ref_lon#114.
		glat = msc_ref_lat#-1.5
		
		# Find index number of central position
		gla = int(round(msc_ind_lon-1.+(glon-msc_ref_lon)/msc_del_lon))
		glb = int(round(msc_ind_lat-1.+(glat-msc_ref_lat)/msc_del_lat))
		
		print glon,gla,lon_array[gla]
		print glat,glb,lat_array[glb]

		exit(0)
		
		db = 224 #214 1 deg
		dl = 450		

		if mosaics_number > 2:
			lat_dim = msc_lat*sqrt(mosaics_number)-db
			lon_dim	= msc_lon*sqrt(mosaics_number)-dl
		else:
			lat_dim = msc_lat*2-db
			lon_dim = msc_lon*2-dl

		skymap = zeros((lat_dim,lon_dim),dtype=float)

		self.logger.info("Number of mosaics: %s"%mosaics_number)
		mosaicList = []
		for m in range(0,mosaics_number):
			msc_name = list[m].split('\n')[0]
			mosaic = path+'CGPS_'+msc_name+'_HI_unabsorbed_column_density.fits'
			checkForFiles(self.logger,[mosaic])
			
			msc_file = pyfits.open(mosaic)
			hdu = msc_file[0]
		
			msc_data = hdu.data
			mosaicList.append(msc_data)
			
			if msc_size > 2:
				self.logger.critical("Mosaic dimension > 2!!")
				os.exit(0)
			
		msc1_array = array(mosaicList[0]) # MV1
		msc2_array = array(mosaicList[1]) # MV2
		msc3_array = array(mosaicList[2]) # MW1
		msc4_array = array(mosaicList[3]) # MW2
		
		b_cpt = msc_lat-db # lat conjunction point
		l_cpt = msc_lon-dl # lon conjunction point
		b_band = 1	   # thickness of the latitudinal band where averaging
		l_band = 1	   # thickness of the longitudinal band where averaging
		
		print b_cpt,l_cpt,lat_dim,lon_dim
		print msc1_array[0:b_cpt,0:lon_dim/2.].shape
		print msc3_array[0:b_cpt,0:lon_dim/2.].shape

		skymap[0:b_cpt,0:lon_dim/2.] = msc1_array[0:b_cpt,0:lon_dim/2.]
		skymap[0:b_cpt,lon_dim/2.:lon_dim] = msc3_array[0:b_cpt,0:lon_dim/2.]

		skymap[b_cpt:lat_dim,0:lon_dim/2.] = msc2_array[0:msc_lat,0:lon_dim/2.]
		skymap[b_cpt:lat_dim,lon_dim/2.:lon_dim] = msc4_array[0:msc_lat,0:lon_dim/2.]

		b_band1 = (msc1_array[b_cpt-b_band:b_cpt+b_band,0:lon_dim/2.]+msc2_array[0:2*b_band,0:lon_dim/2.])/2.
		b_band2 = (msc3_array[b_cpt-b_band:b_cpt+b_band,0:lon_dim/2.]+msc4_array[0:2*b_band,0:lon_dim/2.])/2.
		
		l_band1 = (msc1_array[0:lat_dim/2.,l_cpt-l_band:l_cpt+l_band]+msc2_array[0:lat_dim/2.,0:2*l_band])/2.
		l_band2 = (msc3_array[0:lat_dim/2.,l_cpt-l_band:l_cpt+l_band]+msc4_array[0:lat_dim/2.,0:2*l_band])/2.
		
		#skymap[b_cpt-b_band:b_cpt+b_band,0:lon_dim/2.] = b_band1
		#skymap[b_cpt-b_band:b_cpt+b_band,lon_dim/2.:lon_dim] = b_band2
		
		#skymap[0:lat_dim/2.,l_cpt-l_band:l_cpt+l_band] = l_band1
		#skymap[lat_dim/2.:lat_dim,l_cpt-l_band:l_cpt+l_band] = l_band2
		print skymap.shape

		# Store results
		newheader = pyfits.Header()
		newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
		newheader.update(key="crval1", value=msc_ref_lon, comment="Galactic longitude of reference pixel")
		newheader.update(key="crpix1", value=msc_ind_lon, comment="Reference pixel in lon")
		newheader.update(key="cdelt1", value=msc_del_lon, comment="Longitude increment")
		newheader.update(key="crota1", value=msc_rot_lon, comment="Longitude rotation")
		newheader.update(key="cunit1", value="deg", comment="Unit type")

		newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
		newheader.update(key="crval2", value=msc_ref_lat, comment="Galactic latitude of reference pixel")
		newheader.update(key="crpix2", value=msc_ind_lat, comment="Reference pixel in lat")
		newheader.update(key="cdelt2", value=msc_del_lat, comment="Latitude increment")
		newheader.update(key="crota2", value=msc_rot_lat, comment="Latitude rotation")
		newheader.update(key="cunit2", value="deg", comment="Unit type")

		newheader.update(key="bunit", value="cm-2", comment="Map units")

		newheader.update(key="datamin", value="%e"%amin(skymap))
		newheader.update(key="datamax", value="%e"%amax(skymap))

		newheader.update(key="object", value="CGPS Skymap", comment="GCPS Mosaic")

		results = pyfits.PrimaryHDU(skymap, newheader)
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		results.writeto('skymap_test.fits', output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")

		
