#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits


class combineSMosaics(object):
	
	def __init__(self,mosaicConf,species,type,flag):
		"""
		Allow to combine mosaics in 2D or 3D regions.
		"""
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_CombineMosaics')

		path = '/afs/ifh.de/group/that/work-sf/survey/lab/src_old/info/'
		mosaiclist = path+'mosaic.list'		
		checkForFiles(self.logger,[mosaiclist])
		input = open(mosaiclist,"r")
		list = input.readlines()
		first_mosaic = list[0].split('\n')[0]
		mosaics_number = len(list)
		sqrt_msc_num = sqrt(mosaics_number)

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
		#msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
		msc_ref_lat = float(hdu.header['CRVAL2']) # deg
		msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
		msc_ind_lat = float(hdu.header['CRPIX2']) # 513 px
		#msc_rot_lat = float(hdu.header['CROTA2']) # 0.0

		#lon_array = msc_ref_lon + msc_del_lon*(arange(msc_lon)+1.-msc_ind_lon)
		#lat_array = msc_ref_lat + msc_del_lat*(arange(msc_lat)+1.-msc_ind_lat)

		msc_lside_px = 2*msc_ind_lon
		msc_lside_deg = msc_lside_px*msc_del_lon

		msc_bside_px = 2*msc_ind_lat
		msc_bside_deg = msc_bside_px*msc_del_lat

		band_lon_px = 450 #450
		band_lon_deg = band_lon_px*msc_del_lon

		band_lat_px = 224 #224
		band_lat_deg = band_lat_px*msc_del_lat		

		l0 = msc_ref_lon # deg
		b0 = msc_ref_lat # deg

		skmp_ref_lon = l0 + (sqrt_msc_num-1) * msc_lside_deg/2 - band_lon_deg/2 # deg
		skmp_ref_lat = b0 + (sqrt_msc_num-1) * msc_bside_deg/2 - band_lat_deg/2 # deg
		
		print "Skymap:"
		print "- l0, b0: %s, %s (deg)"%(skmp_ref_lon,skmp_ref_lat)

		# Find index number of central position
		#gla = int(round(msc_ind_lon-1.))#+(glon-msc_ref_lon)/msc_del_lon))
		#glb = int(round(msc_ind_lat-1.))#+(glat-msc_ref_lat)/msc_del_lat))
				
		if mosaics_number > 2:
			lat_dim = msc_lat*sqrt_msc_num - band_lat_px
			lon_dim	= msc_lon*sqrt_msc_num - band_lon_px
		else:
			lat_dim = msc_lat*2 - band_lat_px
			lon_dim = msc_lon*2 - band_lon_px

		skymap = zeros((lat_dim,lon_dim),dtype=float)
		print "- shape: ",skymap.shape

		self.logger.info("Number of mosaics: %s"%mosaics_number)
		mosaicList = []
		for m in range(0,mosaics_number):
			msc_name = list[m].split('\n')[0]
			mosaic = path+'CGPS_'+msc_name+'_HI_unabsorbed_column_density.fits'
			checkForFiles(self.logger,[mosaic])
			
			msc_file = pyfits.open(mosaic)
			hdu = msc_file[0]

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

			msc_data = hdu.data
			mosaicList.append([msc_data,msc_ref_lon,msc_ref_lat])
			
			if msc_size > 2:
				self.logger.critical("Mosaic dimension > 2!!")
				os.exit(0)

			if m <= sqrt_msc_num-1:
				skmp_idxl = int(round(lon_dim/sqrt_msc_num-1+(msc_ref_lon-skmp_ref_lon)/msc_del_lon))
			else:
				skmp_idxl = int(round(lon_dim/sqrt_msc_num-1+(msc_ref_lon-skmp_ref_lon)/msc_del_lon))-band_lon_px/2
			
			skmp_idxb = int(round(lat_dim/sqrt_msc_num-1+(msc_ref_lat-skmp_ref_lat)/msc_del_lat))
		
			skmp_deg_lon = msc_ref_lon + (sqrt_msc_num-1) * msc_lside_deg/2 - band_lon_deg/2 # deg
			skmp_deg_lat = msc_ref_lat + (sqrt_msc_num-1) * msc_bside_deg/2 - band_lat_deg/2 # deg
	


			l1 = skmp_idxl - msc_ind_lon 
			l2 = skmp_idxl + msc_ind_lon
			b1 = skmp_idxb - msc_ind_lat
			b2 = skmp_idxb + msc_ind_lat
			print "Mosaic n %s"%m
			print "- lon0 = %s, lat0 = %s"%(msc_ind_lon,msc_ind_lat)
			print "- lon1 = %s, lon2 = %s"%(msc_ind_lon-msc_ind_lon,2*msc_ind_lon)
			print "- lat1 = %s, lat2 = %s"%(msc_ind_lat-msc_ind_lat,2*msc_ind_lat)
			print "- l0, b0 = %s, %s (deg)"%(msc_ref_lon,msc_ref_lat)
			#print "- lon_del, lat_del = %s, %s"%(msc_del_lon,msc_del_lat)
			print "Skymap"
			print "- lon0 = %s, lat0 = %s"%(skmp_idxl,skmp_idxb)
			print "- lon1 = %s, lon2 = %s"%(l1,l2)
			print "- lat1 = %s, lat2 = %s"%(b1,b2)
			print "- l0, b0 = %s, %s (deg)"%(skmp_deg_lon,skmp_deg_lat)

			print (msc_ref_lon-(skmp_ref_lon-band_lon_deg/2))/msc_del_lon,band_lon_deg/2
			#skymap[b1:b2,l1:l2] = msc_data

		print "Skymap shape: lon = %s, lat = %s"%(lon_dim,lat_dim)
		exit(0)	
		msc1_array = array(mosaicList[0][0]) # MV1
		msc2_array = array(mosaicList[1][0]) # MV2
		msc3_array = array(mosaicList[2][0]) # MW1
		msc4_array = array(mosaicList[3][0]) # MW2
		
		# Find index number of central position
		#gla = int(round(msc_ind_lon-1.))+(glon-msc_ref_lon)/msc_del_lon))
		#glb = int(round(msc_ind_lat-1.))+(glat-msc_ref_lat)/msc_del_lat))
		
		#b_cpt = msc_lat-db # lat conjunction point
		#l_cpt = msc_lon-dl # lon conjunction point
		#b_band = 1	   # thickness of the latitudinal band where averaging
		#l_band = 1	   # thickness of the longitudinal band where averaging
		
		#print b_cpt,l_cpt,lat_dim,lon_dim
		#print msc1_array[0:b_cpt,0:lon_dim/2.].shape
		#print msc3_array[0:b_cpt,0:lon_dim/2.].shape

		#skymap[0:b_cpt,0:lon_dim/2.] = msc1_array[0:b_cpt,0:lon_dim/2.]
		#skymap[0:b_cpt,lon_dim/2.:lon_dim] = msc3_array[0:b_cpt,0:lon_dim/2.]

		#skymap[b_cpt:lat_dim,0:lon_dim/2.] = msc2_array[0:msc_lat,0:lon_dim/2.]
		#skymap[b_cpt:lat_dim,lon_dim/2.:lon_dim] = msc4_array[0:msc_lat,0:lon_dim/2.]

		#b_band1 = (msc1_array[b_cpt-b_band:b_cpt+b_band,0:lon_dim/2.]+msc2_array[0:2*b_band,0:lon_dim/2.])/2.
		#b_band2 = (msc3_array[b_cpt-b_band:b_cpt+b_band,0:lon_dim/2.]+msc4_array[0:2*b_band,0:lon_dim/2.])/2.
		
		#l_band1 = (msc1_array[0:lat_dim/2.,l_cpt-l_band:l_cpt+l_band]+msc2_array[0:lat_dim/2.,0:2*l_band])/2.
		#l_band2 = (msc3_array[0:lat_dim/2.,l_cpt-l_band:l_cpt+l_band]+msc4_array[0:lat_dim/2.,0:2*l_band])/2.


		
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
		#results.writeto('skymap_test.fits', output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")

		
