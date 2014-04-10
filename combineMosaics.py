#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class combineMosaics(object):
	
	def __init__(self,surveyConf,mosaicConf,species,type,flag):
		"""
		Allow to combine mosaics of column density in a 2D map.
		"""
		self.survey = surveyConf['survey']
		self.mosaic = mosaicConf['mosaic']
		self.species = species
		self.type = type
		sur = (self.survey).lower()
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_CombineMosaics')
		
		path = '/afs/ifh.de/group/that/work-sf/survey/mysurvey/lists/'
		mosaiclist = path+flag+'.list'
		checkForFiles(self.logger,[mosaiclist])
		input = open(mosaiclist,"r")
		list = input.readlines()
		first_mosaic = list[0].split('\n')[0]
		mosaics_number = len(list)
		sqrt_msc_num = sqrt(mosaics_number)
		
		if not type==glob_N:
			self.logger.critical("Allowed type is: '"+glob_N+"'. Your entry is: '"+self.type+"'.")
			os.exit(0)
					
		path = getPath(self.logger, key='lustre_'+sur+'_hi_column_density')
		ref_mosaic = path+self.survey+'_'+first_mosaic+'_HI_unabsorbed_column_density.fits'
		checkForFiles(self.logger,[ref_mosaic])
		
		f = pyfits.open(ref_mosaic)
		hdu = f[0]
		#msc_size = hdu.header['NAXIS']	# 4
		msc_lon = hdu.header['NAXIS1']	# 1024
		msc_lat = hdu.header['NAXIS2']	# 1024
		#msc_ref_lon = float(hdu.header['CRVAL1']) # deg (GLON-CAR)
		msc_del_lon = float(hdu.header['CDELT1']) # -4.9e-3 deg
		#msc_ind_lon = float(hdu.header['CRPIX1']) # 513 px
		#msc_rot_lon = float(hdu.header['CROTA1']) # 0.0
		#msc_ref_lat = float(hdu.header['CRVAL2']) # deg
		msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
		msc_bunit = hdu.header['bunit']

		self.logger.info("Number of mosaics: %s"%mosaics_number)
		mosaicList = []
		list1,list2 = [],[]
		skymap = zeros((msc_lat,msc_lon))
		for m in range(0,mosaics_number):
			msc_name = list[m].split('\n')[0]
			mosaic = path+self.survey+'_'+msc_name+'_HI_unabsorbed_column_density.fits'
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
			mosaicList.append([msc_data,m])
			
			list1.append([msc_ind_lon,msc_ind_lat])
			list2.append([msc_ref_lon,msc_ref_lat])
			
			if msc_size > 2:
				self.logger.critical("Mosaic dimension > 2!!")
				os.exit(0)
						
			#print "Mosaic n %s (%s)"%(m,msc_name)
			#print "- lon0 = %s, lat0 = %s"%(msc_ind_lon,msc_ind_lat)
			#print "- lon1 = %s, lon2 = %s"%(msc_ind_lon-msc_ind_lon,2*msc_ind_lon)
			#print "- lat1 = %s, lat2 = %s"%(msc_ind_lat-msc_ind_lat,2*msc_ind_lat)
			#print "- l0, b0 = %s, %s (deg)"%(msc_ref_lon,msc_ref_lat)
			#print "- lon_del, lat_del = %s, %s"%(msc_del_lon,msc_del_lat)
		
		blon_px = 120 # mosaics overlap by 1.2 deg = 240 px
		blat_px = 120 #
			
		msc1 = array(mosaicList[0][0]) # MC1
		lon1 = 2*list1[0][0]
		lat1 = 2*list1[0][1]
		lon1_deg = list2[0][0]
		lat1_deg = list2[0][1]
		msc2 = array(mosaicList[1][0]) # MC2
		lon2 = 2*list1[1][0]
		lat2 = 2*list1[1][1]
		lon2_deg = list2[1][0]
		lat2_deg = list2[1][1]
		msc3 = array(mosaicList[2][0]) # MD1
		lon3 = 2*list1[2][0]
		lat3 = 2*list1[2][1]
		lon3_deg = list2[2][0]
		lat3_deg = list2[2][1]
		msc4 = array(mosaicList[3][0]) # MD2
		lon4 = 2*list1[3][0]
		lat4 = 2*list1[3][1]
		lon4_deg = list2[3][0]
		lat4_deg = list2[3][1]
		
		c1 = concatenate((msc1[:,:lon1-blon_px],msc3[:,blon_px:lon3]),axis=1)
		c2 = concatenate((msc2[:,:lon2-blon_px],msc4[:,blon_px:lon4]),axis=1)
		skymap = concatenate((c1[:lat1-blat_px,:],c2[blat_px:lat1,:]),axis=0)		
		# Header keys
		crpix1,crpix2 = round(skymap.shape[1]/2.),round(skymap.shape[0]/2.)
		crval1 = (lon1_deg+lon3_deg-2*blon_px*msc_del_lon)/2.
		crval2 = (lat1_deg+lat3_deg-2*blat_px*msc_del_lat)/2.
		lonsign = getSign(crval1,string=True)
		latsign = getSign(crval2,string=True)
		
		# Store results
		newheader = pyfits.Header()
		newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
		newheader.update(key="crval1", value=crval1, comment="Galactic longitude of reference pixel")
		newheader.update(key="crpix1", value=crpix1, comment="Reference pixel of lon")
		newheader.update(key="cdelt1", value=msc_del_lon, comment="Longitude increment")
		newheader.update(key="crota1", value=msc_rot_lon, comment="Longitude rotation")
		newheader.update(key="cunit1", value="deg", comment="Unit type")
		
		newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
		newheader.update(key="crval2", value=crval2, comment="Galactic latitude of reference pixel")
		newheader.update(key="crpix2", value=crpix2, comment="Reference pixel of lat")
		newheader.update(key="cdelt2", value=msc_del_lat, comment="Latitude increment")
		newheader.update(key="crota2", value=msc_rot_lat, comment="Latitude rotation")
		newheader.update(key="cunit2", value="deg", comment="Unit type")
		
		newheader.update(key="bunit", value=msc_bunit, comment="Map units")
		
		newheader.update(key="datamin", value=amin(skymap))
		newheader.update(key="datamax", value=amax(skymap))
		
		newheader.update(key="object", value=self.survey+" Skymap", comment=self.survey+" Mosaic")
		
		results = pyfits.PrimaryHDU(skymap, newheader)
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		skymap_name = '%s_G%.2f%s%.2f'%(self.survey,crval1,latsign,crval2)
		results.writeto(skymap_name+'.fits', output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
