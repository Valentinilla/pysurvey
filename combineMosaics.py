#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class combineMosaics(object):
	
	def __init__(self,surveyConf,mosaicConf,species,type,dim):
		"""
		Allow to combine mosaics of column density in 2D or 3D map.
		species = HI, WCO
		"""
		self.survey = surveyConf['survey']
		self.mosaic = mosaicConf['mosaic']
		self.species = species
		self.type = type
		sur = (self.survey).lower()

		flag1, flag2 = '',''
		
		if self.species == 'HI':
			flag1 = 'HI_unabsorbed'
		if self.species == 'CO':
			flag1 = 'WCO'
		
		if dim == '2D':
			flag2 = 'column_density'
		elif dim == '3D':
			flag2 = 'rings'
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_CombineMosaics')
		
		path = getPath(self.logger, key='list_mosaic')
		mosaiclist = path+self.species+'_mosaic.list'
		checkForFiles(self.logger,[mosaiclist])
		input = open(mosaiclist,"r")
		list = input.readlines()
		first_mosaic = list[0].split('\n')[0]
		n_msc = len(list)
		sqrt_n_msc = sqrt(n_msc)
		
		if not type==glob_N:
			self.logger.critical("Allowed type is: '"+glob_N+"'. Your entry is: '"+self.type+"'.")
			sys.exit(0)
		
		path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower()+'_column_density')
		ref_mosaic = path+self.survey+'_'+first_mosaic+'_'+flag1+'_'+flag2+'.fits'
		checkForFiles(self.logger,[ref_mosaic])
		
		f = pyfits.open(ref_mosaic)
		hdu = f[0]
		msc_size = hdu.header['NAXIS']	# 4
		msc_x = hdu.header['NAXIS1']	# 1024
		msc_y = hdu.header['NAXIS2']	# 1024
		
		if dim == '3D':
			msc_z = hdu.header['NAXIS3']
			#msc_rotz = hdu.keyword["crota3"]
		
		msc_dx = float(hdu.header['CDELT1']) # -4.9e-3 deg
		msc_rotx = float(hdu.header['CROTA1']) # 0.0
		msc_dy = float(hdu.header['CDELT2']) # 4.9e-3 deg
		msc_roty = float(hdu.header['CROTA2']) # 0.0

		msc_bunit = hdu.header['bunit']

		self.logger.info("Number of mosaics: %s"%n_msc)
		
		# TODO: Sort mosaics according to their coordinates		

		list1,list2 = [],[]
		for m in xrange(n_msc):
			msc_name = list[m].split('\n')[0]
			mosaic = path+self.survey+'_'+msc_name+'_'+flag1+'_'+flag2+'.fits'
			checkForFiles(self.logger,[mosaic])
			
			msc_file = pyfits.open(mosaic)
			hdu = msc_file[0]
			
			# select mosaics according to their ID number: 1 down, 2 up
			# |2|2|2|2|...
			# |1|1|1|1|...
			num = re.findall(r'\d',hdu.header['object'])
			if num[0] == '1': list1.append(hdu)
			if num[0] == '2': list2.append(hdu)			
						
		if self.species == 'HI':
			overlap_lon_px = 224 # mosaics overlap by 1.12 deg = 224 px
			overlap_lat_px = 224 #
			# needed for indexes
			odx = int(overlap_lon_px/2)
			ody = int(overlap_lat_px/2)
			
		if self.species == 'CO':
			overlap_lon_px = 224 # mosaics overlap by 1.12 deg = 224 px
			overlap_lat_px = 224 #
			# needed for indexes
			odx = int(overlap_lon_px/2)
			ody = int(overlap_lat_px/2)		

		if msc_size == 2:
			nx = msc_x*(n_msc/2) - overlap_lon_px*((n_msc/2)-1)
			ny = 2*msc_y - overlap_lat_px		
			skymap = zeros((ny,nx))
			
		if msc_size == 3:
			nx = msc_x*(n_msc/2) - overlap_lon_px*((n_msc/2)-1)
			ny = 2*msc_y - overlap_lat_px
			nz = msc_z		
			skymap = zeros((nz,ny,nx))

		# Concatenate the lowest mosaics along the longitude axis 
		# |o|o|o|o|... = |o|o|o|o|
		# |x|x|x|x|... = |   x   |
		c1 = concatenateMosaics(list1,dim,odx)
		
		# Concatenate the upmost mosaics along the latitude axis
		# |x|x|x|x|... = |   x   |
		# |x|x|x|x|... = |   x   |
		c2 = concatenateMosaics(list2,dim,odx)
			
		# Concatenate the two raw of mosaics along the latitude axis
		# |   x   |... = |   x   |
		# |   x   |...   |       |
		if dim == '2D':
			skymap = concatenate( (c1[:-ody,:],c2[ody:,:]), axis=0)
		if dim == '3D':		
			skymap = concatenate( (c1[:,:-ody,:],c2[:,ody:,:]), axis=1)

		# Header keys
		# CRVAL1
		crv1_msc1 = list1[((n_msc/4)-1)].header['crval1']
		crv1_msc2 = list1[((n_msc/4))].header['crval1']
		crval1 = (crv1_msc1 + crv1_msc2)/2.
		# CRVAL2
		crv2_msc1 = list1[0].header['crval2'] # -1
		crv2_msc2 = list2[0].header['crval2'] # 3
		crval2 = (crv2_msc1 + crv2_msc2)/2.			

		# CRPIXN
		if dim == '2D':
			crpix1,crpix2 = round(skymap.shape[1]/2.),round(skymap.shape[0]/2.)
		if dim == '3D':
			crpix1,crpix2 = round(skymap.shape[2]/2.),round(skymap.shape[1]/2.)
		
		lonsign = getSign(crval1,string=True)
		latsign = getSign(crval2,string=True)
		
		# Store results
		newheader = pyfits.Header()
		newheader["ctype1"] = ("GLON-CAR","Coordinate type")
		newheader["crval1"] = (crval1,"Galactic longitude of reference pixel")
		newheader["crpix1"] = (crpix1,"Reference pixel of lon")
		newheader["cdelt1"] = (msc_dx,"Longitude increment")
		newheader["crota1"] = (msc_rotx,"Longitude rotation")
		newheader["cunit1"] = ("deg","Unit type")
		
		newheader["ctype2"] = ("GLAT-CAR","Coordinate type")
		newheader["crval2"] = (crval2,"Galactic latitude of reference pixel")
		newheader["crpix2"] = (crpix2,"Reference pixel of lat")
		newheader["cdelt2"] = (msc_dy,"Latitude increment")
		newheader["crota2"] = (msc_roty,"Latitude rotation")
		newheader["cunit2"] = ("deg","Unit type")

		if dim == '3D':
			newheader["ctype3"] = ("Ring","Coordinate type")
			newheader["crval3"] = (1,"Ring of reference pixel")
			newheader["crpix3"] = (1,"Reference pixel of ring")
			newheader["cdelt3"] = (1,"Ring increment")
			#newheader["crota3"] = (msc_rotz,"Ring rotation")

		newheader["bunit"] = (msc_bunit,"Map units")
		newheader["datamin"] = (amin(skymap),"Min value")
		newheader["datamax"] = (amax(skymap),"Max value")
		newheader["object"] = (self.survey+" Skymap",self.survey+" Mosaic")
		
		results = pyfits.PrimaryHDU(skymap, newheader)
		#results.scale('int16', '', bscale=1, bzero=32768)
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower())
		skymap_name = path+'%s_G%.2f%s%.2f'%(self.survey,crval1,latsign,crval2)
		results.writeto(skymap_name+'.fits', output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
