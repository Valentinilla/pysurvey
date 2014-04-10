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
			flag1 = 'HI_unabsorbed_column_density'
		if self.species == 'CO':
			flag1 = 'WCO_intensity_line'
		
		if dim == '2D':
			flag2 = '2D'
		elif dim == '3D':
			flag2 = 'rings'
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_CombineMosaics')

		list = []
		if self.species == 'HI':
			list1 = ['MV1','MV2','MW1','MW2','MX1','MX2','MY1','MY2','MA1','MA2','MB1','MB2']
			list2 = ['MC1','MC2','MD1','MD2','ME1','ME2','MF1','MF2','MG1','MG2','MH1','MH2']
			list3 = ['MJ1','MJ2','MK1','MK2','ML1','ML2','MM1','MM2','MN1','MN2','MO1','MO2']
			list = list1+list2+list3
		if self.species == 'CO':
			# The 20 CO-mosaics
			list1 = ['MW1','MW2','MX1','MX2','MY1','MY2','MA1','MA2','MB1','MB2']
			list2 = ['MC1','MC2','MD1','MD2','ME1','ME2','MF1','MF2','MG1','MG2']
			list = list1+list2
		
		first_mosaic = list[0]
		n_msc = len(list)
		sqrt_n_msc = sqrt(n_msc)
		
		if not type==glob_N:
			self.logger.critical("Allowed type is: '"+glob_N+"'. Your entry is: '"+self.type+"'.")
			sys.exit(0)
		
		path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower()+'_column_density')
		ref_mosaic = path+self.survey+'_'+first_mosaic+'_'+flag1+'_'+flag2+'.fits'
		checkForFiles(self.logger,[ref_mosaic])
		
		f = pyfits.open(ref_mosaic)
		hdu1 = f[0]
		msc_size = hdu1.header['NAXIS']	# 4
		msc_x = hdu1.header['NAXIS1']	# 1024
		msc_y = hdu1.header['NAXIS2']	# 1024
		
		if dim == '3D':
			msc_z = hdu1.header['NAXIS3']
			#msc_rotz = hdu1.keyword["crota3"]
		
		msc_dx = float(hdu1.header['CDELT1']) # -4.9e-3 deg
		msc_rotx = float(hdu1.header['CROTA1']) # 0.0
		msc_dy = float(hdu1.header['CDELT2']) # 4.9e-3 deg
		msc_roty = float(hdu1.header['CROTA2']) # 0.0

		msc_bunit = hdu1.header['bunit']

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
			overlap_lon_px = round(1.12/msc_dy) # mosaics overlap by 1.12 deg = 224 px
			overlap_lat_px = round(1.12/msc_dy) # if dx = 0.005
			# needed for indexes
			odx = int(overlap_lon_px/2)
			ody = int(overlap_lat_px/2)
			
		if self.species == 'CO':
			overlap_lon_px = round(1.12/msc_dy) # mosaics overlap by 1.12 deg = 224 px
			overlap_lat_px = round(1.12/msc_dy) #
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

		# Using Multiprocessing if enough cpus are available
		#import multiprocessing
		#import itertools
		#ncpu = 2
		#self.logger.info("Running on %i cpu(s)"%(ncpu))
		#vec = [dim,odx,msc_x]
		#samples_list = (list1,list2)
		#pool = multiprocessing.Pool(processes=ncpu)
		#c1,c2 = pool.map(concatenateMosaics, itertools.izip(samples_list, itertools.repeat(vec)))
		#pool.close()
		#pool.join()
		#del samples_list
		#del results
		#exit(0)

		vec = [dim,odx,msc_x]
		if n_msc > 2:
			# Concatenate the lowest mosaics along the longitude axis 
			# |o|o|o|o|... = |o|o|o|o|
			# |x|x|x|x|... = |   x   |
			c1 = concatenateMosaics( (list1,vec) )
		
			# Concatenate the upmost mosaics along the latitude axis
			# |x|x|x|x|... = |   x   |
			# |x|x|x|x|... = |   x   |
			c2 = concatenateMosaics( (list2,vec) )
		else:
			c1 = hdu1.data
			c2 = hdu.data

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
			newheader["ctype3"] = ("Rband","Coordinate type")
			newheader["crval3"] = (0,"Ring of reference pixel")
			newheader["crpix3"] = (1.0,"Reference pixel of ring")
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
		
