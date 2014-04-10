#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class combineMosaics(object):
	
	def __init__(self,surveyConf,mosaic,species,type,dim):
		"""
		Allow to combine mosaics of column density in 2D or 3D map.
		species = HI, HI_unabsorbed, HISA, WCO
		"""
		self.survey = surveyConf['survey']
		self.mosaic = mosaic#mosaicConf['mosaic']
		self.species = species
		self.type = type
		sur = self.survey.lower()
		
		flag1, flag2 = '',''
		HI_condition = (self.species == 'HI' or self.species == 'HI_unabsorbed' or self.species == 'HISA')
		
		if HI_condition:
			flag1 = self.species+'_column_density'
		elif self.species == 'CO':
			flag1 = 'WCO_intensity_line'
		
		if dim == '2D':
			flag2 = '2D'
		elif dim == '3D':
			flag2 = 'rings'
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_CombineMosaics')

		if not type==glob_N:
			self.logger.critical("Allowed type is: '"+glob_N+"'. Your entry is: '"+self.type+"'.")
			sys.exit(0)		

		# Read the list of files
		path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower()+'_column_density')
		list = []
		if self.mosaic == 'skymap':
			if self.survey == 'CGPS':
				if HI_condition:
					list1 = ['MV1','MV2','MW1','MW2','MX1','MX2','MY1','MY2','MA1','MA2','MB1','MB2']
					list2 = ['MC1','MC2','MD1','MD2','ME1','ME2','MF1','MF2','MG1','MG2','MH1','MH2']
					list3 = ['MIJ1','MIJ2','MK1','MK2','ML1','ML2','MM1','MM2','MN1','MN2','MO1','MO2']
					tlist = list1+list2+list3
	
				elif self.species == 'CO':
					# The 20 CO-mosaics
					list1 = ['MW1','MW2','MX1','MX2','MY1','MY2','MA1','MA2','MB1','MB2']
					list2 = ['MC1','MC2','MD1','MD2','ME1','ME2','MF1','MF2','MG1','MG2']
					tlist = list1+list2
					
				# this configuration (for m --> for f) sort lists in the right way
				for m in tlist:
					for f in os.listdir(path):
						if f.startswith('CGPS_%s_%s'%(m,flag1)):
							list.append(f)
		else:
			if self.survey == 'CGPS':
				if HI_condition:
					list = [f for f in os.listdir(path) if f.startswith('CGPS_%s_%s_part'%(self.mosaic,flag1))]
					list = sort(list)

		# Total number of mosaics
		n_msc = len(list)
		ref_mosaic = path+list[0]
		checkForFiles(self.logger,[ref_mosaic])
		
		f = pyfits.open(ref_mosaic)
		hdu1 = f[0]
		msc_size = hdu1.header['naxis']	# 4
		msc_x = hdu1.header['naxis1']	# 1024
		msc_y = hdu1.header['naxis2']	# 1024
		
		if dim == '3D':
			msc_z = hdu1.header['naxis3']
			#msc_rotz = hdu1.keyword["crota3"]
		
		msc_dx = float(hdu1.header['cdelt1']) # -4.9e-3 deg
		msc_rotx = float(hdu1.header['crota1']) # 0.0
		msc_dy = float(hdu1.header['cdelt2']) # 4.9e-3 deg
		msc_roty = float(hdu1.header['crota2']) # 0.0

		msc_bunit = hdu1.header['bunit']

		self.logger.info("Number of mosaics: %s"%n_msc)
		
		# TODO: Sort mosaics according to their coordinates		

		list1,list2 = [],[]
		for m in xrange(n_msc):
			mosaic = path+list[m]
			checkForFiles(self.logger,[mosaic])
			
			msc_file = pyfits.open(mosaic)
			hdu = msc_file[0]
			
			if self.mosaic == 'skymap':
				# select mosaics according to their ID number: 1 down, 2 up
				# |2|2|2|2|...
				# |1|1|1|1|...
				num = re.findall(r'\d',hdu.header['object'])
				if num[0] == '1': list1.append(hdu)
				if num[0] == '2': list2.append(hdu)			
			else:
				list1.append(hdu)
		
		if self.species == 'HI':
			if self.mosaic == 'skymap':
				overlap_lon_px = round(1.12/msc_dy) # mosaics overlap by 1.12 deg = 224 px
				overlap_lat_px = round(1.12/msc_dy) # if dx = 0.005
				# needed for indexes
				odx = int(overlap_lon_px/2)
				ody = int(overlap_lat_px/2)
			else:
				odx,ody = 0,0
			
		if self.species == 'CO':
			if self.mosaic == 'skymap':
				overlap_lon_px = round(1.12/msc_dy) # mosaics overlap by 1.12 deg = 224 px
				overlap_lat_px = round(1.12/msc_dy) #
				# needed for indexes
				odx = int(overlap_lon_px/2)
				ody = int(overlap_lat_px/2)		
			else:
				odx,ody = 0,0

		if msc_size == 2:
			nx = msc_x*(n_msc/2) - overlap_lon_px*((n_msc/2)-1)
			ny = 2*msc_y - overlap_lat_px		
			skymap = zeros((ny,nx))
			
		if msc_size == 3:
			if self.mosaic == 'skymap':
				nx = msc_x*(n_msc/2) - overlap_lon_px*((n_msc/2)-1)
				ny = 2*msc_y - overlap_lat_px
			else:
				nx = list1[0].header['naxis1']
				ny = 0
				for m in xrange(len(list1)):
					ny+=list1[m].header['naxis2']
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
		if self.mosaic == 'skymap':
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
			
		else:
			c = []
			index = 0
			for current, next in zip(list1, list1[1:]):
				if index == 0:
					c = concatenate((current.data[:,:,:], next.data[:,:,:]), axis=1)
				elif index > 0 and index < len(list):
					c = concatenate((c[:,:,:], next.data[:,:,:]), axis=1)
				elif index == len(list):
					c = concatenate((c, next.data[:,:,:]), axis=1)
				index += 1
			
			skymap = array(c)
			
			crval1 = list1[0].header['crval1']

			cos1 = 0
			sin1 = 0
			for m in xrange(len(list1)):
				cos1 += cos(list1[m].header['crval2'])
				sin1 += sin(list1[m].header['crval2'])
			crval2 = arctan2(sin1,cos1)-msc_dx
		
		# CRPIXN
		if dim == '2D':
			crpix1,crpix2 = round(skymap.shape[1]/2.),round(skymap.shape[0]/2.)
		if dim == '3D':
			crpix1,crpix2 = round(skymap.shape[2]/2.)+1,round(skymap.shape[1]/2.)+1

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
		if self.mosaic == 'skymap':
			newheader["object"] = (self.survey+" Skymap",self.survey+" Mosaic")
		else:
			newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
			newheader['minfil'] = unravel_index(argmin(skymap),skymap.shape)[0]
			newheader['mincol'] = unravel_index(argmin(skymap),skymap.shape)[1]
			newheader['minrow'] = unravel_index(argmin(skymap),skymap.shape)[2]
			newheader['maxfil'] = unravel_index(argmax(skymap),skymap.shape)[0]
			newheader['maxcol'] = unravel_index(argmax(skymap),skymap.shape)[1]
			newheader['maxrow'] = unravel_index(argmax(skymap),skymap.shape)[2]

		results = pyfits.PrimaryHDU(skymap, newheader)
				
		# Output file
		self.logger.info("Writing data to a fits file in...")
		skymap_name = ''
		path,key = '',''
		if self.mosaic == 'skymap':
			path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower())
			skymap_name = path+'%s_G%.2f%s%.2f.fits'%(self.survey,crval1,latsign,crval2)
		else:
			path = getPath(self.logger, key='lustre_'+sur+'_'+self.species.lower()+'_column_density')
			if HI_condition:
				key = self.species
				flag = key+'_column_density'
			elif self.species == 'CO':
				key = 'WCO'
				flag = key+'_intensity_line'
			dir = path+self.survey+'_'+self.mosaic+'_'+key
			os.system('mkdir %s'%dir)
			os.system('mv %s* %s'%(path+self.survey+'_'+self.mosaic+'_'+flag+'_part_',dir))
			skymap_name = path+self.survey+'_'+self.mosaic+'_'+flag+'_rings.fits'

		rmin,rmax,annuli = getAnnuli(glob_annuli)
		
		# Create a Table with the annuli boundaries
		col1 = pyfits.Column(name='Rmin', format='1E', unit='kpc', array=array(rmin))
		col2 = pyfits.Column(name='Rmax', format='1E', unit='kpc', array=array(rmax))
		cols = pyfits.ColDefs([col1,col2])
		tbl = pyfits.new_table(cols)
		tbl.name = "BINS"
			
		thdulist = pyfits.HDUList([results,tbl])			
		thdulist.writeto(skymap_name, output_verify='fix')

		self.logger.info("%s"%path)
		self.logger.info("Done")

		
