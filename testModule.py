#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*

class testModule(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf):
		"""
		Calculate the column density in galacto-centric rings using the rotation curve of M. Phol et al.
		Boundaries of galacto-centric annuli by M.Ackermann et al.
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.newspec #specified by the user
		self.type = mosaic.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_TestModule')
		file,path,flag,units = '','','',''
		sur = (self.survey).lower()
		
		if self.species == 'HI':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			flag = 'HI_unabsorbed'
			if sur == 'lab':
				flag = 'HI'
			units = '10e+20 H atoms cm-2'
		elif self.species == 'HISA':
			path = getPath(self.logger,'lustre_'+sur+'_hisa_column_density')
			flag = 'HISA'
			units = '10e+20 H atoms cm-2'
		elif self.species == 'CO':
			path = getPath(self.logger,'lustre_'+sur+'_co_column_density')
			flag = 'WCO'
			units = 'K km s-1'
		
		file = path+self.survey+'_'+self.mosaic+'_'+flag+'_rings.fits'
		checkForFiles(self.logger,[file],existence=True)
		print mosaic.filename
		
		self.logger.info("Open file and get data...")
				
		# Get HI emission data
		Tb = mosaic.observation[0,:,:,:]
		print amin(Tb),amax(Tb)
		exit(0)
		lon = mosaic.xarray
		lat = mosaic.yarray
		vel = mosaic.zarray/1000.
		
		nlon = mosaic.nx
		nlat = mosaic.ny
		nvel = mosaic.nz

		# free memory
		del mosaic.observation
		del mosaic.xarray
		del mosaic.yarray
		del mosaic.zarray
		
		Ts = float(utilsConf['tspin'])	  	# [Excitation (or Spin) Temperature] = K (150)
		Tbg = float(utilsConf['tcmb'])	  	# [Cosmic Microwave Background (CMB)] = K
		if not self.species == 'WCO':
			dv = fabs(mosaic.dz/1000.)	# [velocity] = km s-1
		C = float(utilsConf['c'])	  	# [costant] = cm-2
		
		# Corrections for latitude (almost negligible)
		cosdec = cos(radians(lat)) 		# [lat. correction] = rad

		# Boundaries of Galactocentric Annuli - M.Ackermann et al
		# (The Astrophysical Journal, 750:3-38, 2012)
		annuli = 17
		rmin = [0.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.]
		rmax = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.,50.]
		#ann_boundaries = [ [i,rmin[i],rmax[i]] for i in xrange(annuli)]

		
		# Array to store results
		cubemap = zeros((annuli,nlat,nlon),dtype=float)			
		#NHI = zeros((annuli,nlat,nlon),dtype=float)
		#ITb = zeros((annuli,nlat,nlon),dtype=float)	
		
		self.logger.info("Initializing parameters...")
		self.logger.info("1) Ts = %.2f K"%Ts)
		self.logger.info("2) dV = %.2f km/s"%dv)
		self.logger.info("3) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
					
		self.logger.info("Calculating gas distribution...")
		path2 = getPath(self.logger,'rotcurve_mpohl')
		list = [self.species,vel,lat,lon,dv,path2,C,Ts,rmin,rmax]#,self.logger]		
		
		# Using Multiprocessing if enough cpus are available
		import multiprocessing
		import itertools
		
		ncpu = int(ceil(multiprocessing.cpu_count()/1))
		self.logger.info("Running on %i cpu(s)"%(ncpu))
		
		samples_list = array_split(Tb, ncpu, axis = 1)
		pool = multiprocessing.Pool(processes = ncpu)
		results = pool.map(rotation_curve, itertools.izip(samples_list, itertools.repeat(list)))
		
		del samples_list[:]
		del list[:]
		
		cubemap = concatenate(results,axis=1).astype(cubemap.dtype)
		if self.species == 'HI' or self.species == 'HISA':
			cubemap = where(cubemap>0.,cubemap*1e-20,cubemap)
		#cubemap = rotation_curve( (Tb,list) )

		# Store results
		newheader = pyfits.Header()
		newheader["ctype1"] = ("GLON-CAR","Coordinate type")
		newheader["crval1"] = (mosaic.keyword["crval1"],"Galactic longitude of reference pixel")
		newheader["crpix1"] = (mosaic.keyword["crpix1"],"Reference pixel of lon")
		newheader["cdelt1"] = (mosaic.keyword["cdelt1"],"Longitude increment")
		newheader["crota1"] = (mosaic.keyword["crota1"],"Longitude rotation")
		newheader["cunit1"] = ("deg","Unit type")
				
		newheader["ctype2"] = ("GLAT-CAR","Coordinate type")
		newheader["crval2"] = (mosaic.keyword["crval2"],"Galactic latitude of reference pixel")
		newheader["crpix2"] = (mosaic.keyword["crpix2"],"Reference pixel of lat")
		newheader["cdelt2"] = (mosaic.keyword["cdelt2"],"Latitude increment")
		newheader["crota2"] = (mosaic.keyword["crota2"],"Latitude rotation")
		newheader["cunit2"] = ("deg","Unit type")
    		
		newheader["ctype3"] = ("Ring","Coordinate type")
		newheader["crval3"] = (1,"Ring of reference pixel")
		newheader["crpix3"] = (1,"Reference pixel of ring")
		newheader["cdelt3"] = (1,"Ring increment")
		newheader["crota3"] = (mosaic.keyword["crota3"],"Ring rotation")
		
		newheader["bunit"] = (units,"Map units")
		newheader["datamin"] = ("%e"%amin(cubemap),"Min value")
		newheader["datamax"] = ("%e"%amax(cubemap),"Max value")
		newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		
		results = pyfits.PrimaryHDU(cubemap, newheader)
    		#results.scale('int16', '', bscale=mosaic.bscale, bzero=mosaic.bzero)
		
		# Create a Table with the annuli boundaries
		col1 = pyfits.Column(name='Rmin', format='1E', unit='kpc', array=array(rmin))
		col2 = pyfits.Column(name='Rmax', format='1E', unit='kpc', array=array(rmax))
		cols = pyfits.ColDefs([col1,col2])
		tbl = pyfits.new_table(cols)
		tbl.name = "BINS"
		
		thdulist = pyfits.HDUList([results,tbl])
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		thdulist.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
