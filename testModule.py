#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*

class testModule(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf,scale_data=False):
		"""
		Calculate the column density in galacto-centric rings using the rotation curve of M. Phol et al.
		Boundaries of galacto-centric annuli by M.Ackermann et al.
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.newspec #specified by the user
		self.type = mosaic.type
		self.datatype = mosaic.datatype

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_TestModule')
		file,path,flag,units = '','','',''
		sur = (self.survey).lower()
		
		if self.species == 'HI':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			if mosaic.totmsc == 1:
				flag = 'HI_unabsorbed_column_density'
			else:
				flag = 'HI_unabsorbed_column_density_part_%i-%i'%(int(mosaic.nmsc),int(mosaic.totmsc))
			units = '10e+20 H atoms cm-2'
		elif self.species == 'HISA':
			path = getPath(self.logger,'lustre_'+sur+'_hisa_column_density')
			flag = 'HISA_column_density'
			units = '10e+20 H atoms cm-2'
		elif self.species == 'CO':
			path = getPath(self.logger,'lustre_'+sur+'_co_column_density')
			flag = 'WCO_intensity_line'
			units = 'K km s-1'
		
		file = path+self.survey+'_'+self.mosaic+'_'+flag+'_rings.fits'
		checkForFiles(self.logger,[file],existence=True)
		#print mosaic.filename		

		self.logger.info("Open file and get data...")
				
		# Get HI emission data
		Tb = mosaic.observation[0,:,:,:]
		
		lon = mosaic.xarray
		lat = mosaic.yarray
		
		vel = mosaic.zarray/1000.
		dv = fabs(mosaic.dz/1000.)		# [velocity] = km s-1
		
		nlon = mosaic.nx
		nlat = mosaic.ny
		nvel = mosaic.nz

		# free memory
		del mosaic.observation
		del mosaic.xarray
		del mosaic.yarray
		del mosaic.zarray
		
		Ts = float(utilsConf['tspin'])	  	# [Excitation (or Spin) Temperature] = K (150)
		C = float(utilsConf['c'])	  	# [costant] = cm-2
		
		rmin,rmax,annuli = getAnnuli()
		
		# Array to store results
		cubemap = zeros((annuli,nlat,nlon),dtype=float32)
		
		self.logger.info("Initializing parameters...")
		self.logger.info("1) Ts = %.2f K"%Ts)
		self.logger.info("2) dV = %.2f km/s"%dv)
		self.logger.info("3) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
					
		self.logger.info("Calculating gas distribution...")
		path2 = getPath(self.logger,'rotcurve_mpohl')
		list = [self.species,vel,lat,lon,dv,path2,C,Ts,rmin,rmax]#,self.logger]		
		
		# Using Multiprocessing if enough cpus are available
		import multiprocessing
		
		ncpu = glob_ncpu#int(ceil(multiprocessing.cpu_count()/1))
		if ncpu > 16: ncpu = 16
		#print float(Tb.shape[1])/ncpu
		#print Tb.shape
		
		if Tb.shape[1] < ncpu:
			ncpu = Tb.shape[1]
		
		self.logger.info("Running on %i cpu(s)"%(ncpu))
		if ncpu > 1:
			import itertools		
			samples_list = array_split(Tb, ncpu, axis = 1)
			pool = multiprocessing.Pool(processes = ncpu)
			results = pool.map(rotation_curve, itertools.izip(samples_list, itertools.repeat(list)))
			pool.close()
			pool.join()
			cubemap = concatenate(results,axis=1)
			del samples_list
			del list
			del results
		else:
			cubemap = rotation_curve( (Tb,list) )
			
		if self.species == 'HI' or self.species == 'HISA':
			cubemap = cubemap*1e-20

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
    		
		newheader["ctype3"] = ("Rband","Coordinate type")
		newheader["crval3"] = (0,"Ring of reference pixel")
		newheader["crpix3"] = (1.0,"Reference pixel of ring")
		newheader["cdelt3"] = (1,"Ring increment")
		newheader["crota3"] = (mosaic.keyword["crota3"],"Ring rotation")
		
		newheader["bunit"] = (units,"Map units")
		newheader["datamin"] = (amin(cubemap),"Min value")
		newheader["datamax"] = (amax(cubemap),"Max value")

		newheader['minfil'] = unravel_index(argmin(cubemap),cubemap.shape)[0]
		newheader['mincol'] = unravel_index(argmin(cubemap),cubemap.shape)[1]
		newheader['minrow'] = unravel_index(argmin(cubemap),cubemap.shape)[2]
		newheader['maxfil'] = unravel_index(argmax(cubemap),cubemap.shape)[0]
		newheader['maxcol'] = unravel_index(argmax(cubemap),cubemap.shape)[1]
		newheader['maxrow'] = unravel_index(argmax(cubemap),cubemap.shape)[2]

		newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		
		# Output file
		results = pyfits.PrimaryHDU(cubemap, newheader)
		if scale_data:
			self.logger.info("Writing scaled data to a fits file in...")
    			results.scale('int16', '', bscale=mosaic.bscale, bzero=mosaic.bzero)
		else:
			self.logger.info("Writing data to a fits file in...")

		# Create a Table with the annuli boundaries
		col1 = pyfits.Column(name='Rmin', format='1E', unit='kpc', array=array(rmin))
		col2 = pyfits.Column(name='Rmax', format='1E', unit='kpc', array=array(rmax))
		cols = pyfits.ColDefs([col1,col2])
		tbl = pyfits.new_table(cols)
		tbl.name = "BINS"
		
		thdulist = pyfits.HDUList([results,tbl])
		
		thdulist.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
