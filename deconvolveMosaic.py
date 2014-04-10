#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*

class deconvolveMosaic(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf,rotcurve,scale_data=False):
		"""
		Calculate the column density in galacto-centric rings using the rotation curve of M. Phol et al.
		Boundaries of galacto-centric annuli by M.Ackermann et al.
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.newspec #specified by the user
		self.type = mosaic.type
		self.datatype = mosaic.datatype
		self.totmsc = mosaic.totmsc
		self.nmsc = mosaic.nmsc
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_Deconvolution')
		file,flag,units = '','',''
		sur = self.survey.lower()
		HI_all_OR = (self.species == 'HI' or self.species == 'HI_unabsorbed' or self.species == 'HISA')

		path = getPath(self.logger,'lustre_'+sur+'_'+self.species.lower()+'_column_density')		
		if HI_all_OR:
			if self.totmsc == 1:
				flag = self.species+'_column_density'
			else:
				flag = self.species+'_column_density_part_%s-%s'%(self.nmsc,self.totmsc)
			units = '10e+20 H atoms cm-2'
		elif self.species == 'CO':
			flag = 'WCO_intensity_line'
			units = 'K km s-1'
		
		file = path+self.survey+'_'+self.mosaic+'_'+flag+'_rings.fits'
		checkForFiles(self.logger,[file],existence=True)		
		
		self.logger.info("Open file and get data...")
		
		# Get HI emission data
		if self.survey == 'LAB':
			Tb = mosaic.observation[:,:,:]
		else:
			Tb = mosaic.observation[0,:,:,:]
		
		# In case of multiprocessing analysis split the array along the maximum axis
		maxis = 1+argmax(Tb.shape[-2:])
		
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
		
		rmin,rmax,annuli = getAnnuli(glob_annuli)
		if not (rotcurve == 'Bissantz2003' or rotcurve == 'Clemens1985'):
			self.logger.critical("You must enter a right rotation curve! Options are: 'Bissantz2003' or 'Clemens1985'")
			self.logger.critical("Your entry is %s"%rotcurve)
			sys.exit(0)
		
		# Array to store results
		cubemap = zeros((annuli,nlat,nlon),dtype=float32)
		
		self.logger.info("Initializing parameters...")
		self.logger.info("1) Ts = %.2f K"%Ts)
		self.logger.info("2) dV = %.2f km/s"%dv)
		self.logger.info("3) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
		self.logger.info("4) Rotation curve: '%s'"%rotcurve)
		self.logger.info("5) Annuli: '%s'"%glob_annuli)
		
		self.logger.info("Calculating gas distribution...")
		path2 = getPath(self.logger,'rotcurve_mpohl')
		
		list = []
		if maxis == 1:
			list = [self.species,lon,vel,dv,path2,C,Ts,rmin,rmax,rotcurve,maxis]
			coord = lat
		elif maxis == 2:
			list = [self.species,lat,vel,dv,path2,C,Ts,rmin,rmax,rotcurve,maxis]
			coord = lon
		else:
			self.logger.critical("ERROR in splitting Tb!")
			sys.exit(0)
		
		# Using Multiprocessing if enough cpus are available
		import multiprocessing
		
		ncpu = glob_ncpu
		# Maximum number of cpus
		if ncpu > 16: ncpu = 16
		# Minimum number of cpus
		if Tb.shape[maxis] < ncpu:
			ncpu = Tb.shape[maxis]

		#print coord.shape
		#print coord.shape[0]/ncpu
		#exit(0)		
	
		self.logger.info("Running on %i cpu(s)"%(ncpu))
		if ncpu > 1:
			import itertools		
			arrays = array_split(Tb, ncpu, axis=maxis)
			coords = array_split(coord, ncpu, axis=0)
			pool = multiprocessing.Pool(processes = ncpu)
			results = pool.map(Deconvolution, itertools.izip(arrays,coords,itertools.repeat(list)))
			pool.close()
			pool.join()
			cubemap = concatenate(results,axis=maxis)
			del arrays
			del coords
			del list
			del results
		else:
			cubemap = Deconvolution( (Tb,coord,list) )
			
		if HI_all_OR: cubemap = cubemap*1e-20
		
		# Store results
		newheader = pyfits.Header()
		newheader['ctype1'] = ("GLON-CAR","Coordinate type")
		newheader['crval1'] = (mosaic.keyword["crval1"],"Galactic longitude of reference pixel")
		newheader['crpix1'] = (mosaic.keyword["crpix1"],"Reference pixel of lon")
		newheader['cdelt1'] = (mosaic.keyword["cdelt1"],"Longitude increment")
		newheader['crota1'] = (mosaic.keyword["crota1"],"Longitude rotation")
		newheader['cunit1'] = ("deg","Unit type")
				
		newheader['ctype2'] = ("GLAT-CAR","Coordinate type")
		newheader['crval2'] = (mosaic.keyword["crval2"],"Galactic latitude of reference pixel")
		newheader['crpix2'] = (mosaic.keyword["crpix2"],"Reference pixel of lat")
		newheader['cdelt2'] = (mosaic.keyword["cdelt2"],"Latitude increment")
		newheader['crota2'] = (mosaic.keyword["crota2"],"Latitude rotation")
		newheader['cunit2'] = ("deg","Unit type")
    		
		newheader['ctype3'] = ("Rband","Coordinate type")
		newheader['crval3'] = (0,"Ring of reference pixel")
		newheader['crpix3'] = (1.0,"Reference pixel of ring")
		newheader['cdelt3'] = (1,"Ring increment")
		newheader['crota3'] = (mosaic.keyword["crota3"],"Ring rotation")
		
		newheader['bunit'] = (units,"Map units")
		newheader['datamin'] = (amin(cubemap),"Min value")
		newheader['datamax'] = (amax(cubemap),"Max value")

		newheader['minfil'] = unravel_index(argmin(cubemap),cubemap.shape)[0]
		newheader['mincol'] = unravel_index(argmin(cubemap),cubemap.shape)[1]
		newheader['minrow'] = unravel_index(argmin(cubemap),cubemap.shape)[2]
		newheader['maxfil'] = unravel_index(argmax(cubemap),cubemap.shape)[0]
		newheader['maxcol'] = unravel_index(argmax(cubemap),cubemap.shape)[1]
		newheader['maxrow'] = unravel_index(argmax(cubemap),cubemap.shape)[2]
		if self.totmsc == 1:
			newheader['object'] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		else:
			newheader['object'] = ("Mosaic %s (%s/%s)"%(self.mosaic,self.nmsc,self.totmsc),"%s Mosaic (n/tot)"%self.survey)
		newheader.add_history('Rotation curve: %s'%rotcurve)
		newheader.add_history('Annuli: %s'%glob_annuli)

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
		
