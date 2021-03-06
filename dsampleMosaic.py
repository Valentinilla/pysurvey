#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class dsampleMosaic(object):
	
	def __init__(self,obs,scale):
		"""
		Downsample mosaics: scale x dx_msc = new_dx_msc
		Input can be either 'mosaic' or 'skymap'.
		"""
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_DownSample')
		path,flag = '',''
		
		self.filename = obs.filename.split('.fits')[0]+'_lowres.fits'
		path = os.path.dirname(obs.filename)
		
		checkForFiles(self.logger,[self.filename],existence=True)
		
		# Downsample to 0.08deg
		#scale = 16
		
		nxn = obs.nx/scale # 64
		nyn = obs.ny/scale # 64
		dxn = scale*obs.dx # 0.08
		dyn = scale*obs.dy # 0.08
		
		data = zeros((obs.nz,nyn,nxn),dtype=float32)		
		
		Nax = obs.keyword['naxis']
		if Nax == 4:
			Tb = obs.observation[0,:,:,:]
		elif Nax == 3:
			Tb = obs.observation[:,:,:]
		else:
			self.logger.critical('NAXIS Array < 3 or > 4!')
			sys.exit(0)
		
		for i in xrange(nxn):
			xa = scale*i
			xb = xa+(scale-1)
			for j in xrange(nyn):
				ya = scale*j
				yb = ya+(scale-1)
				cosdec = cos(radians(obs.yarray[ya:yb+1]))
				data[:,j,i] = sum(sum(Tb[:,ya:yb+1,xa:xb+1]*cosdec,axis=1),axis=1)/scale**2
		
		# Invert longitude axis to match Galprop standard
		data[:,:,:] = data[:,:,::-1]
		
		# Build a skymap
		nxsky = fabs(round(360./dxn))
		nysky = fabs(round(180./dyn))
		
		skymap = zeros((obs.nz,nysky,nxsky),dtype=float32)
		sxarray = fabs(dxn)*arange(nxsky+1.)
		syarray = fabs(dyn)*arange(nysky+1.)
		
		x1 = amin(obs.xarray)
		y1 = amin(obs.yarray)+90.
		
		idx1 = where(sxarray<=x1)
		ix1 = amax(idx1[0])
		
		idy1 = where(syarray<=y1)
		iy1 = amax(idy1[0])
		
		# Broadcast CO data in the skymap
		skymap[:,iy1:iy1+nyn,ix1:ix1+nxn] = data
		
		obs.keyword['naxis1'] = nxsky #nxn
		obs.keyword['crpix1'] = 1. #round(nxn/2.)+1
		obs.keyword['cdelt1'] = -dxn
		obs.keyword['crval1'] = 0.25 #obs.x
		obs.keyword['naxis2'] = nysky #nyn
		obs.keyword['crpix2'] = 1. #round(nyn/2.+1)
		obs.keyword['cdelt2'] = dyn
		obs.keyword['crval2'] = -89.75 #obs.y

		obs.keyword['datamin'] = amin(skymap)
		obs.keyword['datamax'] = amax(skymap)
		
		obs.keyword['minfil'] = unravel_index(argmin(skymap),skymap.shape)[0]
		obs.keyword['mincol'] = unravel_index(argmin(skymap),skymap.shape)[1]
		obs.keyword['minrow'] = unravel_index(argmin(skymap),skymap.shape)[2]
				
		obs.keyword['maxfil'] = unravel_index(argmax(skymap),skymap.shape)[0]
		obs.keyword['maxcol'] = unravel_index(argmax(skymap),skymap.shape)[1]
		obs.keyword['maxrow'] = unravel_index(argmax(skymap),skymap.shape)[2]
		
		# Output file			
		results = pyfits.PrimaryHDU(skymap,obs.keyword)
		self.logger.info("Writing data to a fits file in...")
		results.writeto(self.filename, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
