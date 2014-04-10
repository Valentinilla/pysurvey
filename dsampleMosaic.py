#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class dsampleMosaic(object):
	
	def __init__(self,obs,scale):
		"""
		Downsample mosaics: scale x dx_msc = new_dx_msc
		"""
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_DownSample')
		path,flag = '',''

		self.filename = obs.filename.split('.fits')[0]+'_lowres.fits'#_x%i.fits'%scale
		path = os.path.dirname(obs.filename)
		
		checkForFiles(self.logger,[self.filename],existence=True)

		#zarray = z+dz*(arange(nz)+1.-px)
		#zarray = (-60.-obs.dz*(arange(obs.nz)-144.))

		# Downsample to 0.08deg
		#scale = 16

		nxn = obs.nx/scale # 64
		nyn = obs.ny/scale # 64
		dxn = scale*obs.dx # 0.08
		dyn = scale*obs.dy # 0.08

		pxn = 0.5*(obs.xarray[7]+obs.xarray[8])
		pyn = 0.5*(obs.yarray[7]+obs.yarray[8])

		data = zeros((1,obs.nz,nyn,nxn),dtype=float32)

		for i in xrange(0,nxn):
			xia=scale*i
			xib=xia+(scale-1)
			for j in xrange(0,nyn):
				yia=scale*j
				yib=yia+(scale-1)
				#print i,j
				#print yia,yib,xia,xib
				data[0,:,j,i] = sum(sum(obs.observation[0,:,yia:yib+1,xia:xib+1],axis=1),axis=1)/scale**2
				#data[0,:,j,i] = mean(sum(observation[0,:,yia:yib,xia:xib]))
		
		# write out
		#keyword = pyfits.Header()
		#sxdelpar,header,['BZERO','BSCALE']
		obs.keyword['naxis1'] = nxn
		obs.keyword['crpix1'] = 1
		obs.keyword['cdelt1'] = dxn
		obs.keyword['crval1'] = pxn
		obs.keyword['naxis2'] = nyn
		obs.keyword['crpix2'] = 1
		obs.keyword['cdelt2'] = dyn
		obs.keyword['crval2'] = pyn

		obs.keyword['datamin'] = amin(data)
		obs.keyword['datamax'] = amax(data)
	
		obs.keyword['minfil'] = unravel_index(argmin(data),data.shape)[1]
		obs.keyword['mincol'] = unravel_index(argmin(data),data.shape)[2]
		obs.keyword['minrow'] = unravel_index(argmin(data),data.shape)[3]
				
		obs.keyword['maxfil'] = unravel_index(argmax(data),data.shape)[1]
		obs.keyword['maxcol'] = unravel_index(argmax(data),data.shape)[2]
		obs.keyword['maxrow'] = unravel_index(argmax(data),data.shape)[3]
		
		#path = ''
		#file = path+'CGPS_'+self.mosaic+'_'+self.species+'_reduced.fits'
		
		# Output file			
		results = pyfits.PrimaryHDU(data,obs.keyword)
		self.logger.info("Writing data to a fits file in...")
		results.writeto(self.filename, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
