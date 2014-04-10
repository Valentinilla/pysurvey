#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class GalpropSkymap(object):
	
	def __init__(self,obs,res):
		"""
		Downsample mosaics: res = dx_msc x scale
		Input can be either 'mosaic' or 'skymap'.
		"""
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_GalpropSkymap')
		path,flag = '',''
		
		self.filename = obs.filename.split('.fits')[0]+'_lowres.fits'
		path = os.path.dirname(obs.filename)
		
		checkForFiles(self.logger,[self.filename],existence=True)
		
		# Downsample to 0.08deg
		scale = res/fabs(obs.dx)
		#print scale
		
		nxn = round(obs.nx/scale) # 64
		nyn = round(obs.ny/scale) # 64
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
		
		#for i in xrange(nxn):
		#	xa = scale*i
		#	xb = xa+(scale-1)
		#	print xa,xb
		#	exit(0)
		#	for j in xrange(nyn):
		#		ya = scale*j
		#		yb = ya+(scale-1)
		#		cosdec = cos(radians(obs.yarray[ya:yb+1]))
		#		data[:,j,i] = sum(sum(Tb[:,ya:yb+1,xa:xb+1]*cosdec,axis=1),axis=1)/scale**2

		def rebin(a, newshape):
			'''
			Rebin an array to a new shape.
			'''
			assert len(a.shape) == len(newshape)
			slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
			coordinates = mgrid[slices]
			indices = coordinates.astype('i')   #choose the biggest smaller integer index
			return a[tuple(indices)]

		for k in xrange(obs.nz):
			data[k,:,:] = rebin(Tb[k,:,:],data[k,:,:].shape)

		# Invert longitude axis to match Galprop standard
		data[:,:,:] = data[:,:,::-1]
		
		# Build a skymap
		nxsky = fabs(round(360./dxn))
		nysky = fabs(round(180./dyn))
		#print nxsky,nysky,dxn,dyn,round(obs.nx/scale)

		skymap = zeros((obs.nz,nysky,nxsky),dtype=float32)
		sxarray = fabs(dxn)*arange(nxsky+1.)
		syarray = fabs(dyn)*arange(nysky+1.)

		x1 = amin(obs.xarray)
		x2 = amax(obs.xarray)
		y1 = amin(obs.yarray)+90.
		y2 = amax(obs.yarray)+90.

		idx1 = where(sxarray<=x1)
		ix1 = amax(idx1[0])
		
		idx2 = where(sxarray<=x2)
		ix2 = amax(idx2[0])
		
		idy1 = where(syarray<=y1)
		iy1 = amax(idy1[0])
		
		idy2 = where(syarray<=y2)
		iy2 = amax(idy2[0])
		
		def indexes(x1,x2,i1,i2,arr,nx):
			
			eps1 = [fabs(arr[i1-1]-x1),i1-1]
			eps2 = [fabs(arr[i1]-x1),i1]
			eps3 = [fabs(arr[i1+1]-x1),i1+1]
			eps4 = [fabs(arr[i2-1]-x2),i2-1]
			eps5 = [fabs(arr[i2]-x2),i2]
			eps6 = [fabs(arr[i2+1]-x2),i2+1]
			eps = min(eps1[0],eps2[0],eps3[0],eps4[0],eps5[0],eps6[0])
			idx1,idx2 = 0,0
			
			if eps == eps1[0]:
				idx1 = eps1[1]
				idx2 = idx1+nx
			elif eps == eps2[0]:
				idx1 = eps2[1]
				idx2 = idx1+nx
			elif eps == eps3[0]:
				idx1 = eps3[1]
				idx2 = idx1+nx
			elif eps == eps4[0]:
		                idx2 = eps4[1]
		                idx1 = idx2-nx
		        elif eps == eps5[0]:
		                idx2 = eps5[1]
		                idx1 = idx2-nx
		        elif eps == eps6[0]:
		                idx2 = eps6[1]
		                idx1 = idx2-nx
		
			return int(idx1),int(idx2)
		
		i1,i2 = indexes(x1,x2,ix1,ix2,sxarray,nxn)
		j1,j2 = indexes(y1,y2,iy1,iy2,syarray,nyn)
				
		data[data==0.]=-999
		# Broadcast data in skymap excluding bad pixels
		if self.species == 'CO':
			skymap[:,j1+5:j2-1,i1+2:i2-15] = data[:,5:nyn-1,2:nxn-15]
		else:	
			if self.survey == 'CGPS':
				#skymap[:,j1+1:j2-1,i1+11:i2] = data[:,1:nyn-1,11:]
				skymap[:,j1+1:j2-1,i1+11+25:i2] = data[:,1:nyn-1,11+25:]
			elif self.survey == 'SGPS':
				skymap[:,j1:j2,i1+4:i2-1] = data[:,:,4:nxn-1]
			else:
				skymap[:,j1:j2,i1:i2] = data
		
		# Open original Galprop map
		if self.species == 'HI' or self.species == 'HI_unabsorbed' or self.species == 'HISA': flag = 'hi'
		if self.species == 'CO': flag = 'co'
		path2 = '/afs/ifh.de/group/that/work-sf/survey/batch/survey_skymaps/'+flag
		file2 = path2+'/skymap_'+flag+'_galprop_rbands_r9_res'+str(res)+'.fits'
		checkForFiles(self.logger,[file2])
		f2 = pyfits.open(file2)
		keyword2 = f2[0].header
		observation2 = f2[0].data

		cubemap = zeros(observation2.shape,dtype=float32)
		cubemap = where(skymap==0.,observation2,skymap)

		skymap[skymap==-999] = 0.
		cubemap[cubemap==-999] = 0.
				
		# Update keywords
		obs.keyword['naxis1'] = nxsky #nxn
		obs.keyword['crpix1'] = 1. #round(nxn/2.)+1
		obs.keyword['cdelt1'] = -dxn
		obs.keyword['crval1'] = fabs(dxn/2)#0.25 #obs.x
		obs.keyword['naxis2'] = nysky #nyn
		obs.keyword['crpix2'] = 1. #round(nyn/2.+1)
		obs.keyword['cdelt2'] = dyn
		obs.keyword['crval2'] = -90+fabs(dyn/2)#-89.75 #obs.y

		obs.keyword['datamin'] = amin(skymap)
		obs.keyword['datamax'] = amax(skymap)
		
		obs.keyword['minfil'] = unravel_index(argmin(skymap),skymap.shape)[0]
		obs.keyword['mincol'] = unravel_index(argmin(skymap),skymap.shape)[1]
		obs.keyword['minrow'] = unravel_index(argmin(skymap),skymap.shape)[2]
				
		obs.keyword['maxfil'] = unravel_index(argmax(skymap),skymap.shape)[0]
		obs.keyword['maxcol'] = unravel_index(argmax(skymap),skymap.shape)[1]
		obs.keyword['maxrow'] = unravel_index(argmax(skymap),skymap.shape)[2]

		if 'history' in obs.keyword:
			for i in xrange(len(obs.keyword['history'])):
				keyword2['history'] = obs.keyword['history'][i]
		
		# Output file			

		# Create a Table with the annuli boundaries
		rmin,rmax,annuli = getAnnuli(glob_annuli)
		col1 = pyfits.Column(name='Rmin', format='1E', unit='kpc', array=array(rmin))
		col2 = pyfits.Column(name='Rmax', format='1E', unit='kpc', array=array(rmax))
		cols = pyfits.ColDefs([col1,col2])
		tbl = pyfits.new_table(cols)
		tbl.name = "BINS"

		# Writing low resolution map
		results = pyfits.PrimaryHDU(skymap,obs.keyword)
		self.logger.info("Writing low resolution map in...")
				
		thdulist = pyfits.HDUList([results,tbl])
		thdulist.writeto(self.filename, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")

		# Writing final Galprop map
		temp = self.filename.split('_')[-2]
		self.logger.info("Writing final Galprop map in...")
		results2 = pyfits.PrimaryHDU(cubemap,keyword2)
		galfile = path+'/skymap_'+self.species.lower()+'_'+self.survey.lower()+'_rbands_r9_res'+str(res)+'_'+temp+'.fits'
		
		thdulist2 = pyfits.HDUList([results2,tbl])
		thdulist2.writeto(galfile, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
