#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class spatialAnalysis(object):
	
	def __init__(self,obs,commonConf,spatialConf):

		self.N_SPATIAL         = int(spatialConf['n_spatial'])
		self.MAX_LOOPS         = int(spatialConf['max_loops'])
		self.HIGH              = int(spatialConf['high'])
		self.RESIDUAL_FRAC     = float(spatialConf['residual_frac'])
		self.CLIP_SPATIAL      = float(spatialConf['clip_spatial'])
		self.GAIN_SPATIAL      = float(spatialConf['gain_spatial'])
		self.FWHM_SPATIAL      = float(spatialConf['fwhm_spatial'])
		self.NOISE_RESOLVE     = float(spatialConf['noise_resolve'])
		self.HISA_F_SPATIAL    = float(spatialConf['hisa_f_spatial'])
		self.TEMP_CRITICAL     = float(spatialConf['temp_critical'])
		self.AMP_MIN_FIRST     = float(spatialConf['amp_min_first'])
		self.FWHM_SPATIAL_HISA = float(spatialConf['fwhm_spatial_hisa'])
		self.MIN_HISA          = float(spatialConf['min_hisa'])

		self.logger = initLogger(obs.msc_area, 'HISA_Spatial_Search')
		Print(self.logger,spatialConf, 'spatial search')
		self.logger.info("Open file and get data...")

		# Array to store results
		result_array=zeros((obs.msc_vel,obs.msc_lat,obs.msc_lon),dtype=float)

		# python convention
		# mosaic[s,v,y,x] i.e., indexing from right to left
		# s = stokes (I=1,Q=2,U=3,V=4), v/k = velocity, y/j = latitude, x/i = longitude
		
		# Set up gaussians to be used for convolution
		self.logger.info("Set up gaussians for convolution...")
		# spatial gaussian
		# (convert to distance in pixels [67px])
		fwhm_spat = (self.FWHM_SPATIAL/abs(obs.msc_del_lon)/60.)
		sigma_spat = fwhm_spat/sqrt(8*log(2))
		npix = 2*fix(1.5*fwhm_spat) # 200px
		x1 = linspace(-npix,npix,num=2*npix+1) # center 0
		r1,r2 = meshgrid(x1,x1)
		gauss_spat = gauss_kern(npix) #gaussian_2d(r1,r2,0.,0.,sigma_spat,sigma_spat)
		#plotFunc2D(r1,r2,gauss_spat)

		#x = arange(npix)
		#x01,x02 = meshgrid(x,x)
		#gauss_spat = gaussian_spat_2d(npix,100.,100.,sigma_spat,sigma_spat)
		#plotFunc2D(x01,x02,gauss_spat)
		#print x
		#print fwhm_spat,sigma_spat,npix,amax(gauss_spat),sum(gauss_spat)
		#exit(0)

		# HISA gaussian
		sigma_spat_HISA = self.FWHM_SPATIAL_HISA/sqrt(8*log(2)) # FWHM = 5px
		npix2 = 2*ceil(1.5*self.FWHM_SPATIAL_HISA) # 16px
		x1 = linspace(-npix2,npix2,num=npix2+1) # center 0
		r1,r2 = meshgrid(x1,x1)
		gauss_spat_HISA = gaussian_2d(r1,r2,0.,0.,sigma_spat_HISA,sigma_spat_HISA)
		#plotFunc2D(r1,r2,gauss_spat_HISA)

		#xi = xrange(obs.msc_lon)
		#plotFunc(xi,residual[700,:],unabsorbed[700,:])

		self.logger.info("Start of spatial search algorithm...")
		kmin = 121#18
		kmax = 122#271

		# Start of spatial search algorithm
		for k in range(kmin,kmax): #18,270

			self.logger.info("k = %s"%k)
			
			# O(i,j) is the observed spectrum
			observed = obs.observation[0,k,:,:]
			# S(i,j) is a smoothed version of O(i,j)
			#smoothed = blur_image(observed,self.N_SPATIAL) #15px
			smoothed = spatialAverage2D(observed,self.N_SPATIAL) #15px
			# U(i,j) is the unabsorbed spectrum (initially set to zero)
			unabsorbed = zeros(observed.shape,dtype=float)
			# R(k) is the residual spectrum (initially set equal to S(k))
			residual = smoothed
			
			#ar = random.rand(9,9) #arange(81).reshape((9,9))
			#br = old_spatialAverage2D(ar,3)
			#br = plu_spatialAverage2D(ar,3)
			#xi = xrange(obs.msc_lon)
			#plotFunc(xi,observed[700,:])
			
			# N(i,j) is the smoothed noise on a large angular scale [20' = 66.7 px]
			noise_box = self.NOISE_RESOLVE/(60.*abs(obs.msc_del_lon)) #67px
			#noise = blur_image(observed,noise_box)
			noise = spatialAverage2D(observed,noise_box)
			
			# Estimatation of the rms noise in O(i,j)
			sigma_obs = rms_estimation2D(observed,noise,noise_box)
			
			# Estimation of the average rms noise in S(i,j)	
			sigma_sm = zeros(observed.shape,dtype=float)
			sigma_sm1 = rms_estimation2D(smoothed,noise,noise_box)
			
			imin = int(noise_box/2.)
			imax_lat = int(obs.msc_lat-(imin+1))
			imax_lon = int(obs.msc_lon-(imin+1))
			
			for y in range(imin,imax_lat):
				for x in range(imin,imax_lon):
					min1 = min(sigma_sm1[y-imin,x-imin],sigma_sm1[y-imin,x+imin])
					min2 = min(sigma_sm1[y+imin,x-imin],sigma_sm1[y+imin,x+imin])
					sigma_sm[y,x] = min(min1,min2)
			
			den = (obs.msc_lon-noise_box)**2
			sigma_sm_avg = (sigma_sm[imin:obs.msc_lat-(imin+1),imin:obs.msc_lon-(imin+1)]).sum()/den
			
			# Clean loop
			s_max = amax(smoothed)
			self.logger.info("- Start of CLEAN loop...")
			
			for loop in range(0,self.MAX_LOOPS-1):

				# Find the 10th highest residual
				rmax = get_nth_maxvalue(residual,self.HIGH)
				
				# 1. 
				if(rmax < self.RESIDUAL_FRAC*s_max or rmax == 0):
					self.logger.info("- loop %s exit (1st condition)"%loop)
					self.logger.info("- rmax = %s"%rmax)
					break
				
				# 2.
				correction = zeros(observed.shape,dtype=float)
				index = where(residual > self.CLIP_SPATIAL*rmax) 
				correction[index] = self.GAIN_SPATIAL*residual[index]
				
				# 3.
				unabsorbed += fftconvolve(correction,gauss_spat,"same")

				# 4.
				# calculate the rms of positive values of residual
				residual = smoothed - unabsorbed
				
				x_pos = residual[residual > 0.]
				if( len(x_pos) > 0):
					sigma_pos = sqrt(mean( power(x_pos,2) ))
				else:
					sigma_pos = 0.
				
				if(sigma_pos < sigma_sm_avg):
					self.logger.info("- loop %s exit (2nd condition)"%loop)
					self.logger.info("- (sig_pos, sig_sm) = (%.2f K, %.2f K)"%(sigma_pos,sigma_sm_avg))
					break
				
				if not loop%10 or loop == self.MAX_LOOPS-1:
					self.logger.info("- (sig_pos, sig_sm) = (%.2f K, %.2f K)"%(sigma_pos,sigma_sm_avg))
					self.logger.info("- rmax = %s"%rmax)
					self.logger.info("- loop %s done"%loop)

				if 0:
					print "correction sum=%.3f, max=%.3f"%(correction.sum(),amax(correction))
					print "unabsorbed sum=%.3f, max=%.3f"%(unabsorbed.sum(),amax(unabsorbed))
					print "residual   sum=%.3f, max=%.3f"%(residual.sum(),amax(residual))
					print "rmax=%.3f, sig_pos=%.3f, sig_sm=%.3f\n"%(rmax,sigma_pos,sigma_sm_avg)
					
				if loop == -2:
					exit(0)

			# Process results to look for potential HISA
			observed_dif = observed-unabsorbed

			# eliminate pixels with no valid noise calculations from consideration
			omit = [arange(noise_box),obs.msc_lon-arange(noise_box)-1]
			omit2 = [arange(noise_box),obs.msc_lat-arange(noise_box)-1]
			residual[:,omit] = 0.
			residual[omit2,:] = 0.
			observed_dif[:,omit] = 0.
			observed_dif[omit2,:] = 0.

			res_check = where(residual < self.HISA_F_SPATIAL*sigma_sm)
			obs_check = where(observed_dif < self.HISA_F_SPATIAL*sigma_obs)
			
			hisa_detected = zeros(observed.shape,dtype=float)
			if(len(res_check) > 0): 
				hisa_detected[res_check] = unabsorbed[res_check]-observed[res_check]
			if(len(obs_check) > 0):
				hisa_detected[obs_check] = unabsorbed[obs_check]-observed[obs_check]
			
			# perform initial filtering of detected HISA
			filter1 = where(unabsorbed < self.TEMP_CRITICAL)
			filter2 = where(hisa_detected < self.AMP_MIN_FIRST)
			if(len(filter1) > 0):
				hisa_detected[filter1] = 0.
			if(len(filter2) > 0):
				hisa_detected[filter2] = 0.
			
			# Do spatial smoothing of HISA
			#result_array[k,:,:] = hisa_detected
	 		result_array[k,:,:] = fftconvolve(hisa_detected,gauss_spat_HISA,"same")

		# Do spatial smoothing of HISA
		#self.logger.info("Spatial smoothing (convolution)...")
		#for k in range(0,obs.msc_vel-1):
		#	result_array[k,:,:] = fftconvolve(result_array[k,:,:],gauss_spat_HISA,"same")
		#	if not k%10:
		#		self.logger.info("- done with k = %s"%k)
			
		self.logger.info("Only take smoothed HISA > 2K...")
		# Turn results into a map with 1 where smoothed HISA > 2K and 0 everywhere else
		result_array[result_array < self.MIN_HISA] = 0.
	
		obs.observation[0,:,:,:] = result_array
		
		#results = result_array.astype(uint8)
		results = pyfits.PrimaryHDU(obs.observation,obs.header)
		results.scale('int16', '',bscale=obs.msc_bscale, bzero=obs.msc_bzero)

		# Output file
		self.logger.info("Write scaled data (int16) to a fits file...")
		results.writeto('/lustre/fs4/group/that/sf/HISA/CGPS_'+commonConf['mosaic']+'_spatial_search.fits', output_verify='fix')
		self.logger.info("Done")


