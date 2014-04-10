#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class spectralAnalysis(object):
	
	def __init__(self,obs,spectralConf):

		self.N_SPECTRAL        = int(spectralConf['n_spectral'])
		self.MAX_LOOPS         = int(spectralConf['max_loops'])
		self.RESIDUAL_FRAC     = float(spectralConf['residual_frac'])
		self.CLIP_SPECTRAL     = float(spectralConf['clip_spectral'])
		self.GAIN_SPECTRAL     = float(spectralConf['gain_spectral'])
		self.FWHM_SPECTRAL     = float(spectralConf['fwhm_spectral'])
		self.HISA_F_SPECTRAL   = float(spectralConf['hisa_f_spectral'])
		self.TEMP_CRITICAL     = float(spectralConf['temp_critical'])
		self.FIT_NARROW        = float(spectralConf['fit_narrow'])
		self.FIT_BROAD         = float(spectralConf['fit_broad'])
		self.FIT_QUAL          = float(spectralConf['fit_qual'])
		self.DIP_CUT           = float(spectralConf['dip_cut'])
		self.FWHM_SPATIAL_HISA = float(spectralConf['fwhm_spatial_hisa'])
		self.MIN_HISA          = float(spectralConf['min_hisa'])
		
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type
		
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_HISA_SpectralSearch')
		self.logger.info("Getting the data from...")
		self.logger.info("%s"%obs.filename)

		# Array to store results
		result_array=zeros((obs.nz,obs.ny,obs.nx),dtype=float32)

		# python convention
		# mosaic[s,v,y,x] i.e., indexing from right to left
		# s = stokes (I=1,Q=2,U=3,V=4), v/k = velocity, y/j = latitude, x/i = longitude
		
		# Set up gaussians to be used for convolution
		self.logger.info("Set up gaussians for convolution...")
		# spectral convilution
		sigma_spec = self.FWHM_SPECTRAL/sqrt(8*log(2))
		L = 3*self.FWHM_SPECTRAL
		dx = fabs(obs.dz/1000.)
		n_half = floor( (L/2)/dx )
		xi = linspace(-n_half*dx,n_half*dx,num=1+2*n_half)
		gauss_spec = gaussian(xi,[0,0,sigma_spec],normalized=True)
		gauss_spec *= 1/(gauss_spec).sum()
		#plotFunc(xi,gauss_spec)
		
		# spatial convolution (HISA gaussian)
		sigma_HISA = self.FWHM_SPATIAL_HISA/sqrt(8*log(2))
		n = 5
		rsqr = zeros((2*n+1,2*n+1),dtype=float32)
		for i in range(-n,n+1):
			for j in range(-n,n+1):
		    		rsqr[i+n,j+n] = i**2+j**2
		gauss_HISA = 1/(2*pi*sigma_HISA**2)*exp(-rsqr/(2*sigma_HISA**2))
		gauss_HISA *= 1/(gauss_HISA).sum()

		# Find sigma values for gaussian fits, in units of pixels
		sigma_narrow = self.FIT_NARROW/(sqrt(8*log(2))*abs(obs.dz/1000.))
		sigma_broad = self.FIT_BROAD/(sqrt(8*log(2))*abs(obs.dz/1000.))

		self.logger.info("Start of spectral search algorithm...")
		# needed for rms (see below)
		a=round((obs.nz-obs.zmin)/3)+17
		b=round((obs.nz-obs.zmin)/3)+18
		c=round(2*(obs.nz-obs.zmin)/3)+17
		d=round(2*(obs.nz-obs.zmin)/3)+18
		#print a,b,c,d #a=101, b=102, c=186, d=187

		printing = 0

		for i in xrange(7,obs.nx-9):#min:7 #max:1014
			for j in xrange(7,obs.ny-9):#7
				
				# O(k) is the observed spectrum
				observed = obs.observation[0,:,j,i]
				# S(k) is a smoothed version of O(k)
				smoothed = spatialAverage1D(obs.observation[0,:,:,:],i,j,self.N_SPECTRAL)
				# U(k) is the unabsorbed spectrum (initially set to zero)
				unabsorbed = zeros(observed.shape,dtype=float32)
				# R(k) is the residual spectrum (initially set equal to S(k))
				residual = smoothed
				# Spectral rms noise in O(k)
				noise = observed - smoothed
				sigma_array = zeros(3,dtype=float32)
				# Compute the standard deviation along the specified axis:
				# std = sqrt(mean(abs(x - x.mean())**2))
				sigma_array[0] = std(noise[18:a])
				sigma_array[1] = std(noise[b:c])
				sigma_array[2] = std(noise[d:])
				sigma_obs = amin(sigma_array)
				
				# Initialize the HISA container
				HISA_merged = zeros(obs.msc_vel,dtype=float)
				
				# Initialize the rms of positive values of residual
				sigma_pos = 0
				
				# Clean loop
				smax = amax(smoothed)
				
				for loop in range(0,self.MAX_LOOPS-1):

					rmax = amax(residual)
				
					# 1. 
					if(rmax < self.RESIDUAL_FRAC*smax or rmax == 0.):
						break
					
					# 2.
					correction = zeros(obs.msc_vel,dtype=float)
					correct = where(residual > self.CLIP_SPECTRAL*rmax) 
					correction[correct] = self.GAIN_SPECTRAL*residual[correct]

					# 3.
					# two convolution algorithms (the second is faster)
					#unabsorbed += convolve(correction,gauss_spec,"same")
					unabsorbed += fftconvolve(correction,gauss_spec,"same")

					# 4.
					# calculate the rms of positive values of residual
					residual = smoothed - unabsorbed
					x_pos = residual[residual > 0.]
				
					if(len(x_pos) > 0):
						sigma_pos = sqrt(mean( power(x_pos,2) ))
					else:
						sigma_pos = 0.
					
					if(sigma_pos < sigma_obs):
						break
					
				# Process results to look for potential HISA
				if(rmax > 0):
					suspected_hisa_index = where(residual < sigma_pos*self.HISA_F_SPECTRAL)
					# the counting starts from right to left of the graph (remember the 18 initial 0-events)
					cnt1 = size(suspected_hisa_index)
				else: 
					cnt1 = 0

				# Skip HISA search if HICA present
				if(amin(observed) < -6.*sigma_obs):
					cnt1 = 0				

				counter = 0
				segments=[]
				# Loop over the elements of suspected_hisa_index
				if(printing):
					self.logger.info("Suspected HISA index_array: %s"%suspected_hisa_index)
					self.logger.info("Suspected HISA value_array: %s"%residual[residual < sigma_pos*self.HISA_F_SPECTRAL])
					self.logger.info("Number of suspected HISA: %i"%cnt1)

				# Build consecutive HISA candidates into segments
				while(counter < cnt1):
					segments = [residual[suspected_hisa_index[0][counter]]]
					kmin = [suspected_hisa_index[0][counter],0]
					count_start = counter
					counter += 1
					while(counter < cnt1 and (suspected_hisa_index[0][counter]-1) == suspected_hisa_index[0][counter-1]):
						segments.append(residual[suspected_hisa_index[0][counter]])
						if(segments[counter-count_start] < segments[kmin[1]] ):
							kmin = [suspected_hisa_index[0][counter],counter-count_start]
						counter += 1
					# kmin[0][0] = index in residual, kmin[0][1] = index in segments, kmin[1] = size of segment
					kmin = [kmin, counter-count_start]
	
					# Fit Gaussians to candidate segments
					#perform initial filtering on candidate segments
					if(smoothed[kmin[0][0]] > 2.*self.HISA_F_SPECTRAL*sigma_pos and unabsorbed[kmin[0][0]] > self.TEMP_CRITICAL):

						# Pad array (segments) with 0s to make sure fit possible
						num_of_zeros = 6 #6
						data = zeros(kmin[1]+num_of_zeros)
						i_inf = 0.5*num_of_zeros
						i_sup = kmin[1] + 0.5*num_of_zeros
						data[i_inf:i_sup] = abs(array(segments))
				
						amp = amax(data)
						mu = kmin[0][0]
						params = [amp,mu,sigma_narrow,sigma_broad]

						x = xrange( kmin[1] + num_of_zeros ) + ( kmin[0][0] - (kmin[0][1]+ 0.5*num_of_zeros) )

						if(printing):
							self.logger.info("[%s] Fit:"%counter)
							self.logger.info("- (amp, mu, sigma)_narrow = (%s, %s, %s)"%(params[0],params[1],params[2]))
							self.logger.info("- (amp, mu, sigma)_broad  = (%s, %s, %s)\n"%(params[0],params[1],params[3]))

						paramsN = Parameters()
						paramsB = Parameters()

						if((data.size-num_of_zeros) > 3): #3 (at least 3 points to fit)
							paramsN.add('amp', value=params[0], vary=True)
							paramsN.add('mu', value=params[1], vary=False)
							paramsN.add('sigma', value=params[2], vary=False)
							paramsB.add('amp', value=params[0], vary=True)
							paramsB.add('mu', value=params[1], vary=False)
							paramsB.add('sigma', value=params[3], vary=False)
						else:
							paramsN.add('amp', value=params[0], vary=False)
							paramsN.add('mu', value=params[1], vary=True)
							paramsN.add('sigma', value=params[2], vary=True)
							paramsB.add('amp', value=params[0], vary=False)
							paramsB.add('mu', value=params[1], vary=True)
							paramsB.add('sigma', value=params[3], vary=True)

						fitN = minimize(residualG, paramsN, args = (x, data), engine = "leastsq")
						fitB = minimize(residualG, paramsB, args = (x, data), engine = "leastsq")

						if(printing):
							printResults(self.logger,fitN,paramsN,"Narrow")
							printResults(self.logger,fitB,paramsB,"Broad")

						bestfit_paramsN = [float(paramsN['amp'].value),float(paramsN['mu'].value),float(paramsN['sigma'].value)]
						bestfit_paramsB = [float(paramsB['amp'].value),float(paramsB['mu'].value),float(paramsB['sigma'].value)]

						#p = (amp,mu,sigma_narrow)
						#gfit = gaussian(x,p,normalized=False)
						#gfitN = gaussian(x,bestfit_paramsN,normalized=False)
						#gfitB = gaussian(x,bestfit_paramsB,normalized=False)
						#plotFunc(x,data,gfitN,gfitB)
						
						# Check quality of fit
						sigma_fitN = sqrt((fitN.residual**2).sum()/data.size)
						sigma_fitB = sqrt((fitB.residual**2).sum()/data.size)

						condition = min([bestfit_paramsN[0]/max([sigma_fitN,sigma_pos]),bestfit_paramsB[0]/max([sigma_fitB,sigma_pos])])

						
						if( (condition > self.FIT_QUAL) and fitN.success and fitB.success):
							fit_quality = True
						else:
							fit_quality = False

						if(printing):
							self.logger.info("Fit_quality %s"%fit_quality)
						
						# Process and combine fitted Gaussians
						if(fit_quality):
							HISA_narrow = gaussian(arange(obs.msc_vel),bestfit_paramsN,normalized=False)
							HISA_broad  = gaussian(arange(obs.msc_vel),bestfit_paramsB,normalized=False)

							# Only take values > 5% of peak (cut the gaussian tails)
							filterN = 0.05 * HISA_narrow[kmin[0][0]]
							filterB = 0.05 * HISA_broad[kmin[0][0]]

							HISA_narrow[HISA_narrow < filterN] = 0
							HISA_broad[HISA_broad < filterB] = 0

							HISA_merged_temp = where(HISA_narrow > HISA_broad, HISA_narrow, HISA_broad)

							# Check for "dip" in O(k) at k_min
							dip_frac = dipFilter(observed,kmin[0][0],num_of_zeros)
							if(dip_frac > self.DIP_CUT):
								HISA_merged += HISA_merged_temp

							#steps = xrange(obs.msc_vel)
							#ar = -1*ones(observed.shape,dtype=float64)
							#trh1 = ones(obs.msc_vel)*(sigma_pos*self.HISA_F_SPECTRAL-10)
							#trh2 = 30*ones(obs.msc_vel)
							#f1 = [observed,unabsorbed,smoothed,10*ar+residual]
							#f2 = [60*ar+HISA_narrow,90*ar+HISA_broad,120*ar+HISA_merged]
							# U(k),S(k),R(k),HISA_N,HISA_B,HISA_M
							#plotFunc(steps,f1[1],f1[0],f1[3],f2[0],f2[1],f2[2],trh1,trh2)
							#self.logger.info("Done with (i,j) = (%s,%s) and counter %s"%(i+1,j+1,counter))

				# Store HISA in result_array
				result_array[:,j,i] = HISA_merged 
					
			self.logger.info("Done with (i,j) = (%s,%s)"%(i,j))
		
		# Do spatial smoothing of HISA
		self.logger.info("Spatial smoothing (convolution)...")
		for k in range(0,obs.msc_vel-1):
			result_array[k,:,:] = fftconvolve(result_array[k,:,:],gauss_HISA,"same")
			if not k%10:
				self.logger.info("- done with k = %s"%k)

		self.logger.info("Only take smoothed HISA > 2K...")
		# Turn results into a map with 1 where smoothed HISA > 2K and 0 everywhere else
		result_array[result_array < self.MIN_HISA] = 0.

		obs.observation[0,:,:,:] = result_array

		#results = result_array.astype(uint8)
		results = pyfits.PrimaryHDU(obs.observation,obs.header)
		results.scale('int16', '',bscale=obs.msc_bscale, bzero=obs.msc_bzero)

		# Output file
		self.logger.info("Write scaled data (int16) to a fits file...")
		file = '/lustre/fs4/group/that/sf/Survey/CGPS/HISA/CGPS_'+commonConf['mosaic']+'_spectral_search_pass'+commonConf['spectral_iteration']+'.fits'
		results.writeto(file, output_verify='fix')
		self.logger.info("Done")


