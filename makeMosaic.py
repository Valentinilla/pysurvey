#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class makeMosaic(object):
	
	def __init__(self,obs,mosaicConf):
		"""
		Generate mosaics of 'unabsorbed HI','HISA' (not HI+HISA becasue is the equivalent of the CGPS data),
		and velocity integrated 'CO' (Wco)
		"""
		self.survey = obs.survey
		self.mosaic = obs.mosaic
		self.species = obs.species
		self.type = obs.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_GenerateMosaic')
		path,flag = '',''
		if self.survey == 'Galprop':
			if self.species == 'HI':
				path = getPath(self.logger, key="lustre_galprop_hi")
				self.mosaic = mosaicConf['mosaic']
				flag = 'HI'
			if self.species == 'WCO':
				path = getPath(self.logger, key="lustre_galprop_co")
				self.mosaic = mosaicConf['mosaic']
				flag = 'WCO'
		elif self.survey == 'LAB':
			path = getPath(self.logger, key="lustre_lab")
			if self.species == 'HI':
				self.mosaic = mosaicConf['mosaic']
				flag = 'HI'
		elif self.survey == 'Dame':
			path = getPath(self.logger, key="lustre_dame")
			if self.species == 'WCO':
				self.mosaic = mosaicConf['mosaic']
				flag = 'WCO'

		else:
			if self.species == 'HI':
				if self.survey == 'CGPS':
					path = getPath(self.logger, key="lustre_cgps_hi")
				if self.survey == 'SGPS':
					path = getPath(self.logger, key="lustre_sgps_hi")
				flag = 'HI_unabsorbed'
			elif self.species == 'HISA':
				if self.survey == 'CGPS':
					path = getPath(self.logger, key="lustre_cgps_hisa")
				if self.survey == 'SGPS':
					path = getPath(self.logger, key="lustre_sgps_hisa")
				flag = 'HISA'
			elif self.species == 'CO':
				path = getPath(self.logger, key="lustre_cgps_co")
				flag = 'CO'

		file = path+self.survey+'_'+self.mosaic+'_'+flag+'_line.fits'
		checkForFiles(self.logger,[file],existence=True)

		self.logger.info("Open file and get data...")

		if not self.survey == 'CGPS':

			self.mosaic = mosaicConf['mosaic']
			lon = mosaicConf['lon'] # deg
			lat = mosaicConf['lat'] # deg
			v1 = mosaicConf['vel1'] # km/s
			v2 = mosaicConf['vel2'] # km/s
			side = mosaicConf['side']
			
			if self.survey == 'Galprop':
				if v1 == 'INDEF' and v2 == 'INDEF':
					v1 = amin(obs.ring_array)
					v2 = amax(obs.ring_array)
				else:
					v1 = float(v1) # kpc
					v2 = float(v2) # kpc
				v1_px = int(ceil(obs.msc_ind_ring-1.+(v1-obs.msc_ref_ring)/obs.msc_del_ring))
				v2_px = int(ceil(obs.msc_ind_ring-1.+(v2-obs.msc_ref_ring)/obs.msc_del_ring))
				
			elif self.survey == 'LAB' or self.survey == 'SGPS':
				if v1 == 'INDEF' and v2 == 'INDEF':
					v1 = amin(obs.vel_array)
					v2 = amax(obs.vel_array)
				else:
					v1 = float(v1) # km/s
					v2 = float(v2) # km/s
				v1_px = int(ceil(obs.msc_ind_vel-1.+(v1-obs.msc_ref_vel)/obs.msc_del_vel))
				v2_px = int(ceil(obs.msc_ind_vel-1.+(v2-obs.msc_ref_vel)/obs.msc_del_vel))
			
			if lon == 'INDEF' and lat == 'INDEF':
				lon = obs.msc_ref_lon
				lat = obs.msc_ref_lat
			else:
				lon = float(lon) # deg
				lat = float(lat) # deg
			
			if side == 'INDEF':
				side = fabs(obs.msc_lon*obs.msc_del_lon)
			else:
				side = float(side)/2.  # deg
			
			l,b,sign = getMosaicCoordinate(self.logger,obs,self.survey,lon,lat,side)
			
			if not self.survey == 'Dame':
				v = [v1_px,v2_px]
				# Order indexes	
				if v1_px>v2_px:
					v = [v2_px,v1_px]
			
			crpix1 = int(round((l[1]-l[0])/2.))
			crpix2 = int(round((b[1]-b[0])/2.))
			
			self.logger.info("Mosaic properties:")
			self.logger.info("- (l0,b0) = (%.3f,%.3f) deg, [%s,%s] px"%(lon,lat,l[0]+crpix1,b[0]+crpix2))
			if self.survey == 'Galprop':
				self.logger.info("- (r1,r2) = (%.2f,%.2f) px"%(v[0],v[1]))
			elif self.survey == 'LAB' or self.survey == 'SGPS':
				self.logger.info("- (v1,v2) = (%.3f,%.3f) km/s, [%s,%s] px"%(v1,v2,v[0],v[1]))
			self.logger.info("- (h0,w0) = (%s,%s) deg, [%s,%s] px"%(side,side,crpix1,crpix2))
			
			if not self.survey == 'Dame':
				# The third axes can not be negative then I need to increase the index on the right v[1] by 1
				# in order to match the size of the origianl axes
				if v[0] == 0:
					v[0] = 1
					v[1] = v[1]+1
			
			newmosaic = ''
			if self.survey == 'Dame':
				newmosaic = obs.observation[b[0]:b[1],l[0]:l[1]]
			if self.survey == 'LAB':
				newmosaic = obs.observation[v[0]-1:v[1],b[0]:b[1],l[0]:l[1]] 
			if self.survey == 'SGPS':
				newmosaic = obs.observation[:,v[0]-1:v[1],b[0]:b[1],l[0]:l[1]]
			if self.survey == 'Galprop':
				newmosaic = obs.observation[v[0]-1:v[1],b[0]:b[1],l[0]:l[1]]
				# Revert the longitude
				newmosaic = newmosaic[:,:,::-1] #sintax [::-1] = [start:stop:step]
				obs.msc_del_lon=sign[0]*obs.msc_del_lon
				obs.msc_del_lat=sign[1]*obs.msc_del_lat
			if not self.survey == 'LAB':	
				obs.msc_rot_lon = 0
				obs.msc_rot_lat = 0
				obs.msc_rot_ring = 0
			
			# Write new header
			newheader = pyfits.Header()
			newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
			newheader.update(key="crval1", value=lon, comment="Galactic longitude of reference pixel")
			newheader.update(key="crpix1", value=crpix1, comment="Reference pixel of lon")
			newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
			newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
			newheader.update(key="cunit1", value="deg", comment="Unit type")
			newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
			newheader.update(key="crval2", value=lat, comment="Galactic latitude of reference pixel")
			newheader.update(key="crpix2", value=crpix2, comment="Reference pixel of lat")
			newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
			newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
			newheader.update(key="cunit2", value="deg", comment="Unit type")

			if self.survey == 'Dame':
				newheader.update(key="bunit", value="K km s-1", comment="Map units")
			
			if self.survey == 'Galprop':
				newheader.update(key="ctype3", value="RingBand", comment="Coordinate type")
				newheader.update(key="crval3", value=obs.ring_array[v[0]-1], comment="Ring of reference pixel")
				newheader.update(key="crpix3", value=1, comment="Reference pixel of ring")
				newheader.update(key="cdelt3", value=obs.msc_del_ring, comment="Ring increment")
				newheader.update(key="crota3", value=obs.msc_rot_ring, comment="Ring rotation")
				newheader.update(key="cunit3", value="ring num", comment="Unit type")
				newheader.update(key="bunit", value="K km s-1", comment="Map units")
			elif self.survey == 'LAB' or self.survey == 'SGPS':
				newheader.update(key="ctype3", value="VELO-LSR", comment="Coordinate type")
				newheader.update(key="crval3", value=obs.vel_array[v[0]-1]*1e3, comment="Velocity of reference pixel")
				newheader.update(key="crpix3", value=1, comment="Reference pixel of vel")
				newheader.update(key="cdelt3", value=obs.msc_del_vel*1e3, comment="Velocity increment")
				newheader.update(key="crota3", value=obs.msc_rot_vel, comment="Velocity rotation")
				newheader.update(key="cunit3", value="m/s", comment="Unit type")
				newheader.update(key="bunit", value="K", comment="Map units")
			
			newheader.update(key="system", value="GALACTIC", comment="Coordinate system")
			newheader.update(key="equinox", value=2000., comment="Equinox of ref. coord.")
			
			newheader.update(key="datamin", value=amin(newmosaic))
			newheader.update(key="datamax", value=amax(newmosaic))
			
			newheader.update(key="object", value=self.survey+" Mosaic "+self.mosaic, comment=self.survey+" Mosaic")
			if self.survey == 'LAB':
				newheader.update(key="freq0", value=obs.msc_band, comment="Rest frequency in Hz")
			if self.survey == 'SGPS':
				newheader.update(key="restfreq", value=obs.msc_band, comment="Rest frequency in Hz")
			
			print newheader
			results = pyfits.PrimaryHDU(newmosaic,newheader)
					
		else:
			
			# Get emission data and velocity interval
			Tb = obs.observation[:,:,:,:]
			dv = fabs(obs.msc_del_vel) # [velocity] = km s-1
			
			if self.species == 'HI' or self.species == 'HISA':
				
				#filter1 = where(Tb < 0.)
				#Tb[filter1] = 0.
					
				# Usefull velocity channels
        	       		kmin = 18
        	       		kmax = 271				
				
				path_data = getPath(self.logger, key="cgps_hisa_dat")
				datafile = path_data+self.survey+'_'+self.mosaic+'_HISA.dat'		
				checkForFiles(self.logger,[datafile])
				input = open(datafile,"r")
				lines = input.readlines()
							
				if self.species == 'HISA':
					result_array = zeros(Tb.shape,float)
				
				for line in lines:
					na = float(line.split('\n')[0].split()[0]) # HISA region
					nb = float(line.split('\n')[0].split()[1]) # HISA location
					nc = float(line.split('\n')[0].split()[2]) # Unabsorbed brightness (Tu) Tb only HI
					nd = float(line.split('\n')[0].split()[3]) # Delta T (always negative)  Tb only HISA
					#print "%.2f\t%.2f\t%.2f\t%.2f"%(na,nb,nc,nd)
					
					m = floor(nb/1024/1024)
					ha = nb-m*1024*1024
					l = floor(ha/1024)
					k = ha-l*1024
					
					if self.species == 'HISA':
						#if type == 'unabsorbed brightness temperature':
						#	result_array[0,m,l,k] = nc #Tb[0,m,l,k]-nd
						#elif type == 'amplitude':
						result_array[0,m,l,k] = abs(nd)
					elif self.species == 'HI':
						#if type == 'unabsorbed brightness temperature':
						Tb[0,m,l,k] = nc #Tb[0,m,l,k]-nd
						#if type == 'amplitude':
						#	Tb[0,m,l,k] = nd #Tb[0,m,l,k]-nd
				
			if self.species == 'HISA':
				obs.header['DATAMIN'] = amin(result_array)
				obs.header['DATAMAX'] = amax(result_array)
				results = pyfits.PrimaryHDU(result_array,obs.header)
				results.scale('int16', '', bscale=obs.msc_bscale, bzero=obs.msc_bzero)
			
			elif self.species == 'HI':
				obs.header['DATAMIN'] = amin(Tb)
				obs.header['DATAMAX'] = amax(Tb)
				results = pyfits.PrimaryHDU(Tb,obs.header)
				results.scale('int16', '', bscale=obs.msc_bscale, bzero=obs.msc_bzero)
			
			elif self.species == 'CO':
				
				self.logger.info("Applying Moment Mask method (T.M.Dame)...")
				T = Tb[0,:,:,:]
				
				# Calculate the rms for raw data
				rms_t = getRMS(self.logger,T)
				
				# Degrading the resolution spatially and in velocity by a factor of 2
				fwhm_spat = 2 #px
				sig = fwhm_spat/sqrt(8*log(2))
				Tsmooth = ndimage.gaussian_filter(T,sigma=(sig,sig,sig),order=(0,0,0))
				
				# Calculate the rms for smoothed data
				rms_ts = getRMS(self.logger,Tsmooth)
				
				# Set the clipping level equal 5 times the rms noise in Tsmooth
				Tclipping = 5*rms_ts
				
				# Generate a masking cube initially filled with zeros with the same dimensions as Tb
				Mask = zeros(Tsmooth.shape)
				
				# Unmask the pixel with a value > Tclipping
				index = where(Tsmooth>Tclipping)
				Mask[index] = 1
				
				# Calculate the moment-masked cube
				Tb[0,:,:,:] = Mask*T
				
				obs.header['DATAMIN'] = amin(Tb)
				obs.header['DATAMAX'] = amax(Tb)
				
				results = pyfits.PrimaryHDU(Tb,obs.header)
				results.scale('int16', '', bscale=obs.msc_bscale, bzero=obs.msc_bzero)		
		
		# Output file
		self.logger.info("Write scaled data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
