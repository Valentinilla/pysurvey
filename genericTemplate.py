#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class genericMosaic(object):
	
	def __init__(self,obs,mosaicConf):
		"""
		Generate mosaics of 'unabsorbed HI','HISA' (not HI+HISA becasue is the equivalent of the CGPS data),
		and velocity integrated 'CO' (Wco)
		"""
		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_GenerateMosaic')
		
		f = pyfits.open(filename)
		hdu = f[0]
		self.msc_hdu = hdu
		
		if self.survey == 'Galprop':
						
			self.msc_size = hdu.header['NAXIS']	# 4
			self.msc_lon = hdu.header['NAXIS1']	# 1024
			self.msc_lat = hdu.header['NAXIS2']	# 1024
			self.msc_rband = hdu.header['NAXIS3']	# 1024
			self.msc_ref_lon = float(hdu.header['CRVAL1']) # deg (GLON-CAR)
			self.msc_del_lon = float(hdu.header['CDELT1']) # -4.9e-3 deg
			self.msc_ind_lon = float(hdu.header['CRPIX1']) # 513 px
			self.msc_ref_lat = float(hdu.header['CRVAL2']) # deg
			self.msc_del_lat = float(hdu.header['CDELT2']) # 4.9e-3 deg
			self.msc_ind_lat = float(hdu.header['CRPIX2']) # 513 px
			self.msc_ref_rband = float(hdu.header['CRVAL3']) # deg
			self.msc_del_rband = float(hdu.header['CDELT3']) # 4.9e-3 deg
			self.msc_ind_rband = float(hdu.header['CRPIX3']) # 513 px
			
			self.msc_bunit = float(hdu.header['BUNIT']) # 'K km s-1' or '10e20 cm-2'
			
			self.lon_array=self.msc_ref_lon+self.msc_del_lon*(arange(self.msc_lon)+1.-self.msc_ind_lon)
			self.lat_array=self.msc_ref_lat+self.msc_del_lat*(arange(self.msc_lat)+1.-self.msc_ind_lat)
			self.rband_array=self.msc_ref_rband+self.msc_del_rband*(arange(self.msc_rband)+1.-self.msc_ind_rband)
			self.observation = hdu.data



		file = path+self.survey+'_'+self.mosaic+'_'+flag+'_line.fits'
		checkForFiles(self.logger,[file],existence=True)

		self.logger.info("Open file and get data...")

		if self.survey == 'Dame':
			
			self.mosaic = mosaicConf['mosaic']
			lon = mosaicConf['lon'] # deg
			lat = mosaicConf['lat'] # deg
			side = mosaicConf['side']
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
				
			# Check if the coordinates in .cfg are inside the mosaic object
			obj_lon_max = obs.msc_ref_lon+((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
			obj_lon_min = obs.msc_ref_lon-((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
			obj_lat_max = obs.msc_ref_lat+((obs.msc_lat-1.)*obs.msc_del_lat)/2.
			obj_lat_min = obs.msc_ref_lat-((obs.msc_lat-1.)*obs.msc_del_lat)/2.
			
			lon_max = lon+side
			lon_min = lon-side
			lat_max = lat+side
			lat_min = lat-side
			
			if lon_max>obj_lon_max or lon_min<obj_lon_min:
				self.logger.critical("The longitude within the .cfg file")
				self.logger.critical("doesn't match that of the mosaic!")
				sys.exit(0)
			if lat_max>obj_lat_max or lat_min<obj_lat_min:
				self.logger.critical("The latitude within the .cfg file")
				self.logger.critical("doesn't match that of the mosaic!")	
				sys.exit(0)
			
			l1,l2,b1,b2 = getMosaicCoordinate(obs,lon,lat,side)
			crpix1 = int(ceil((l2-l1)/2.))
			crpix2 = int(ceil((b2-b1)/2.))
			
			newmosaic = obs.observation[b1:b2,l1:l2]
			obs.msc_rot_lon = 0
			obs.msc_rot_lat = 0
			obs.msc_rot_vel = 0
			
			# Write new header
			newheader = pyfits.Header()
			newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
			newheader.update(key="crval1", value=lon, comment="Galactic longitude of reference pixel")
			newheader.update(key="crpix1", value=crpix1, comment="Reference pixel in lon")
			newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
			newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
			newheader.update(key="cunit1", value="deg", comment="Unit type")
			newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
			newheader.update(key="crval2", value=lat, comment="Galactic latitude of reference pixel")
			newheader.update(key="crpix2", value=crpix2, comment="Reference pixel in lat")
			newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
			newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
			newheader.update(key="cunit2", value="deg", comment="Unit type")
			newheader.update(key="bunit", value="K km s-1", comment="Map units")
			newheader.update(key="datamin", value="%e"%amin(newmosaic))
			newheader.update(key="datamax", value="%e"%amax(newmosaic))
			newheader.update(key="object", value="Mosaic "+self.mosaic, comment=self.survey+" Mosaic")
			results = pyfits.PrimaryHDU(newmosaic, newheader)	

		if self.survey == 'LAB' or self.survey == 'SGPS':

			self.mosaic = mosaicConf['mosaic']
			lon = mosaicConf['lon'] # deg
			lat = mosaicConf['lat'] # deg
			v1 = mosaicConf['vel1'] # km/s
			v2 = mosaicConf['vel2'] # km/s
			side = mosaicConf['side']
			
			#mosaic = MW2
			#lon = 140.75
			#lat = 3.0
			#vel1= 286
			#vel2= 503
			#side = 5.12

			if v1 == 'INDEF' and v2 == 'INDEF':
				v1 = amin(obs.vel_array)
				v2 = amax(obs.vel_array)
			else:
				v1 = float(v1) # km/s
				v2 = float(v2) # km/s
			
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

			v1_px = int(round(obs.msc_ind_vel-1.+(v1-obs.msc_ref_vel)/obs.msc_del_vel))
			v2_px = int(round(obs.msc_ind_vel-1.+(v2-obs.msc_ref_vel)/obs.msc_del_vel))
			#print v1_px,v2_px
			#print obs.vel_array[v1],obs.vel_array[v2]
			#print obs.vel_array.shape
			side = float(mosaicConf['side'])/2. # deg
	
			# Check if the coordinates in .cfg are inside the mosaic object
			obj_lon_max = obs.msc_ref_lon+((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
			obj_lon_min = obs.msc_ref_lon-((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
			if self.survey == 'LAB':
				obj_lat_max = obs.msc_ref_lat+((obs.msc_lat-1.)*obs.msc_del_lat)
				obj_lat_min = obs.msc_ref_lat
			else:
				obj_lat_max = obs.msc_ref_lat+((obs.msc_lat-1.)*obs.msc_del_lat)/2.
				obj_lat_min = obs.msc_ref_lat-((obs.msc_lat-1.)*obs.msc_del_lat)/2.
			
			lon_max = lon+side
			lon_min = lon-side
			lat_max = lat+side
			lat_min = lat-side
			
			if lon_max>obj_lon_max or lon_min<obj_lon_min:
				self.logger.critical("The longitude within the .cfg file")
				self.logger.critical("doesn't match that of the mosaic!")
				sys.exit(0)
			if lat_max>obj_lat_max or lat_min<obj_lat_min:
				self.logger.critical("The latitude within the .cfg file")
				self.logger.critical("doesn't match that of the mosaic!")	
				sys.exit(0)

			#print lon_min,lon_max,lat_min,lat_max
			#print obj_lon_min,obj_lon_max,obj_lat_min,obj_lat_max
			#sys.exit(0)
			l1,l2,b1,b2 = getMosaicCoordinate(obs,lon,lat,side)
		
			crpix1 = int(ceil((l2-l1)/2.))
			crpix2 = int(ceil((b2-b1)/2.))

			self.logger.info("Mosaic properties:")
			self.logger.info("- (l0,b0) = (%.3f,%.3f) deg"%(lon,lat))
			self.logger.info("- (x0,y0) = (%s,%s) px"%(l1+crpix1,b1+crpix2))
			self.logger.info("- (v1,v2) = (%.3f,%.3f) km/s"%(v1,v2))
			self.logger.info("- (h0,w0) = (%s,%s) deg"%(side,side))
			self.logger.info("- (w0,h0) = (%s,%s) px"%(crpix1,crpix2))

			#print l1,l2, b1,b2,v1_px,v2_px
			#print obs.vel_array[min(v1_px,v2_px)-1]*1e3
			#sys.exit(0)
			
			newmosaic = ''
			if v1_px == 0:
				v1_px = 1
			if self.survey == 'LAB':
				newmosaic = obs.observation[v1_px-1:v2_px,b1:b2,l1:l2]

			if self.survey == 'SGPS':
				newmosaic = obs.observation[:,v1_px-1:v2_px,b1:b2,l1:l2]
				obs.msc_rot_lon = 0
				obs.msc_rot_lat = 0
				obs.msc_rot_vel = 0
			#print newmosaic.shape
			#print obs.vel_array[min(v1_px,v2_px)+1]*1e3
				
			newheader = pyfits.Header()
			
			newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
			newheader.update(key="crval1", value=lon, comment="Galactic longitude of reference pixel")
			newheader.update(key="crpix1", value=crpix1, comment="Reference pixel in lon")
			newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
			newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
			newheader.update(key="cunit1", value="deg", comment="Unit type")
	
			newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
			newheader.update(key="crval2", value=lat, comment="Galactic latitude of reference pixel")
			newheader.update(key="crpix2", value=crpix2, comment="Reference pixel in lat")
			newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
			newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
			newheader.update(key="cunit2", value="deg", comment="Unit type")
	
			newheader.update(key="ctype3", value="VELO-LSR", comment="Coordinate type")
			newheader.update(key="crval3", value=obs.vel_array[min(v1_px,v2_px)-1]*1e3, comment="Velocity of reference pixel")
			#newheader.update(key="crval3", value=obs.msc_ref_vel*1e3, comment="Velocity of reference pixel")
			#newheader.update(key="crpix3", value=obs.msc_ind_vel, comment="Reference pixel in vel")
			newheader.update(key="crpix3", value=1, comment="Reference pixel in vel")
			newheader.update(key="cdelt3", value=obs.msc_del_vel*1e3, comment="Velocity increment")
			newheader.update(key="crota3", value=obs.msc_rot_vel, comment="Velocity rotation")
			newheader.update(key="cunit3", value="m/s", comment="Unit type")
			
			newheader.update(key="system", value="GALACTIC", comment="Coordinate system")
			newheader.update(key="bunit", value="K", comment="Map units")
			newheader.update(key="equinox", value=2000., comment="Equinox of ref. coord.")
		
			newheader.update(key="datamin", value=amin(newmosaic))
			newheader.update(key="datamax", value=amax(newmosaic))
	
			newheader.update(key="object", value=self.survey+" Mosaic "+self.mosaic, comment=self.survey+" Mosaic")
			if self.survey == 'LAB':
				newheader.update(key="freq0", value=obs.msc_band, comment="Rest frequency in Hz")
			if self.survey == 'SGPS':
				newheader.update(key="restfreq", value=obs.msc_band, comment="Rest frequency in Hz")
			
			results = pyfits.PrimaryHDU(newmosaic,newheader)
			
		elif self.survey == 'CGPS':

			# Get emission data and velocity interval
			Tb = obs.observation[:,:,:,:]
			dv = fabs(obs.msc_del_vel) # [velocity] = km s-1
			
			if self.species == 'HI' or self.species == 'HISA':
				
				filter1 = where(Tb < 0.)
				Tb[filter1] = 0.
					
				# Usefull velocity channels
        	       		kmin = 18
        	       		kmax = 271				
				
				path = getPath(self.logger, key="cgps_hisa_dat")
				datafile = path+self.survey+'_'+self.mosaic+'_HISA.dat'		
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
				
				#wco = zeros((obs.msc_lat,obs.msc_lon),dtype=float)
				#wco = sum(Tb[0,kmin:kmax,:,:],axis=0)*dv
				# Write new header
				#newheader = pyfits.Header()
				#newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
				#newheader.update(key="crval1", value=obs.msc_ref_lon, comment="Galactic longitude of reference pixel")
				#newheader.update(key="crpix1", value=obs.msc_ind_lon, comment="Reference pixel in lon")
				#newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
				#newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
				#newheader.update(key="cunit1", value="deg", comment="Unit type")
				#newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
				#newheader.update(key="crval2", value=obs.msc_ref_lat, comment="Galactic latitude of reference pixel")
				#newheader.update(key="crpix2", value=obs.msc_ind_lat, comment="Reference pixel in lat")
				#newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
				#newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
				#newheader.update(key="cunit2", value="deg", comment="Unit type")
				#newheader.update(key="bunit", value="K km s-1", comment="Map units")
				#newheader.update(key="datamin", value="%e"%amin(wco))
				#newheader.update(key="datamax", value="%e"%amax(wco))
				#newheader.update(key="object", value="Mosaic "+self.mosaic, comment=self.survey+" Mosaic")
				#results = pyfits.PrimaryHDU(wco, newheader)			
		
		# Output file
		self.logger.info("Write scaled data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
		
