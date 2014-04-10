#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits


class makeCGPSMosaic(object):
	
	def __init__(self,obs,mosaicConf):
		"""
		Generate mosaics of 'unabsorbed HI','HISA' (not HI+HISA becasue is the equivalent of the CGPS data),
		and velocity integrated 'CO' (Wco)
		"""
		self.species = obs.species
		self.logger = initLogger(mosaicConf['mosaic'], 'CGPS_GenerateMosaic')
		self.logger.info("Open file and get data...")
		
		# Get emission data and velocity interval
		Tb = obs.observation[:,:,:,:]
		dv = fabs(obs.msc_del_vel) # [velocity] = km s-1
		filter1 = where(Tb < 0.)
		Tb[filter1] = 0.

		# Usefull velocity channels
               	kmin = 18
               	kmax = 271
		
		if self.species == 'HI' or self.species == 'HISA':
		
			path = getPath(self.logger, key="cgps_hisa_dat")
			datafile = path+'CGPS_'+mosaicConf['mosaic']+'_HISA.dat'		
			checkForFiles(self.logger,[datafile])
			input = open(datafile,"r")
			lines = input.readlines()
						
			if self.species == 'HISA':
				result_array = zeros(Tb.shape,float)
			
			for line in lines:
				#print line.split('\n')[0]
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
			wco = zeros((obs.msc_lat,obs.msc_lon),dtype=float)
			wco = sum(Tb[0,kmin:kmax,:,:],axis=0)*dv
			
			# Write new header
			newheader = pyfits.Header()
			newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
			newheader.update(key="crval1", value=obs.msc_ref_lon, comment="Galactic longitude of reference pixel")
			newheader.update(key="crpix1", value=obs.msc_ind_lon, comment="Reference pixel in lon")
			newheader.update(key="cdelt1", value=obs.msc_del_lon, comment="Longitude increment")
			newheader.update(key="crota1", value=obs.msc_rot_lon, comment="Longitude rotation")
			newheader.update(key="cunit1", value="deg", comment="Unit type")
			newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
			newheader.update(key="crval2", value=obs.msc_ref_lat, comment="Galactic latitude of reference pixel")
			newheader.update(key="crpix2", value=obs.msc_ind_lat, comment="Reference pixel in lat")
			newheader.update(key="cdelt2", value=obs.msc_del_lat, comment="Latitude increment")
			newheader.update(key="crota2", value=obs.msc_rot_lat, comment="Latitude rotation")
			newheader.update(key="cunit2", value="deg", comment="Unit type")
			newheader.update(key="bunit", value="K km s-1", comment="Map units")
			newheader.update(key="datamin", value="%e"%amin(wco))
			newheader.update(key="datamax", value="%e"%amax(wco))
			newheader.update(key="object", value="CGPS Mosaic %s"%mosaicConf['mosaic'], comment="GCPS Mosaic")
			results = pyfits.PrimaryHDU(wco, newheader)
		
		# Output file
		self.logger.info("Write scaled data (int16) to a fits file in...")
		path = ''
		if self.species == 'HISA':
			path = getPath(self.logger, key="lustre_cgps_hisa")
			#if type == 'unabsorbed brightness temperature':
			results.writeto(path+'CGPS_'+mosaicConf['mosaic']+'_HISA_line.fits', output_verify='fix')
			#if type == 'amplitude':
			#	results.writeto(path+'CGPS_'+mosaicConf['mosaic']+'_HISA_amplitude.fits', output_verify='fix')
		elif self.species == 'HI':
			path = getPath(self.logger, key="lustre_cgps_hi")
			#if type == 'unabsorbed brightness temperature':
			results.writeto(path+'CGPS_'+mosaicConf['mosaic']+'_HI_unabsorbed_line.fits', output_verify='fix')
			#if type == 'amplitude':
			#	results.writeto(path+'CGPS_'+mosaicConf['mosaic']+'_HI_amplitude.fits', output_verify='fix')
		elif self.species == 'CO':
			path = getPath(self.logger, key="lustre_cgps_co")
			results.writeto(path+'CGPS_'+mosaicConf['mosaic']+'_WCO_line.fits', output_verify='fix')
		
		self.logger.info("%s"%path)
		self.logger.info("Done")

