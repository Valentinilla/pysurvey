#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*

class splitMosaic(object):
	
	def __init__(self,mosaic,ntot):
		"""
		Calculate the column density in galacto-centric rings using the rotation curve of M. Phol et al.
		Boundaries of galacto-centric annuli by M.Ackermann et al.
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.species
		self.type = mosaic.type
		self.ntot = ntot

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_SplitModule')
		file,path,flag = '','',''
		sur = (self.survey).lower()
		
		if self.species == 'HI':
			path = getPath(self.logger,'lustre_'+sur+'_hi_split')
			flag = 'HI_unabsorbed_line'
		elif self.species == 'HISA':
			path = getPath(self.logger,'lustre_'+sur+'_hisa_split')
			flag = 'HISA_line'
		elif self.species == 'CO':
			self.logger.critical("Split module not implemented for CO.")
			sys.exit(0)
				
		self.logger.info("Open file and get data...")
				
		# Get HI emission data
		Tb = mosaic.observation[:,:,:,:]
		lat = mosaic.yarray
		
		alist = array_split(Tb, self.ntot, axis=2)
		ind_lat = 0
		
		for z in xrange(0,self.ntot):

			file = "%s%s_%s_%s_part_%i-%i.fits"%(path,self.survey,self.mosaic,flag,z+1,self.ntot)
			checkForFiles(self.logger,[file],existence=True)
			
			# Store results
			crpix_lat = (Tb.shape[2]/self.ntot)/2
			i = z+(z+1)
			ind_lat = crpix_lat*i
			#print i, ind_lat	
			mosaic.keyword['crval2'] = lat[ind_lat-1]
			mosaic.keyword['crpix2'] = crpix_lat
			
			#crpix_vel = (Tb.shape[1]/self.ntot)/2
			#i = z+(z+1)
			#ind_vel = crpix_vel*i			
			#mosaic.keyword["crval3"] = vel[ind_vel-1]
			#mosaic.keyword["crpix3"] = crpix_vel

			mosaic.keyword['datamin'] = (amin(alist[z]),"Min value")
			mosaic.keyword['datamax'] = (amax(alist[z]),"Max value")
			
			mosaic.keyword['minfil'] = unravel_index(argmin(alist[z]),alist[z].shape)[0]
			mosaic.keyword['mincol'] = unravel_index(argmin(alist[z]),alist[z].shape)[1]
			mosaic.keyword['minrow'] = unravel_index(argmin(alist[z]),alist[z].shape)[2]
			mosaic.keyword['maxfil'] = unravel_index(argmax(alist[z]),alist[z].shape)[0]
			mosaic.keyword['maxcol'] = unravel_index(argmax(alist[z]),alist[z].shape)[1]
			mosaic.keyword['maxrow'] = unravel_index(argmax(alist[z]),alist[z].shape)[2]
			
			mosaic.keyword['object'] = ("Mosaic %s (%i/%i)"%(self.mosaic,z+1,self.ntot),"%s Mosaic (part/tot)"%self.survey)
			
			# Output file
			results = pyfits.PrimaryHDU(alist[z], mosaic.keyword)
			self.logger.info("Writing fits file %i out of %s in..."%(z+1,self.ntot))
			results.writeto(file,output_verify='fix')
			self.logger.info("%s"%path)
			self.logger.info("Done")
