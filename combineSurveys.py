#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class combineSurveys(object):
	
	def __init__(self,surveylist,mosaiclist,species,res):
		"""
		Allow to combine surveys of column density in 2D or 3D map.
		species = HI, HI_unabsorbed+HISA
		"""
		self.species = species
		self.logger = initLogger(self.species+'_CombineSurveys')

		if not (self.species == 'HI' or self.species == 'HI_unabsorbed+HISA'):
			self.logger.critical("Wrong species! Only 'HI' and 'HI_unabsorbed+HISA' allowed.")
			sys.exit(0)

		nsurveys = len(surveylist)		
		self.logger.info("Number of survey: %s"%nsurveys)
		if nsurveys == 0:
			self.logger.critical("Survey list empty!")
			sys.exit(0)
		elif nsurveys == 1:
			self.logger.info("Summing HI_unabsorbed and HISA...")
			# you want to sum HI_unabs and HISA
			survey = surveylist[0]
			mosaic = mosaiclist[0]
			path1 = getPath(self.logger, key='lustre_'+survey.lower()+'_hi_unabsorbed')+'skymaps/'
			path2 = getPath(self.logger, key='lustre_'+survey.lower()+'_hisa')+'skymaps/'
		
			file1 = path1+'skymap_hi_unabsorbed_'+survey.lower()+'_rbands_r9_res'+str(res)+'_'+mosaic+'.fits'
			checkForFiles(self.logger,[file1])
			f1 = pyfits.open(file1)
			keyword1 = f1[0].header
			skymap1 = f1[0].data
		
			file2 = path2+'skymap_hisa_'+survey.lower()+'_rbands_r9_res'+str(res)+'_'+mosaic+'.fits'
			checkForFiles(self.logger,[file2])
			f2 = pyfits.open(file2)
			keyword2 = f2[0].header
			skymap2 = f2[0].data
		
			skymap = zeros(skymap1.shape,dtype=float32)
			skymap = where(skymap1==skymap2,skymap1,skymap1+skymap2)

		else:
			self.logger.info("Combining %s surveys..."%surveylist)
			
			key = [[] for i in range(nsurveys)]
			sky = [[] for i in range(nsurveys)]
			for i,srv in enumerate(surveylist):
				path = getPath(self.logger, key='lustre_'+srv.lower()+'_hi')+'skymaps/'
				file = path+'skymap_'+self.species.lower()+'_'+srv.lower()+'_rbands_r9_res'+str(res)+'_'+mosaiclist[i]+'.fits'
				f = pyfits.open(file)
				key[i] = f[0].header
				sky[i] = f[0].data	
			
			path = '/afs/ifh.de/group/that/work-sf/survey/batch/survey_skymaps/hi/' 
			file = path+'skymap_hi_galprop_rbands_r9_res'+str(res)+'.fits'
			f = pyfits.open(file)
			okey = f[0].header
			osky = f[0].data

			skymap = osky
			i = where(osky<>sky[0])
			j = where(osky<>sky[1])
			skymap[i] = sky[0][i]
			skymap[j] = sky[1][j]
			#skymap = zeros(sky[0].shape,dtype=float32)
			#for i,srv in enumerate(surveylist-1):
			#	skymap = where(skymap==sky[i+1],0,1)
			keyword1 = okey
			keyword1['object'] = ("Skymap","Surveys: %s"%surveylist)
		# Store results
		#newheader = pyfits.Header()
		#keyword1["ctype1"] = ("GLON-CAR","Coordinate type")
		#newheader["crval1"] = (crval1,"Galactic longitude of reference pixel")
		#newheader["crpix1"] = (crpix1,"Reference pixel of lon")
		#newheader["cdelt1"] = (msc_dx,"Longitude increment")
		#newheader["crota1"] = (msc_rotx,"Longitude rotation")
		#newheader["cunit1"] = ("deg","Unit type")
		
		#keyword1["ctype2"] = ("GLAT-CAR","Coordinate type")
		#newheader["crval2"] = (crval2,"Galactic latitude of reference pixel")
		#newheader["crpix2"] = (crpix2,"Reference pixel of lat")
		#newheader["cdelt2"] = (msc_dy,"Latitude increment")
		#newheader["crota2"] = (msc_roty,"Latitude rotation")
		#newheader["cunit2"] = ("deg","Unit type")

		#newheader['ctype3'] = ("Rband","Coordinate type")
		#newheader['crval3'] = (0,"Ring of reference pixel")
		#newheader['crpix3'] = (1.0,"Reference pixel of ring")
		#newheader['cdelt3'] = (1,"Ring increment")
		#newheader['crota3'] = (msc_rotz,"Ring rotation")

		#newheader['bunit'] = (msc_bunit,"Map units")
		keyword1['datamin'] = (amin(skymap),"Min value")
		keyword1['datamax'] = (amax(skymap),"Max value")
		#if self.mosaic == 'skymap':
		#	newheader['object'] = (self.survey+" Skymap",self.survey+" Mosaic")
		#else:
		#	newheader['object'] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		keyword1['minfil'] = unravel_index(argmin(skymap),skymap.shape)[0]
		keyword1['mincol'] = unravel_index(argmin(skymap),skymap.shape)[1]
		keyword1['minrow'] = unravel_index(argmin(skymap),skymap.shape)[2]
		keyword1['maxfil'] = unravel_index(argmax(skymap),skymap.shape)[0]
		keyword1['maxcol'] = unravel_index(argmax(skymap),skymap.shape)[1]
		keyword1['maxrow'] = unravel_index(argmax(skymap),skymap.shape)[2]

		#if 'history' in hdu1.header:
		#	for i in xrange(len(hdu1.header['history'])):
		#		newheader['history'] = hdu1.header['history'][i]
		
		results = pyfits.PrimaryHDU(skymap,keyword1)
				
		# Output file
		self.logger.info("Writing data to a fits file in...")
		if nsurveys == 1:
			path = getPath(self.logger, key='lustre_'+survey.lower()+'_hi')+'skymaps/'
			skymap_name = path+'skymap_hi_unabsorbed+hisa_'+survey.lower()+'_rbands_r9_res'+str(res)+'_'+mosaic+'.fits'
		else:
			srv_name = '_'.join(surveylist).lower()
			path = '/afs/ifh.de/group/that/work-sf/survey/batch/survey_skymaps/hi/'
			skymap_name = path+'skymap_'+self.species.lower()+'_'+srv_name+'_rbands_r9_res'+str(res)+'.fits'
		
		rmin,rmax,annuli = getAnnuli(glob_annuli)
		
		# Create a Table with the annuli boundaries
		col1 = pyfits.Column(name='Rmin', format='1E', unit='kpc', array=array(rmin))
		col2 = pyfits.Column(name='Rmax', format='1E', unit='kpc', array=array(rmax))
		cols = pyfits.ColDefs([col1,col2])
		tbl = pyfits.new_table(cols)
		tbl.name = "BINS"
			
		thdulist = pyfits.HDUList([results,tbl])			
		thdulist.writeto(skymap_name, output_verify='fix')

		self.logger.info("%s"%path)
		self.logger.info("Done")

		
