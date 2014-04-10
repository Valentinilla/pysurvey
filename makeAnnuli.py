#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*


class makeAnnuli(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf):
		"""
		Calculate the column density of a mosaic (applying corrections) for the following species:
		'HI'(default), and 'HISA','HI+HISA','HI+CO','HI+HISA+CO' only for CGPS/SGPS
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.newspec #specified by the user
		self.type = mosaic.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_GenGalpropMap')
		file = ''
		path = ''
		flag = ''
		sur = (self.survey).lower()

		if self.species == 'HI':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			flag = 'HI_unabsorbed'
			if sur == 'lab':
				flag = 'HI'
		elif self.species == 'HISA':
			path = getPath(self.logger,'lustre_'+sur+'_hisa_column_density')
			flag = 'HISA'
		elif self.species == 'HI+HISA':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			flag = 'HI+HISA'
		elif self.species == 'CO':
			path = getPath(self.logger,'lustre_'+sur+'_co_column_density')
			flag = 'H2'
		elif self.species == 'WCO':
			path = getPath(self.logger,'lustre_'+sur+'_co_column_density')
			flag = 'H2'
		elif self.species == 'HI+CO':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			flag = 'HI+H2'
		elif self.species == 'HI+HISA+CO':
			path = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
			flag = 'HI+HISA+H2'

		file = path+self.survey+'_'+self.mosaic+'_galprop_'+flag+'_column_density.fits'
		checkForFiles(self.logger,[file],existence=True)
		
		self.logger.info("Open file and get data...")
				
		# Boundaries of Galactocentric Annuli - M.Ackermann et al
		# (The Astrophysical Journal, 750:3-38, 2012)
		annuli = 17
		rmin = [0.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.]
		rmax = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.,50.]
		ann_boundaries = []
		
		for a in range(0,annuli):
			ann_boundaries.append([a,rmin[a],rmax[a]])
		
		lon = mosaic.keyword['naxis1']
		lat = mosaic.keyword['naxis2']
		vel = mosaic.keyword['naxis3']

		# Array to store results
		N = zeros((annuli,lat,lon),dtype=float)
		
		Ts = float(utilsConf['tspin'])	  # [Excitation (or Spin) Temperature] = K (150)
		Tbg = float(utilsConf['tcmb'])	  # [Cosmic Microwave Background (CMB)] = K
		if not self.species == 'WCO':
			dv = fabs(mosaic.keyword['cdelt3']/1000.)	  # [velocity] = km s-1
		C = float(utilsConf['c'])	  # [costant] = cm-2
		
		# Corrections for latitude (almost negligible)
		cosdec = cos( radians(mosaic.yarray) ) # [lat. correction] = rad
		
		if self.survey == 'Galprop':
			if self.species == 'WCO':
				# Get emission data and velocity interval
				wco = sum(mosaic.observation,axis=0)
				# 2 accounts for the convertion from molecules to atoms
				N = 2*wco*float(utilsConf['xfactor'])*cosdec

		elif self.survey == 'LAB':
						
			# Get data		
			Tb = mosaic.observation[:,:,:]
			
			# Tb must be < Ts, otherwise problems with log
			index = where(Tb>=Ts)
			Tb[index] = 0.999*Ts
			
			# Optical depth correction
			# Without the continuum component
			cTb = log( (Ts)/(Ts-Tb) ) * Ts
			#cTb = -log(1-(Tb/(Ts-Tbg))) * Ts
			
			# Integrated brightness temperature over velocity (axis = 0)
			ITb = sum(cTb,axis=0)
			
			self.logger.info("Calculating NHI...")			
			# Column density
			N = C * ITb * dv	# [N] = cm-2
			# Corrected column density
			N = N*cosdec

		elif self.survey == 'Dame':
			if self.species == 'WCO':
				# Get emission data and velocity interval
				wco = mosaic.observation
				# 2 accounts for the convertion from molecules to atoms
				N = 2*wco*float(utilsConf['xfactor'])*cosdec
				
		elif self.survey == 'CGPS' or self.survey == 'SGPS':
			
			# Usefull velocity channels
			kmin = 0
			kmax = 0
			if self.survey == 'CGPS':
	               		kmin = 18
	               		kmax = 271
			if self.survey == 'SGPS':
	               		kmin = 1
        	       		kmax = 410
			
			# Used to skip calculation (see below) - do not change it!
			flagHI,flagCO,flagHISA = True, True, True
			
			# If the two column density files (or one of them) already exist, 
			# than add them up and skip the calculation (or skip the unnecessary one)		
			if self.species == 'HI+HISA':
				path1 = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
				file1 = path1+self.survey+'_'+self.mosaic+'_galprop_HI_unabsorbed_column_density.fits'
				
				path2 = getPath(self.logger,'lustre_'+sur+'_hisa_column_density')
				file2 = path2+self.survey+'_'+self.mosaic+'_galprop_HISA_column_density.fits'
				
				N,flagHI,flagHISA = checkToGetData(self.logger,file1,file2)
				
			if self.species == 'HI+CO':
				path1 = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
				file1 = path1+self.survey+'_'+self.mosaic+'_galprop_HI_unabsorbed_column_density.fits'
				
				path2 = getPath(self.logger,'lustre_'+sur+'_co_column_density')
				file2 = path2+self.survey+'_'+self.mosaic+'_galprop_H2_column_density.fits'
				
				N,flagHI,flagCO = checkToGetData(self.logger,file1,file2)
							
			if self.species == 'HI+HISA+CO':
				path1 = getPath(self.logger,'lustre_'+sur+'_hi_column_density')
				file1 = path1+self.survey+'_'+self.mosaic+'_galprop_HI_unabsorbed_column_density.fits'
				
				path2 = getPath(self.logger,'lustre_'+sur+'_hisa_column_density')
				file2 = path2+self.survey+'_'+self.mosaic+'_galprop_HISA_column_density.fits'
				
				path3 = getPath(self.logger,'lustre_'+sur+'_co_column_density')
				file3 = path3+self.survey+'_'+self.mosaic+'_galprop_H2_column_density.fits'

				N,flagHI,flagHISA,flagCO = checkToGetData(self.logger,file1,file2,file3)
						
			# Computing Column Density
			if(not self.species == 'CO'):
				# Get HI emission data
				Tb = mosaic.observation[0,:,:,:]
				# Setting the negative/0-values to Tcmb 
				#Tb = where( (Tb<0.) | (Tb==0.),Tbg,Tb)
				
				# HI Column Density						
				if(not self.species == 'HISA') and (flagHI == True):

					NHI = zeros((annuli,lat,lon),dtype=float)
					ITb = zeros((annuli,lat,lon),dtype=float)
					
					self.logger.info("Initializing parameters...")
					self.logger.info("1) Ts = %.2f K"%Ts)
					self.logger.info("2) dV = %.2f km/s"%dv)
					self.logger.info("3) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
					#self.logger.info("4) Tc(min) = %.2f K, Tc(max) = %.2f K"%(amin(Tc),amax(Tc)))
					
					self.logger.info("Calculating NHI...")
					# Optical depth correction
					# With the continuum component
					#Tfunc = (Ts-Tc)/(Ts-Tc-Tb)
					# Without the continuum component
					Tfunc = (Ts)/(Ts-Tb)
					
					Tfunc[Tfunc<1.] = 1. # <------ TO JUSTIFY
					cTb = log(Tfunc) * Ts
					
					# Integrated brightness temperature over each annulus
					for l in range(0,lon):
						for b in range(0,lat):
							for v in range(18,vel):
								glon = mosaic.xarray[l]
								glat = mosaic.yarray[b]
								vlsr = mosaic.zarray[v]/1000.
					 			#print glon,glat,vlsr
								d = rotCurveMPohl(self.logger,glon,glat,vlsr) #kpc
								i = find_ge(self.logger,rmax,d)
								annulus = ann_boundaries[i][0]
								ITb[annulus,b,l] += cTb[v,b,l]
								if not v%150:
									self.logger.info("(l,b) = (%i,%i) - lon=%.3f lat=%.3f d=%.1f ring=%i i=%i"%(l,b,glon,glat,d,annulus,i))
					
					# Column density
					NHI = C*ITb*dv # [NHI] = cm-2
					N = NHI*cosdec
					#exit(0)
				# HISA Column Density
				if(not self.species == 'HI') and (not self.species == 'HI+CO') and (flagHISA == True):
					NHISA = zeros((mosaic.observation.shape),dtype=float)
					
					# Get HI continuum data
					pathc = getPath(self.logger, sur+'_hi_continuum')
					continuum = pathc+self.survey+'_'+self.mosaic+'_1420_MHz_I_image.fits'
					checkForFiles(self.logger,[continuum])
					data, header = pyfits.getdata(continuum, 0, header = True)
					data[isnan(data)] = 0.
					Tc = data
					if self.survey == 'CGPS':
						Tc = data[0,0,:,:]
					if self.survey == 'SGPS':
						Tc = data[:,:]
						#Tc = where( (Tc<0.) | (Tc==0.),Tbg,Tc)

					# Get HISA data
					path_hisa_dat = getPath(self.logger,sur+'_hisa_dat')
					amplitude = path_hisa_dat+self.survey+'_'+self.mosaic+'_HISA.dat'
					checkForFiles(self.logger,[amplitude])
					input = open(amplitude,'r')
					lines = input.readlines()
					
					self.logger.info("Initializing parameters...")
					self.logger.info("1) dV = %.2f km/s"%dv)
					self.logger.info("2) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
					self.logger.info("3) Tc(min) = %.2f K, Tc(max) = %.2f K"%(amin(Tc),amax(Tc)))
					galax_dist = open("galaxy_"+self.mosaic+"_"+self.species+".dat","w")
					string = ""
					for line in lines:
						na = float(line.split('\n')[0].split()[0]) # HISA region
						nb = float(line.split('\n')[0].split()[1]) # HISA location
						Tu = float(line.split('\n')[0].split()[2]) # Unabsorbed brightness (Tu)
						dT = float(line.split('\n')[0].split()[3]) # Delta T (always negative)
						#print "%.2f\t%.2f\t%.2f\t%.2f"%(na,nb,nc,nd)
						
						m = floor(nb/1024/1024)
						ha = nb-m*1024*1024
						l = floor(ha/1024)
						k = ha-l*1024
						
						# Params
						glon = mosaic.xarray[k]
						glat = mosaic.yarray[l] #<--- 160? (it's in the martin's idl code)
						vlsr = mosaic.zarray[m]/1000.
												
						d = rotCurveMPohl(self.logger,glon,glat,vlsr) #kpc
						string = "%s %s %s %s\n"%(d,glon,glat,vlsr)
						galax_dist.write(string)
						
						theta = radians(mosaic.keyword['cdelt2']) #rad
						ds = d*tan(theta)*1e3 #pc
						
						A1 = float(utilsConf['pc2cm'])*float(utilsConf['poverk'])
						A2 = float(utilsConf['fn'])*ds/(C*dv)
						A = A1*A2
						
						B = Tc[l,k]+float(utilsConf['p'])*Tu
						
						init = [1.,10.]
						def equations(xx):
							'''
							Ts and tau functions - S.J.Gibson et al
							(The Astrophysical Journal, 540:852-862, 2000)
							'''
							tt, T = xx
							# Equation (6)
							Tfunc = (T-B)/(T-B-dT)
							if Tfunc<1.:
								Tfunc = 1.# Tbg # <------ TO JUSTIFY
							f1 = log(Tfunc)-tt
							# Equation (9)
							ttfunc = A/tt	
							if ttfunc<0.:
								ttfunc = 0.# <------ TO JUSTIFY
							f2 = sqrt(ttfunc)-T
														
							return array([f1, f2], dtype=float)
						
						(tau,Ts),infodict,ier,mesg = fsolve(equations,init, full_output=1)
						#plotFunc(tau,Ts)
						TsMin = Tbg
						if Ts < TsMin:
							Ts = TsMin
							Tfunc = (Ts-B)/(Ts-B-dT)
							tau = log(Tfunc)
						TsMax = Tc[l,k]+(dT+Tu)+(float(utilsConf['p'])-1.)*Tu
						if Ts > TsMax:
							# For Ts = TsMax, tau --> +oo
							Ts = TsMax
							tau = 1e4
						if tau < 0.:
							Ts = 0.
							tau = log(B/(B+dT))
							
						#print Ts,tau
						solution = False
						if ier == 1:
							solution = True
							NHISA[0,m,l,k] = tau*Ts*dv*C
							NHISA[0,m,l,k] = NHISA[0,m,l,k]*cos(radians(mosaic.yarray[l]))
						if not solution:
							#print "Could not find a valid solution:"
							#print mesg
							NHISA[0,m,l,k] = 0.
					
					galax_dist.close()
					
					self.logger.info("Calculating NHI...")
					# Corrected column density
					N = N + sum(NHISA[0,kmin:kmax,:,:],axis=0)

			if(self.species == 'CO' or self.species == 'HI+CO' or self.species == 'HI+HISA+CO') and (flagCO == True):
				# Get CO emission data
				pathco = getPath(self.logger, 'lustre_'+sur+'_co')
				co = pathco+self.survey+'_'+self.mosaic+'_CO_line.fits'
				checkForFiles(self.logger,[co])
				Wco, header = pyfits.getdata(co, 0, header = True)
				# 2 accounts for the convertion from molecules to atoms
				N = N + 2*Wco*float(utilsConf['xfactor'])*cosdec
				
			if self.species == 'CO':
	        	       	kmin = 30
	        	       	kmax = 230
				
				# Get emission data and velocity interval
				Tb = mosaic.observation[0,:,:,:]
				wco = zeros((mosaic.keyword['naxis2'],mosaic.keyword['naxis1']),dtype=float)
				wco = sum(Tb[kmin:kmax,:,:],axis=0)*dv
				
				# 2 accounts for the convertion from molecules to atoms
				N = 2*wco*float(utilsConf['xfactor'])*cosdec
		
		# Store results
		newheader = pyfits.Header()
		newheader["ctype1"] = ("GLON-CAR","Coordinate type")
		newheader["crval1"] = (mosaic.keyword["crval1"],"Galactic longitude of reference pixel")
		newheader["crpix1"] = (mosaic.keyword["crpix1"],"Reference pixel of lon")
		newheader["cdelt1"] = (mosaic.keyword["cdelt1"],"Longitude increment")
		newheader["crota1"] = (mosaic.keyword["crota1"],"Longitude rotation")
		newheader["cunit1"] = ("deg","Unit type")
		
		newheader["ctype2"] = ("GLAT-CAR","Coordinate type")
		newheader["crval2"] = (mosaic.keyword["crval2"],"Galactic latitude of reference pixel")
		newheader["crpix2"] = (mosaic.keyword["crpix2"],"Reference pixel of lat")
		newheader["cdelt2"] = (mosaic.keyword["cdelt2"],"Latitude increment")
		newheader["crota2"] = (mosaic.keyword["crota2"],"Latitude rotation")
		newheader["cunit2"] = ("deg","Unit type")

		newheader["ctype3"] = ("Ring","Coordinate type")
		newheader["crval3"] = (1,"Ring of reference pixel")
		newheader["crpix3"] = (1,"Reference pixel of ring")
		#newheader["cdelt3"] = (,"Ring increment")
		#newheader["crota3"] = (mosaic.keyword["crota3"],"Ring rotation")

		newheader["bunit"] = ("atoms cm-2","Map units")
		newheader["datamin"] = ("%e"%amin(N),"Min value")
		newheader["datamax"] = ("%e"%amax(N),"Max value")
		newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
		
		results = pyfits.PrimaryHDU(N, newheader)
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")


