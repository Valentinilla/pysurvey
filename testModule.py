#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*

class testModule(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf):
		"""
		Calculate the column density of a mosaic (applying corrections) for the following species:
		'HI'(default), and 'HISA','HI+HISA','HI+CO','HI+HISA+CO' only for CGPS/SGPS
		"""
		self.survey = mosaic.survey
		self.mosaic = mosaic.mosaic
		self.species = mosaic.newspec #specified by the user
		self.type = mosaic.type

		self.logger = initLogger(self.survey+'_'+self.mosaic+'_'+self.species+'_TestModule')
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

		file = path+self.survey+'_'+self.mosaic+'_test_'+flag+'_ring_column_density.fits'
		checkForFiles(self.logger,[file],existence=True)
		
		self.logger.info("Open file and get data...")
				
		# Boundaries of Galactocentric Annuli - M.Ackermann et al
		# (The Astrophysical Journal, 750:3-38, 2012)
		annuli = 17
		rmin = [0.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.]
		rmax = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.,50.]
		ann_boundaries = []
		
		for a in xrange(0,annuli):
			ann_boundaries.append([a,rmin[a],rmax[a]])
		
		nlon = mosaic.nx
		nlat = mosaic.ny
		nvel = mosaic.nz

		lon = mosaic.xarray
		lat = mosaic.yarray
		vel = mosaic.zarray/1000.
		
		# Array to store results
		cubemap = zeros((annuli,nlat,nlon),dtype=float)
		
		Ts = float(utilsConf['tspin'])	  # [Excitation (or Spin) Temperature] = K (150)
		Tbg = float(utilsConf['tcmb'])	  # [Cosmic Microwave Background (CMB)] = K
		if not self.species == 'WCO':
			dv = fabs(mosaic.dz/1000.)	  # [velocity] = km s-1
		C = float(utilsConf['c'])	  # [costant] = cm-2
		
		# Corrections for latitude (almost negligible)
		cosdec = cos(radians(lat)) # [lat. correction] = rad
		
		# Usefull velocity channels
	        kmin = 18
	    	kmax = 271
			
		# Computing Column Density
		# Get HI emission data
		Tb = mosaic.observation[0,:,:,:]
		# free memory
		del mosaic.observation

		# Setting the negative/0-values to Tcmb 
		#Tb = where( (Tb<0.) | (Tb==0.),Tbg,Tb)
				
		# HI Column Density						
		#NHI = zeros((annuli,nlat,nlon),dtype=float)
		#ITb = zeros((annuli,nlat,nlon),dtype=float)
					
		self.logger.info("Initializing parameters...")
		self.logger.info("1) Ts = %.2f K"%Ts)
		self.logger.info("2) dV = %.2f km/s"%dv)
		self.logger.info("3) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
		#self.logger.info("4) Tc(min) = %.2f K, Tc(max) = %.2f K"%(amin(Tc),amax(Tc)))
					
		self.logger.info("Calculating gas distribution...")
					
		#Limits: galactic longitude < |165 deg|, galactic latitude < |5 deg|.
		path2 = getPath(self.logger,'rotcurve_mpohl')	
	
		# Read in SPH model results
		# Get Bissantz's data
		file2 = path2+'testvr1.dat'
		input = open(file2,'r')
		bissantz = input.read().split()
		n1 = 200
		vrs = array(bissantz).reshape(n1,n1)
		vrs = vrs.astype(float)
		# free memory
		del bissantz
		#print vrs[100,98]  #-79.1474
		#print vrs[100,99]  #-56.3561
		#print vrs[100,100] #-25.6225
		
		R = 0.5
		xa = -10.+0.1*(R+arange(n1))
		ya = -10.+0.1*(R+arange(n1))
		rb = zeros((n1,n1),dtype=float)
		for i in xrange(0,n1-1):
			rb[:,i] = sqrt(xa[i]**2+ya**2)
		ia = where(rb > 8.0)
		vrs[ia] = 0.
	
		# Position of sun in SPH model and 
		# unit vectors (e) to GC (tangential) 
		# and l = 90 (normal) direction
		x_sun_sph = 7.518  #kpc
		y_sun_sph = -2.735 #kpc
		r_sun_sph = sqrt(x_sun_sph**2+y_sun_sph**2)
		ex_tan = -x_sun_sph/r_sun_sph
		ey_tan = -y_sun_sph/r_sun_sph
		ex_norm = -ey_tan
		ey_norm = ex_tan
		xha = zeros(3800,dtype=int)
		yha = zeros(3800,dtype=int)

		# Line properties
		sigma_co = 3.       #co velocity dispersion [km s-1]
		sigma_co_inner = 5. #co velocity dispersion inner galaxy [km s-1]
		
		ivzero = int(floor(20./dv))  #24 CO CGPS; 17 LAB
		iv_vec = arange(2*ivzero+1)  #49 CO CGPS

		# Define dummy line profile and its weight for W_co
		line = exp(-0.5*(dv*(iv_vec-ivzero)/sigma_co)**2)
		sigma_line = sigma_co*sqrt(2.*pi)
		lip = exp(-0.5*(dv*(iv_vec-ivzero)/(sigma_co/2))**2)
		wlip = sigma_line/2.
		lim = line/sigma_line

		# Line profile for GC region
		line_inner = exp(-0.5*(dv*(iv_vec-ivzero)/sigma_co_inner)**2)
		sigma_line_inner = sigma_co_inner*sqrt(2.*pi)
		lgp = exp(-0.5*(dv*(iv_vec-ivzero)/(sigma_co_inner/2))**2)
		wlgp = sigma_line_inner/2.
		lgm=lgp/wlgp
	
		# Warp parameters
		#**** rx=r-10
		rx = 0.1*arange(400)
		phia = (-8.+4.*rx-0.182*rx**2)/exp(rx**2/400.)/57.3
		warpb = (9.+197.*rx-3.1*rx**2)/1.e3
		phib = (27.13-4.65*rx+0.125*rx**2)/57.3
		warpa = (-66.+150.*(rx-5.)-0.47*(rx-5.)**2)/1.e3
		warpc = (-70.+171.*(rx-5.)-5.3*(rx-5.)**2)/1.e3
		warpa[0:53] = 0.
		warpc[0:53] = 0.
	
		# Read in rotation curve
		file3 = path2+'rotcurv4.dat'
		input = open(file3,'r')
		rotation = input.readlines() #4700 lines
		rotc = zeros((len(rotation)),dtype=float)
		drot = zeros((len(rotation)),dtype=float)
		i=0
		for value in rotation:
			ha = float(value.split('\n')[0].split()[0])
			rotc[i] = float(value.split('\n')[0].split()[1])
			drot[i] = float(value.split('\n')[0].split()[2])
			i=i+1
		# free memory
		del rotation

		# Physical variables
		r_sun = 8.    #kpc
		z_sun = 0.015 #kpc
		v_sun = 210.  #km s-1
		gal_radius = 20. #kpc
		gal_thick = 1. #kpc
		dbin = 1/gal_radius #50 pc
		r_scale = 10. # radial scale [kpc]
		
		N = 760
		
		# Cuts
		v_offset = 10.     #velocity offset [10 km/s]
		lon_inner = 20.    #inner Galaxy longitude (|l|<=20) [deg]
		residual_line = 1. #threshold of remaining  CO line spectrum [K km/s]
		amp_frac = 0.2     #percentage of the peak value [x100 %]
		
		# Array definition
		true_dis = dbin*(0.5+arange(N))
		vbgr = zeros(N,dtype=float)
		
		for l in xrange(0,nlon):
			glo_deg = lon[l]
			#vlsr = vel[v]/1000.
			self.logger.info("%i) longitude: %.3f"%(l,lon[l]))
			if (abs(glo_deg) < 165.):
				glon = radians(glo_deg)
				pmean = r_sun*cos(glon)
				dismin = floor(r_sun*abs(cos(glon))/dbin)
				radmin = floor(r_sun*abs(sin(glon))/dbin)
				hzp = 12+radmin
				hzc = dismin-hzp
				hzd = dismin+hzp
				for b in xrange(0,nlat):
					#print "  %i) latitude: %.3f"%(b,lat[b])
					gla_deg = lat[b]
					glat = radians(gla_deg)
					vbgr[:] = 0.
					dimax = N-1
				
					# Calculate peculiar velocity of the sun: Equation (6)
					vpec = 10.*cos(glon)*cos(glat)+5.2*sin(glon)*cos(glat)+7.2*sin(glat)
					
					# Define intervals and heights: Equations (4)   
					z = z_sun+true_dis*sin(glat)
					proj_dis = true_dis*cos(glat)
					radi = sqrt(r_sun**2+proj_dis**2-2*r_sun*proj_dis*cos(glon)) # distance btw r_sun and proj_dis
					
					# Bissantz & Gerhard velocities
					xdir = ex_tan*cos(glon)+ex_norm*sin(glon)
					ydir = ey_tan*cos(glon)+ey_norm*sin(glon)
					xha = 100+floor(10*(x_sun_sph+proj_dis*xdir))
					yha = 99-floor(10*(y_sun_sph+proj_dis*ydir))
					ix = where((xha >= 0.) & (xha <= 199.) & (yha >= 0.) & (yha <= 199.))
					cx = size(ix[0])
					xha = xha.astype(int)
					yha = yha.astype(int)
					if(cx > 0.5):
						vbgr[ix[0]] = vrs[yha[ix[0]],xha[ix[0]]]
					
					# Remove data holes by linear interpolation
					dmax = floor(2.*r_sun*cos(glon)/dbin)
					if(dmax>6):
						vba = zeros(vbgr.shape)
						if(vbgr[0]==0.): vbgr[0] = vpec
						idx = where(vbgr[1:dmax]==0.)
						cnt = size(idx[0])
						
						while(cnt>0):
							vba[:] = 0.
							for k in xrange(0,cnt):
								ia = idx[0][k]+1
								if(vbgr[ia-1] != 0.):
									if(vbgr[ia+1] != 0.):
										vba[ia] = 0.5*(vbgr[ia-1]+vbgr[ia+1])
									else:
										vba[ia] = vbgr[ia-1]
								else:
									if(vbgr[ia+1] != 0.):
										vba[ia] = vbgr[ia+1]
								vbgr[ia] = vba[ia]
							idx = where(vbgr[1:dmax]==0.)
							cnt = size(idx[0])
					
					# Radial velocity
					# First express rotc at function of radi
					hna = floor(100.*radi)
					hnc = hna+1
					hnb = 100.*radi-hna
					hnd = 1.-hnb
					hnc = hnc.astype(int)
					hna = hna.astype(int)
					rot_curve = hnb*rotc[hnc]+hnd*rotc[hna]
					# Uncorrected radial velocity: Equation (5)
					v_lsr = (r_sun*rot_curve/radi-v_sun)*sin(glon)*cos(glat)
					#plotFunc(proj_dis,v_lsr)
					
					# Extrapolate vbgr
					if(cos(glon)>0.1):
						phn = int(round(40.*pmean-5.))
						vmean = sum( vbgr[phn-5:phn] )/6.
						vrmean = sum( v_lsr[phn-5:phn] )/6.
						vbgr[phn+1:N-1] = vmean-vrmean+v_lsr[phn+1:N-1]
					# Merging
					wradi = where((9-radi)/2 < 0.,0.,(9-radi)/2)
					veff = v_lsr+(vbgr-v_lsr)*where(wradi>1.,1.,wradi)		
					# Corrected, effective velocity: Equation (7)
					veff = movingaverage1D(veff,7)-vpec
					veff[-3:] = veff[-4]	
					#plotFunc(proj_dis,veff)

					# Weights from delta veff
					dveff = array(veff)
					dveff[-1] = fabs(veff[-2]-veff[-1])
					for i in range(0,N-1):
						dveff[i] = fabs(veff[i+1]-veff[i])
					weight_veff = zeros(veff.shape)
					# Equation (14)
					weight_veff = where((dveff+1.e-8)>dv,dv,dveff+1.e-8)
					#plotFunc(proj_dis,veff)
					#print "dveff"
					#print dveff[0:5]
					#print "weight"
					#print weight_veff[0:5]
					#exit(0)
					
					# Line spectrum
					spec = array(nvel)
					spec = Tb[:,b,l]
					spec[0] = 0.
					spec[nvel-1] = 0.
					
					rspec = fftconvolve(spec,lim,'same')
					#plotFunc(vel,spec,rspec)
					
					wco = abs(dv*sum(spec))
					wcb = wco/sigma_line

					#plotFunc(proj_dis,veff)
					
					# Start deconvolution
					while(wco > residual_line):
						
						ivpeak = argmax(rspec)
						vpeak = vel[ivpeak]
						
						amp = amp_frac*rspec[ivpeak]
						amp = where(wcb>amp,amp,wcb)
						ivlow = ivpeak-ivzero
						ivhigh = ivpeak+ivzero
						
						if(ivlow > -0.5):
							iv1 = 0 # hae = iv1; haf = iv2
							if(ivhigh < (nvel-0.5)):
								iv2 = size(iv_vec)-1
								sigma_line = sigma_co*sqrt(2.*pi)
								sigma_line_inner = sigma_co_inner*sqrt(2.*pi)
							else:
								iv2 = size(iv_vec)-ivhigh-2+nvel
								ivhigh = nvel-1
								sigma_line = fabs(dv*sum(line[iv1:iv2]))
								sigma_line_inner = fabs(dv*sum(line_inner[iv1:iv2]))
						else:
							iv1 = -ivlow
							iv2 = size(iv_vec)-1
							ivlow = 0
							sigma_line = fabs(dv*sum(line[iv1:iv2]))
							sigma_line_inner = fabs(dv*sum(line_inner[iv1:iv2]))
						
						# Finding a match between gas velocity and rotation curve
						ivgood = where((vpeak > veff) & (vpeak < (veff+dveff)))
						cnt_ivgood = size(ivgood[0])
						#print cnt_ivgood
						#print ivgood
						#print veff[ivgood],vpeak,proj_dis[ivgood]
						
						linevpeak = ones(veff.shape)*vpeak
						#plotFunc(proj_dis[0:30],veff[0:30],veff[0:30]+dveff[0:30],linevpeak[0:30])
						#exit(0)
						
						# Standard selection of locations
						
						vlist = fabs(veff[0:dimax]-vpeak)
						#plotFunc(proj_dis[0:dimax],veff[0:dimax],vlist[0:dimax],linevpeak[0:dimax])
						
						# Checking if the matching gas velocity exceeds v offset or is in the inner Galaxy 
						if((min(vlist)>v_offset) and (abs(glo_deg)<lon_inner)):
							ivmatch = argsort(vlist[hzc:hzd])+hzc
							#print,'corrected',v_match(0)
						else:
							ivmatch = argsort(vlist)
							#print,'non corrected',v_match(0)
				
						# The line signal is distributed among 8 solutions with weights
						if(cnt_ivgood < 1):
							roots = 8 # eight kinematically best-fitting location
							ilocation = zeros(roots,dtype=float)
							ilocation[0:roots] = ivmatch[0:roots]
							ika = zeros(roots,dtype=int)
							ika = (0.5*ilocation-0.25).round()
						else:
							roots = cnt_ivgood+8
							ilocation = zeros(roots,dtype=float)
							ilocation[0:cnt_ivgood] = ivgood[0][0:cnt_ivgood]
							ilocation[cnt_ivgood:roots] = ivmatch[0:8]
							ika = zeros(roots,dtype=int)
							ika = (0.5*ilocation-0.25).round()
						
						# Weights from height above plane
						wa = zeros(roots,dtype=float)
						
						for i in xrange(0,roots):
							j = ilocation[i]
							# Thickness of the gas layer - equation (15)
							sigma_z = 1.204*((0.06-0.04*radi[j]/r_sun)+0.095*(radi[j]/r_sun)**2) # pc
							zc = 0.
							# Warp in the outer region of the Galaxy (r >= 11 kpc)
							if(radi[j] > 10.):
								sphi = proj_dis[j]*sin(glon)/radi[j]
								cphi = (proj_dis[j]*cos(glon)-r_sun)/radi[j]
								nrx = floor(10.*(radi[j]-10.))
								sphia = sphi*cos(phia[nrx])-cphi*sin(phia[nrx])
								sphib = 2.*sphi*cphi*cos(phia[nrx])+(sphi**2-cphi**2)*sin(phia[nrx])
								# equation (16)
								zc = warpa[nrx]+warpb[nrx]*sphia+warpc[nrx]*sphib

							arg_wz = 0.5*((z[j]-zc)/sigma_z)**2
							dz = (true_dis[j]/sigma_z)*(pi/360)
							
							# Weights from height above plane
							weight_z = exp(-where(arg_wz<20,arg_wz,20))
							# Minimize the kinematically allowed but physically unlikely placing of gas
							weight_k = exp(-0.5*(radi[j]/r_scale)**2)
							
							wa[i] = weight_veff[j]*weight_k*(1.+dz**2*(2.*arg_wz-1)/12.)*weight_z

						wgn = 0.
						wtot = sum(wa)
						#ika.astype(int)
						
						#print roots
						#print ilocation[0:roots]
						#print ivmatch[0:roots]
						#exit(0)

						# add W_co (= sigma_line*amp) to density vector
						for i in xrange(0,roots):
							#j = ika[i]
							k = ilocation[i]
							#print "k = %.3f j = %.3f"%(k,j)
							#print "td = %.3f pd = %.3f"%(true_dis[k],proj_dis[k])
							if(radi[k] < 1.): wgn += wa[i]/wtot
							wga = wa[i]/wtot

							for a in xrange(0,annuli):
								if(proj_dis[k] > rmin[a]) and (proj_dis[k] < rmax[a]):
									cubemap[a,b,l] += wga*amp*sigma_line
									#print k,proj_dis[k],rmin[a],rmax[a],a
								#else:
									#print k
						
						#plotFunc(proj_dis[0:dimax],veff[0:dimax],vlist[0:dimax],linevpeak[0:dimax])
						#exit(0)
						
						wgo = 1.-wgn
						spec[ivlow:ivhigh] = spec[ivlow:ivhigh]-wgo*amp*line[iv1:iv2]-wgn*amp*line_inner[iv1:iv2]*sigma_line/sigma_line_inner
						rspec = fftconvolve(spec,lim,'same')
						wco = fabs(dv*sum(spec))
						wcb = wco/sigma_line
						
						#if not wco%10:
						#print wco,wcb,sigma_line
						#plotFunc(vel,rspec)
						#exit(0)
				
				#densi[b,l,379] = 0.

				
					# Integrated brightness temperature over each annulus
					#i = find_ge(self.logger,rmax,d)
					#annulus = ann_boundaries[i][0]
					#ITb[annulus,b,l] += cTb[v,b,l]
					#if not v%150:
					#		self.logger.info("(l,b) = (%i,%i) - lon=%.3f lat=%.3f d=%.1f ring=%i i=%i"%(l,b,glon,glat,d,annulus,i))
				
		# Column density
		#NHI = C*ITb*dv # [NHI] = cm-2
		#cubemap = NHI*cosdec
		#exit(0)

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
		newheader["datamin"] = ("%e"%amin(cubemap),"Min value")
		newheader["datamax"] = ("%e"%amax(cubemap),"Max value")
		newheader["object"] = ("Mosaic "+self.mosaic,self.survey+" Mosaic")
			
		results = pyfits.PrimaryHDU(cubemap, newheader)
    			
		# Output file
		self.logger.info("Write data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")
	
