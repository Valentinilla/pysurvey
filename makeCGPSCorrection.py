#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import*
import pyfits
import math
from scipy.optimize import fsolve


class makeCGPSCorrection(object):
	
	def __init__(self,mosaic,mosaicConf,utilsConf):
		
		self.logger = initLogger(mosaicConf['mosaic'], 'CGPS_ColumnDensity')
		self.species = mosaic.species
		file = ''
		path = ''
		if self.species == 'HISA':
			path = getPath(self.logger,'lustre_cgps_hisa_column_density')
			file = path+'CGPS_'+mosaicConf['mosaic']+'_HISA_column_density.fits'
		elif self.species == 'HI':
			path = getPath(self.logger,'lustre_cgps_hi_column_density')
			file = path+'CGPS_'+mosaicConf['mosaic']+'_HI_unabsorbed_column_density.fits'
		elif self.species == 'HI+HISA':
			path = getPath(self.logger,'lustre_cgps_hi_column_density')
			file = path+'CGPS_'+mosaicConf['mosaic']+'_HI+HISA_column_density.fits'
		elif self.species == 'CO':
			path = getPath(self.logger,'lustre_cgps_co_column_density')
			file = path+'CGPS_'+mosaicConf['mosaic']+'_H2_column_density.fits'
		checkForFiles(self.logger,[file],existence=True)
		
		self.logger.info("Open file and get data...")
		
		# Array to store results
		N = zeros((mosaic.msc_lat,mosaic.msc_lon),dtype=float)
		# Usefull velocity channels
               	kmin = 18
               	kmax = 271
		# Used to skip calculation (see below) - do not change it!
		flag = True

		if self.species == 'HI' or self.species == 'HISA' or self.species == 'HI+HISA': 
			# Get emission data and velocity interval
			Tb = mosaic.observation[0,:,:,:]
			dv = fabs(mosaic.msc_del_vel) # [velocity] = km s-1
			# Setting the negative/0-values to Tcmb 
			Tb = where( (Tb<0.) | (Tb==0.),float(utilsConf['tcmb']),Tb)
		
			# If the two column density files already exist, than add them up and skip the calculation
			if self.species == 'HI+HISA':
				path = getPath(self.logger,'lustre_cgps_hisa_column_density')
				file1 = path+'CGPS_'+mosaicConf['mosaic']+'_HISA_column_density.fits'
				
				path = getPath(self.logger,'lustre_cgps_hi_column_density')
				file2 = path+'CGPS_'+mosaicConf['mosaic']+'_HI_unabsorbed_column_density.fits'
				
				if os.path.isfile(file1) and os.path.isfile(file2):
					flag = False
					data1, header = pyfits.getdata(file1, 0, header = True)
					data2, header = pyfits.getdata(file2, 0, header = True)
					N = data1+data2
		
			if self.species == 'HI' or self.species == 'HI+HISA' and flag == True:
			
				NHI = zeros((mosaic.msc_lat,mosaic.msc_lon),dtype=float)
				
				self.logger.info("Initializing parameters...")
				#Tmp = Tb[:928,:928]
				Ts = 150#ceil(amax(Tmp))	  # [Excitation (or Spin) Temperature] = K (150)
				
				#index = where(Tb>=Ts)
				#Tb[index] = 0.999*Ts # Tb must be < Ts, otherwise problems with log
				
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
				
				Tfunc[Tfunc<1.] = float(utilsConf['tcmb'])# <------ TO JUSTIFY
				cTb = log(Tfunc) * Ts
				# Integrated brightness temperature over velocity (axis = 0)
				ITb = sum(cTb[kmin:kmax,:,:],axis=0)
				# Column density
				NHI = float(utilsConf['c'])*ITb*dv # [NHI] = cm-2
				N = NHI*cos( (pi/180.)*mosaic.lat_array)
			
			if self.species == 'HISA' or self.species == 'HI+HISA' and flag == True:
			
				NHISA = zeros((mosaic.observation.shape),dtype=float)
				
				# Get HI continuum data
				pathc = getPath(self.logger, 'cgps_hi_continuum')
				continuum = pathc+'CGPS_'+mosaicConf['mosaic']+'_1420_MHz_I_image.fits'
				data, header = pyfits.getdata(continuum, 0, header = True)
				data[isnan(data)] = 0.
				Tc = data[0,0,:,:]
				
				# Get HISA data
				path = getPath(self.logger,'cgps_hisa_dat')
				amplitude = path+'CGPS_'+mosaicConf['mosaic']+'_HISA.dat'
				input = open(amplitude,'r')
				lines = input.readlines()
				
				self.logger.info("Initializing parameters...")
				self.logger.info("1) dV = %.2f km/s"%dv)
				self.logger.info("2) Tb(min) = %.2f K, Tb(max) = %.2f K"%(amin(Tb),amax(Tb)))
				self.logger.info("3) Tc(min) = %.2f K, Tc(max) = %.2f K"%(amin(Tc),amax(Tc)))
				
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
					glon = mosaic.lon_array[k]
					glat = mosaic.lat_array[l+160] #<--- 160? (it's in the martin's idl code)
					vlsr = mosaic.vel_array[m]
					#print k,l,glon,glat,vlsr
					#exit(0)
					d = rotCurveMPohl(self.logger,glon,glat,vlsr) #kpc
					theta = mosaic.msc_del_lat*pi/180 #rad
					ds = d*tan(theta)*1e3 #pc #(theta*pi/10800) from arcm to rad
					
					A1 = float(utilsConf['pc2cm'])*float(utilsConf['poverk'])
					A2 = float(utilsConf['fn'])*ds/(float(utilsConf['c'])*dv)
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
							Tfunc = 1.#float(utilsConf['tcmb'])# <------ TO JUSTIFY
						f1 = log(Tfunc)-tt
						# Equation (9)
						ttfunc = A/tt	
						if ttfunc<0.:
							ttfunc = 0.# <------ TO JUSTIFY
						f2 = sqrt(ttfunc)-T
													
						return array([f1, f2], dtype=float)
					
					(tau,Ts),infodict,ier,mesg = fsolve(equations,init, full_output=1)
					#plotFunc(tau,Ts)
					TsMin = float(utilsConf['tcmb'])
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
						NHISA[0,m,l,k] = tau*Ts*dv*float(utilsConf['c'])
						NHISA[0,m,l,k] = NHISA[0,m,l,k]*cos( (pi/180.)*mosaic.lat_array[l])
					if not solution:
						#print "Could not find a valid solution:"
						#print mesg
						NHISA[0,m,l,k] = 0.
				
				self.logger.info("Calculating NHI...")
				# Corrected column density
				N = N + sum(NHISA[0,kmin:kmax,:,:],axis=0)

		if self.species == 'CO':
			# Get emission data and velocity interval
			Tb = mosaic.observation[:,:]
			N = Tb*float(utilsConf['xfactor'])
		
		# Store results
		newheader = pyfits.Header()
		newheader.update(key="ctype1", value="GLON-CAR", comment="Coordinate type")
		newheader.update(key="crval1", value=mosaic.msc_ref_lon, comment="Galactic longitude of reference pixel")
		newheader.update(key="crpix1", value=mosaic.msc_ind_lon, comment="Reference pixel in lon")
		newheader.update(key="cdelt1", value=mosaic.msc_del_lon, comment="Longitude increment")
		newheader.update(key="crota1", value=mosaic.msc_rot_lon, comment="Longitude rotation")
		newheader.update(key="cunit1", value="deg", comment="Unit type")

		newheader.update(key="ctype2", value="GLAT-CAR", comment="Coordinate type")
		newheader.update(key="crval2", value=mosaic.msc_ref_lat, comment="Galactic latitude of reference pixel")
		newheader.update(key="crpix2", value=mosaic.msc_ind_lat, comment="Reference pixel in lat")
		newheader.update(key="cdelt2", value=mosaic.msc_del_lat, comment="Latitude increment")
		newheader.update(key="crota2", value=mosaic.msc_rot_lat, comment="Latitude rotation")
		newheader.update(key="cunit2", value="deg", comment="Unit type")

		newheader.update(key="bunit", value="cm-2", comment="Map units")

		newheader.update(key="datamin", value="%e"%amin(N))
		newheader.update(key="datamax", value="%e"%amax(N))

		newheader.update(key="object", value="CGPS Mosaic %s"%mosaicConf['mosaic'], comment="GCPS Mosaic")

		results = pyfits.PrimaryHDU(N, newheader)
		
		# Output file
		self.logger.info("Write data to a fits file in...")
		results.writeto(file, output_verify='fix')
		self.logger.info("%s"%path)
		self.logger.info("Done")

