#!/usr/bin/env python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

import os
import sys
import logging
import ConfigParser

# lmfit info: http://newville.github.com/lmfit-py/index.html
from lmfit import minimize, Parameters
from lmfit.printfuncs import*

from numpy import *
from numpy.lib.stride_tricks import as_strided

from scipy.signal import fftconvolve
from scipy.optimize import fsolve
import scipy.ndimage as ndimage

import pyfits

# Avoid math module if possible because many functions 
# conflict with those in numpy (math functions operate in 
# 1D arrays, numpy's in multidimensional arrays)
# If it is strictly necessary to use math then 
# import functions explicitly, i.e. 
# from math import func1,func2,...

#################################################
# GLOBAL VARIABLES
#################################################
glob_Tb = 'brightness temperature'
glob_ITb = 'integrated brightness temperature'
glob_N = 'column density'

#################################################
#	START GENERAL FUNCTIONS
#################################################
def getMosaicCoordinate(surveyLogger,obs,survey,lon,lat,side):
	"""
	This function allows to generate 'mosaics' like those
	defined in CGPS, SGPS, and VGPS.
	Input parameters:
	- longitude: galactic longitude of the mosaic's center
	- latitude:  galactic latitude of the mosaic's center
	- side:      length in degrees of the mosaic's side
	Outout parameters:
	- l: list of longitude coordinates in pixels (l = [l1,l2])
	- b: list of latitude coordinates in pixels (b = [b1,b2])
	- sign: list (sign = [lon_sign,lat_sign])
	"""
	
	lon1 = lon+side
	lon2 = lon-side
	lat1 = lat-side
	lat2 = lat+side
	
	# Check if the coordinates in .cfg are inside the mosaic object
	obj_lon_max = obs.msc_ref_lon+((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
	obj_lon_min = obs.msc_ref_lon-((obs.msc_lon-1.)*abs(obs.msc_del_lon))/2.
	if survey == 'LAB':
		obj_lat_max = obs.msc_ref_lat+((obs.msc_lat-1.)*obs.msc_del_lat)
		obj_lat_min = obs.msc_ref_lat
	else:
		obj_lat_max = obs.msc_ref_lat+((obs.msc_lat-1.)*obs.msc_del_lat)/2.
		obj_lat_min = obs.msc_ref_lat-((obs.msc_lat-1.)*obs.msc_del_lat)/2.

	lon_max = max(lon1,lon2)
	lon_min = min(lon1,lon2)
	lat_max = max(lat1,lat2)
	lat_min = min(lat1,lat2)
	
	if not survey == 'Galprop':
		if lon_max>obj_lon_max or lon_min<obj_lon_min:
			surveyLogger.critical("The longitude within the .cfg file")
			surveyLogger.critical("doesn't match that of the mosaic!")
			sys.exit(0)
		if lat_max>obj_lat_max or lat_min<obj_lat_min:
			surveyLogger.critical("The latitude within the .cfg file")
			surveyLogger.critical("doesn't match that of the mosaic!")	
			sys.exit(0)

	l1 = int(round(obs.msc_ind_lon+(lon1-obs.msc_ref_lon)/obs.msc_del_lon))
	l2 = int(round(obs.msc_ind_lon+(lon2-obs.msc_ref_lon)/obs.msc_del_lon))
	b1 = int(round(obs.msc_ind_lat+(lat1-obs.msc_ref_lat)/obs.msc_del_lat))
	b2 = int(round(obs.msc_ind_lat+(lat2-obs.msc_ref_lat)/obs.msc_del_lat))

	l = [l1,l2]
	b = [b1,b2]

	# Get the sign of delta_l and delta_b
	lsign = int((l[1]-l[0])/fabs(l[1]-l[0]))
	bsign = int((b[1]-b[0])/fabs(b[1]-b[0]))
	sign = [lsign,bsign]
	
	# Order the indexes
	if l1>l2:
		l = [l2,l1]
	if b1>b2:
		b = [b2,b1]
	
	length = rint(side/fabs(obs.msc_del_lat)*2)
	
	# Make an equal sides (Lx=Ly) mosaic
	if length%2 == 0:
		length = length+1
	if (l[1]-l[0]) < length:
		l[0] = l[0]-1
	if (b[1]-b[0]) < length:
		b[0] = b[0]-1
	
	return l,b,sign

#################################################
#	START CONVERTING FUNCTIONS
#################################################
def ga2equ(lon,lat,ref='J2000'):
	"""
	Convert Galactic to Equatorial coordinates (J2000.0)
	Input: [l,b] in decimal degrees
	Returns: [ra,dec] in decimal degrees
	Source:
	- Book: "Practical astronomy with your calculator" (Peter Duffett-Smith)
	- Wikipedia "Galactic coordinates"
	Tests (examples given on the Wikipedia page):
	>>> ga2equ([0.0, 0.0]).round(3)
	array([ 266.405, -28.936])
	>>> ga2equ([359.9443056, -0.0461944444]).round(3)
	array([ 266.417, -29.008])
	"""
	l = radians(lon) # == ga*pi/180.
	b = radians(lat)

	if ref=='J2000':
		# North galactic pole (J2000)
		pole_ra = radians(192.859508)
		pole_dec = radians(27.128336)
		posangle = radians(122.932-90.0)
	if ref=='B1950':
		# North galactic pole (B1950)
		pole_ra = radians(192.25)
		pole_dec = radians(27.4)
		posangle = radians(123.0-90.0)
	
	
	y = (cos(b)*cos(l-posangle))
	x = (sin(b)*cos(pole_dec)-cos(b)*sin(pole_dec)*sin(l-posangle))
	
	ra_rad = atan2(y,x) + pole_ra
	dec_rad = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
	
	ra_deg = degrees(ra_rad)-360*int(degrees(ra_rad)/360)
	dec_deg = degrees(dec_rad)

	return ra_deg,dec_deg

def eq2gal(ra,dec,ref='J2000'):
	"""
	Convert Equatorial to Galactic coordinates (J2000.0/B1950)
	Input: [l,b] in decimal degrees
	Returns: [ra,dec] in decimal degrees
	"""
	ra = radians(ra)
	dec = radians(dec)
	
	if ref=='J2000':
		# North galactic pole (J2000)
		pole_ra = radians(192.859508)
		pole_dec = radians(27.128336)
		posangle = radians(122.932-90.0)
	if ref=='B1950':
		# North galactic pole (B1950)
		pole_ra = radians(192.25)
		pole_dec = radians(27.4)
		posangle = radians(123.0-90.0)
	
	b = asin(cos(dec)*cos(pole_dec)*cos(ra-pole_ra)+sin(dec)*sin(pole_dec))
	
	A = sin(dec)-sin(b)*sin(pole_dec)
	B = cos(dec)*sin(ra-pole_ra)*cos(pole_dec)
	l = atan2(A,B) + posangle
	
	return degrees(l),degrees(b)

def sex2dec(deg,min,sec):
	"""
	Convert from sexagesimal (i.e., 182 deg 31'27'') to decimal degrees (i.e., 182.524167 deg)
	"""
	A = fabs(sec/60.)
	B = (fabs(min)+A)/60.
	newdeg = (fabs(deg)+B)
	if deg<0. or min<0. or sec<0.:
		newdeg = -1*newdeg
	return newdeg

def dec2sex(deg):
	"""
	Convert from decimal (i.e., 182.524167 deg) to sexagesimal degrees, minutes, seconds (i.e., 182 deg 31'27'')
	"""
	A = fabs(deg) # unsigned decimal
	B = A*3600 # total seconds
	C = round((B%60.),2) # seconds (2 dp)
	if C==60.:
		sec=0.  # corrected seconds
		D=B+60. # corrected remainder
	else:
		sec=C
		D=B
	min = int(D/60.)%60. # minutes
	newdeg = (deg/A)*int(D/3600.) # signed degrees
		
	return int(newdeg),int(min),int(sec)

#################################################
#	END CONVERTING FUNCTIONS
#################################################

def get_nth_maxvalue(a, nth):
	b = a.flatten()
	res_sort = sort(b)
	c = res_sort[-nth:]
	return c[0]

def getSign(number, string=False):
	if string:
		Nsign = int(number/fabs(number))
		Ssign = ('%s'%Nsign).split('1')[0]
		if Ssign == '':
			return '+'
		elif Ssign == '-':
			return ''
	else:
		return int(number/fabs(number))

def test(func):
	#print test(get_nth_maxvalue(residual,self.HIGH))
	#print test(get_nth_maxvalue2(residual,self.HIGH))
	from timeit import Timer
	t = Timer(lambda: func)
	return t.timeit()

def movingaverage1D(array,w):
        window = ones(w,dtype=int)/float(w)
        return fftconvolve(array,window,'same')

def spatialAverage1D(array,i,j,n_spec): #7px
	# smoothing with spatial average
  	box_size = fix(n_spec/2) #3
	s = array.shape[0]
	smoothed = zeros(s,dtype=float)

	if i < box_size or j < box_size:
		for k in range(18,s-1):
			smoothed[k] = array[k,j,i]
	else:
		for k in range(18,s-1):
			j1 = (j-box_size)-1
			j2 = (j+box_size)
			i1 = (i-box_size)-1
			i2 = (i+box_size)
			smoothed[k] = mean(array[k,j1:j2,i1:i2])
			#print "(%i, %i) [%i,%i:%i,%i:%i] mean = %f" %(j,i,k,j1,j2,i1,i2,smoothed[k])
			#print array[k,j1:j2,i1:i2].shape
	return smoothed

def extrap(x, xp, yp):
	"""
	numpy.interp function with linear extrapolation
	"""
	y = interp(x, xp, yp)
	y[x < xp[0]] = yp[0] + (x-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
	y[x > xp[-1]] = yp[-1] + (x-xp[-1]) * (yp[-1]-yp[-2]) / (xp[-1]-xp[-2])
	return y

# Moment Masking Routins
def getRMS(surveyLogger,T):
	
	sigma_array = zeros(3,dtype=float)
	# Compute the standard deviation along the specified axis:
	# std = sqrt(mean(abs(x - x.mean())**2))
	sigma_array[0] = mean(T[235,:,:]**2)#std(T[235,:,:])
	sigma_array[1] = mean(T[240,:,:]**2)#std(T[240,:,:])
	sigma_array[2] = mean(T[245,:,:]**2)#std(T[245,:,:])
	rms = amin(sigma_array)
	
	surveyLogger.info("RMS = %s"%rms)
	return rms

		
#################################################
#	END GENERAL FUNCTIONS
#################################################

#################################################
#	START ROTATION CURVE
#################################################
def rotCurveMPohl(surveyLogger,glo_deg,gla_deg,vlsr):
	"""
	Rotation curve of the Galaxy - M.Pohl, P.Englmaier, and N.Bissantz
	(The Astrophysical Journal, 677:283-291, 2008)
	Limits: galactic longitude < |165 deg|, galactic latitude < |5 deg|.
	"""
	path = getPath(surveyLogger,'rotcurve_mpohl')	
	
	# Read in SPH model results
	# Get Bissantz's data
	file = path+'testvr1.dat'
	input = open(file,'r')
	bissantz = input.read().split()
	n1 = 200
	vrs = array(bissantz).reshape(n1,n1)
	vrs = vrs.astype(float)
	#print vrs[100,98]  #-79.1474
	#print vrs[100,99]  #-56.3561
	#print vrs[100,100] #-25.6225
	
	R = 0.5
	xa = -10.+0.1*(R+arange(n1))
	ya = -10.+0.1*(R+arange(n1))
	rb = zeros((n1,n1),dtype=float)
	for i in range(0,n1-1):
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
	
	# Read in rotation curve
	file = path+'rotcurv4.dat'
	input = open(file,'r')
	rotation = input.readlines() #4700 lines
	rotc = zeros((len(rotation)),dtype=float)
	drot = zeros((len(rotation)),dtype=float)
	i=0
	for value in rotation:
		ha = float(value.split('\n')[0].split()[0])
		rotc[i] = float(value.split('\n')[0].split()[1])
		drot[i] = float(value.split('\n')[0].split()[2])
		i=i+1
	
	# Array definition
	true_dis = 0.05*(0.5+arange(760))
	vbgr = zeros((760),dtype=float)
	
	r_sun = 8.    #kpc
	z_sun = 0.015 #kpc
	v_sun = 210.  #km s-1
	
	#glo_deg = 30.
	#glat_deg = 0.
	#vlsr = -10.

	if (abs(glo_deg) < 165.):
		glon = radians(glo_deg)
		glat = radians(gla_deg)
		pmean = r_sun*cos(glon)
		dismin = floor(r_sun*abs(cos(glon))/0.05)
		radmin = floor(r_sun*abs(sin(glon))/0.05)
		hzp = 12+radmin
		hzc = dismin-hzp
		hzd = dismin+hzp
		vbgr[:] = 0.
		
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
		
		# Remove data holes
		dmax = floor(16.*cos(glon)/0.05)
		if(dmax>6):
			vba = zeros(vbgr.shape)
			if(vbgr[0]==0.):
				vbgr[0] = vpec
			ib = where(vbgr[1:dmax]==0.)
			ca = size(ib[0])
			while(ca>0):
				vba[:] = 0.
				for k in range(0,ca):
					ia = ib[0][k]+1
					if(vbgr[ia-1] != 0.):
						if(vbgr[ia+1] != 0.):
							vba[ia] = 0.5*(vbgr[ia-1]+vbgr[ia+1])
						else:
							vba[ia] = vbgr[ia-1]
					else:
						if(vbgr[ia+1] != 0.):
							vba[ia] = vbgr[ia+1]
					vbgr[ia] = vba[ia]
				ib = where(vbgr[1:dmax]==0.)
				ca = size(ib[0])
		
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
		#if(abs(cos(glon))>0.1):
			#print 'Here'
			phn = int(round(40.*pmean-5.))
			vmean = sum( vbgr[phn-5:phn] )/6.
			vrmean = sum( v_lsr[phn-5:phn] )/6.
			vbgr[phn+1:759] = vmean-vrmean+v_lsr[phn+1:759]
		# Merging
		wradi = where((9-radi)/2 < 0.,0.,(9-radi)/2)
		vr = v_lsr+(vbgr-v_lsr)*where(wradi>1.,1.,wradi)		
		# Corrected, effective velocity: Equation (7)
		vr = movingaverage1D(vr,7)-vpec
		vr[-3:] = vr[-4]
		#extrap(vr, proj_dis, yp)
		#plotFunc(proj_dis,vr)
		
		try:
			i = 0
			diff=[]
			for vel in vr:
				if abs(vel-vlsr) < 10.: # km/s
					diff.append([abs(vel-vlsr),i])
					#print "1) vlsr = %s, v = %s, r = %s, i = %s"%(vlsr,vel,proj_dis[i],i)
				i=i+1			
			radius = proj_dis[min(diff)[1]]
			return radius # kpc
		except ValueError:
			surveyLogger.critical("The rotation curve doesn't contain the vlsr value of the mosaic!!")
			surveyLogger.critical("i) [vrot_min, vrot_max] = [%.2f,%.2f], vlsr_msc = %.2f"%(amin(vr),amax(vr),vlsr))
			sys.exit(0)	

#################################################
# END ROTATION CURVE
#################################################

#################################################
# START PLOTTING FUNCTIONS
#################################################
def plotFunc(x,*func):
	import matplotlib.pyplot as plt
	n = len(func)

	if n==1:
		plt.plot(x,func[0],label='func1',color='black',lw=1)#,'x',xi2,yi2)
		plt.legend()# ('func1'),loc='upper left', shadow=False, fancybox=True)
	if n==2:
		plt.plot(x,func[0],'o',x,func[1])
		plt.legend( ('func1','func2'),loc='upper left',shadow=False,fancybox=True)
	if n==3:
		plt.plot(x,func[0],x,func[1],x,func[2])
		plt.legend( ('func1','func2','func3'),loc='upper left',shadow=False,fancybox=True)
	if n==4:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3])
		plt.legend( ('func1','func2','func3','func4'),loc='upper left',shadow=False,fancybox=True)
	if n==5:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4])
		plt.legend( ('f1','f2','f3','f4','f5'),loc='upper left',shadow=False,fancybox=True)
	if n==6:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4],x,func[5])
		plt.legend( ('f1','f2','f3','f4','f5','f6'),loc='upper left',shadow=False,fancybox=True)
	if n==7:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4],x,func[5],x,func[6])
		plt.legend( ('f1','f2','f3','f4','f5','f6'),loc='upper left',shadow=False,fancybox=True)
	if n==8:
		plt.plot(x,func[0],x,func[1],'o',x,func[2],x,func[3],x,func[4],x,func[5],x,func[6],'--',x,func[7],'--')
		#unabsorbed,smoothed,residual,HISA_narrow,HISA_broad,HISA_merged
		plt.legend( ('unabsorbed','smoothed','residual','narrow','broad','merged'),loc='upper right',shadow=False,fancybox=True)
	if n>8:
		print "Max number of functions is 8, %i found!"%n
		exit(0)
	#plt.axis([0,1.1,3.0,5.5])
	plt.xlabel('$v_{LSR}(km/s)$')
	plt.ylabel('$T_b(K)$')
	plt.title('Velocity profile')
	plt.grid(True)
	plt.show()

def plotFunc2D(x,y,f):

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#ax.plot_wireframe(x, y, f, color='0.3', lw=0.5)
	ax.plot_surface(x, y, f,  rstride=4, cstride=4, color='b')
	#cax = ax.imshow(funcs, interpolation='nearest')
	ax.set_title('2D Gaussian')
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])
	plt.show()

def polarMap():
	from matplotlib.pyplot import figure, show
	
	fig = figure()
	ax = fig.add_subplot(111, polar=True)
	#r = arange(0,1,0.001)
	#theta = 2*2*pi*r
	r = 0.2249
	theta = 122.98
	line, = ax.plot(theta, r, color='#ee8d18', lw=3)
	
	#ind = 800
	thisr, thistheta = r,theta#r[ind], theta[ind]
	ax.plot([thistheta], [thisr], 'o')
	ax.set_rmax(20.0)
	ax.annotate('a polar annotation',
		xy=(thistheta, thisr),  # theta, radius
		xytext=(0., 35),    # fraction, fraction
		textcoords='figure fraction',
		arrowprops=dict(facecolor='black', shrink=0.05),
		horizontalalignment='left',
		verticalalignment='bottom',
	)
	show()

def plotNvsTs(surveyLogger,obs,utilsConf,plot,lon,lat):
	"""
	Compute the column density as a function of Ts for a certain LOS.
	The output is a datafile.
	"""
	if plot == 'NH vs Ts':
		NHI = zeros((obs.msc_lat,obs.msc_lon),dtype=float)
		if obs.survey == 'LAB':
			Tb = obs.observation[:,:,:]
		else:
			Tb = obs.observation[0,:,:,:]
		
		#Ts = float(utilsConf['tspin'])	  # [Excitation (or Spin) Temperature] = K (150)
		Tbg = float(utilsConf['tcmb'])	  # [Cosmic Microwave Background (CMB)] = K
		dv = fabs(obs.msc_del_vel)	  # [velocity] = km s-1
		C = float(utilsConf['c'])	  # [costant] = cm-2
		
		pxsize = obs.msc_del_lat*(pi/180.)	 # [pixel size] = rad
		Apx = power(pxsize, 2)	 		 # [Area pixel] = sr
		cosdec = cos((pi/180.)*obs.lat_array)     # [lat. correction] = rad
		
		Tb = where( (Tb<0.)|(Tb==0.),Tbg,Tb)
		# Usefull velocity channels
		kmin = 0
		kmax = 0
		if obs.survey == 'CGPS':
	        	kmin = 18
	        	kmax = 271
		if obs.survey == 'SGPS':
	        	kmin = 1
        		kmax = 410
		
		l = int(round(obs.msc_ind_lon-1.+(lon-obs.msc_ref_lon)/obs.msc_del_lon))
		b = int(round(obs.msc_ind_lat-1.+(lat-obs.msc_ref_lat)/obs.msc_del_lat))
		#print l,b
		
		datafile = obs.survey+'_'+obs.mosaic+'_'+obs.species+'_nh_vs_ts.dat'
		f = os.path.isfile(datafile)
		
		if f == False:	
			outfile = open(datafile,"w")
			string = ""
			
			surveyLogger.info("Calculating NHI for different Ts at the LOS (%s,%s) deg"%(lon,lat))
			
			for Ts in range(100,201,5):
				NHI[:,:] = 0.
				# Without the continuum component
				index = where(Tb>=Ts)
				Tb[index] = 0.999*Ts # Tb must be < Ts, otherwise problems with log
				
				cTb = log( (Ts)/(Ts-Tb) ) * Ts
				ITb = sum(cTb[kmin:kmax,:,:],axis=0)
				NHI = C*ITb*dv
				
				string += "%s\t%s\n" %(Ts,NHI[b,l]*cosdec[b])
			outfile.write(string)
			outfile.close()
		
			surveyLogger.info("Saving data on '%s'"%datafile)

		surveyLogger.info("Opening data of NHI as a function of Ts at the LOS (%s,%s) deg"%(lon,lat))
		data = open(datafile,"r")
		lines = data.readlines()
		ts = []
		nh = []
		for i in range(len(lines)):
			ts.append(lines[i].split()[0])
			nh.append(lines[i].split()[1])
		
		import matplotlib.pyplot as plt
		plt.plot(ts,nh,label='Column Density',color='blue',lw=1)#,'x',xi2,yi2)
		plt.legend(loc='upper right', shadow=False, fancybox=True)
		plt.xlabel('$T_s (K)$')
		plt.ylabel('$N_H (cm^{-2})$')
		plt.title(obs.survey+': '+obs.mosaic+' ('+obs.species+'), LOS = (%s,%s) deg'%(lon,lat))
		plt.grid(True)
		plt.show()
		#fig = plt.figure()
		#ax = fig.add_subplot(2,1,1)
		#line, = ax.plot(ts,nh, color='blue', lw=2)
		#ax.set_yscale('log')
		#plt.show()

#################################################
# END PLOTTING FUNCTIONS
#################################################


class FileNotFound: pass
class CommandNotFound: pass

def getPath(surveyLogger, key="cgps_hi"):
		
	path = False
	# Rot.Curve
	if key == 'rotcurve_mpohl':
		path = True
		return '/afs/ifh.de/group/that/work-sf/survey/rotcurve/'
	
	# TODO
	#survey = survey_name
	#sur = survey_name.lower()
	#
	#if key == 'lustre_'+sur:
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/'
	#if key == 'lustre_'+sur+'_hi':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/HI/'
	#if key == 'lustre_'+sur+'_hi_column_density':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/HI/col_den/'
	#if key == 'lustre_'+sur+'_hisa':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/HISA/'
	#if key == 'lustre_'+sur+'_hisa_column_density':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/HISA/col_den/'
	#if key == 'lustre_'+sur+'_co':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/CO/'
	#if key == 'lustre_'+sur+'_co_column_density':
	#	path = True
	#	return '/lustre/fs4/group/that/sf/Surveys/'+survey+'/CO/col_den/'

	# Galprop
	if key=="galprop_hi":
		path = True
		return '/afs/ifh.de/group/that/soft/opt/Galprop/FITS/'
	if key=="galprop_co":
		path = True
		return '/afs/ifh.de/group/that/soft/opt/Galprop/FITS/'
	
	if key=="lustre_galprop":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/'
	if key=="lustre_galprop_hi":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/HI/'
	if key=="lustre_galprop_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/HI/col_den/'
	if key=="lustre_galprop_hisa":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/HISA/'
	if key=="lustre_galprop_hisa_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/HISA/col_den/'
	if key=="lustre_galprop_co":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/CO/'
	if key=="lustre_galprop_co_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Galprop/CO/col_den/'

	# Dame
	if key=="dame_co":
		path = True
		return '/afs/ifh.de/group/that/data/DiffuseEmission/co/'
	if key=="lustre_dame":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Dame/'
	if key=="lustre_dame_co_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/Dame/col_den/'
	
	# LAB
	if key=="lab_hi":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/hi/LAB/'
	if key=="lustre_lab":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/LAB/'
	if key=="lustre_lab_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/LAB/col_den/'

	# CGPS
	if key=="cgps_hi":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/hi/CGPS/'
	if key=="cgps_hi_continuum":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/continum/cgps/'
	if key=="cgps_hisa_dat":
		path = True
		return '/afs/ifh.de/group/that/data/HISA/cgps/results/'
	if key=="cgps_co":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/co/cgps/'

	if key=="lustre_cgps":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/'
	if key=="lustre_cgps_hi":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/HI/'
	if key=="lustre_cgps_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/HI/col_den/'
	if key=="lustre_cgps_hisa":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/HISA/'
	if key=="lustre_cgps_hisa_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/HISA/col_den/'
	if key=="lustre_cgps_co":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/CO/'
	if key=="lustre_cgps_co_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/CGPS/CO/col_den/'

	# SGPS
	if key=="sgps_hi":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/hi/SGPS/'
	if key=="sgps_hi_continuum":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/continum/sgps/'
	if key=="sgps_hisa_dat":
		path = True
		return '/afs/ifh.de/group/that/data/HISA/sgps/results/'

	if key=="lustre_sgps":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/SGPS/'
	if key=="lustre_sgps_hi":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/SGPS/HI/'
	if key=="lustre_sgps_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/SGPS/HI/col_den/'
	if key=="lustre_sgps_hisa":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/SGPS/HISA/'
	if key=="lustre_sgps_hisa_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/Surveys/SGPS/HISA/col_den/'
		
	if(not path):
		surveyLogger.critical("Path '%s' doesn't exist."%key)
		raise FileNotFound

#--------------------------------
# ERROR MSG FUNCTION
#--------------------------------
def typeErrorMsg(surveyLogger,type,typeEntry='HI'):

	if typeEntry=='HI' or typeEntry=='HISA' or typeEntry=='CO':
		msg = "Allowed types are: '"+glob_Tb+"', and '"+glob_N+"'. Your entry is: '"+type+"'."

	elif typeEntry=='WCO':
		msg = "Allowed types are: '"+glob_ITb+"', '"+glob_N+"'. Your entry is: '"+type+"'."

	elif typeEntry=='HI+HISA' or typeEntry=='HI+CO' or typeEntry=='HI+HISA+CO':
		msg = "Allowed type is only '"+glob_N+"'. Your entry is: '"+type+"'."
	
	surveyLogger.critical(msg)

def checkToGetData(surveyLogger,f1,f2,f3=None):
	
	file1 = os.path.isfile(f1)
	file2 = os.path.isfile(f2)
	file3 = os.path.isfile(f3)

	flag1,flag2,flag3 = True,True,True

	if f3 == None:
		# if both files already exist		
		if file1 and file2:
			flag1, flag2 = False, False
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			N = data1+data2
		#file1 already exists
		if file1 and not file2:
			flag1 = False
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			N = data1
		# file2 already exists
		if file2 and not file1:
			flag2 = False
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			N = data2
		return N,flag1,flag2
	else:
		# All of the files already exist
		if(file1 and file2 and file3):
			flag1,flag2,flag3 = False,False,False
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			data3, header3 = pyfits.getdata(f3, 0, header = True)
			N = data1+data2+data3
		# file1 and file2 already exist	
		if(file1 and file2) and not file3:
			flag1,flag2,flag3 = False,False,True
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			N = data1+data2
		# file1 and file3 already exist	
		if(file1 and file3) and not file2:
			flag1,flag2,flag3 = False,True,False
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			data3, header3 = pyfits.getdata(f3, 0, header = True)
			N = data1+data3
		# file2 and file3 already exist	
		if(file2 and file3) and not file1:
			flag1,flag2,flag3 = True,False,False
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			data3, header3 = pyfits.getdata(f3, 0, header = True)
			N = data2+data3
		# file1 already exists
		if file1 and not(file2 and file3):
			flag1,flag2,flag3 = False,True,True
			data1, header1 = pyfits.getdata(f1, 0, header = True)
			N = data1
		# file2 already exists
		if file2 and not(file1 and file3):
			flag1,flag2,flag3 = True,False,True
			data2, header2 = pyfits.getdata(f2, 0, header = True)
			N = data2
		# file3 already exists
		if file3 and not(file1 and file2):
			flag1,flag2,flag3 = True,True,False
			data3, header3 = pyfits.getdata(f3, 0, header = True)
			N = data3
		return N,flag1,flag2,flag3

def checkForFiles(surveyLogger, fileList, existence=False):
	"""
	Checks for the existence of needed files in the list.
	"""
	for filename in fileList:
		if not os.path.exists(filename) and not existence:
			surveyLogger.critical(filename+" doesn't exist.")
			raise FileNotFound
		elif os.path.exists(filename) and existence:
			surveyLogger.critical(filename+" already exists.")
			raise FileNotFound

def checkForCommand(surveyLogger, commandList):
	"""
	Checks for the existence of a certain command.
	"""
	for command in commandList:
		cmd = "which -s " + command + " > " + os.devnull + " 2>&1"
		retcode = os.system(cmd)
	
		if(retcode):
			surveyLogger.critical("unix command "+command+" not found.")
			raise CommandNotFound     
        
def writeConfig(surveyLogger, surveyDictionary, mosaicDictionary = {}, utilsDictionary = {}, spectralDictionary = {}, spatialDictionary = {}):
	"""
	Writes all of the needed information to the config file called <mosaicname>.cfg
	"""
	configfilename = surveyDictionary['survey']+'_'+mosaicDictionary['mosaic']

	config = ConfigParser.RawConfigParser()
	config.read(configfilename+'.cfg')
	if(not config.has_section('survey')):
		config.add_section('survey')

	for variable, value in surveyDictionary.iteritems():
		config.set('survey', variable, value)
	surveyLogger.info('wrote common config to '+configfilename+'.cfg.')

	if(mosaicDictionary):
		if(config.has_section('mosaic')):
			surveyLogger.info("mosaic config exists, overwriting...")        
		else:
			config.add_section('mosaic')
		for variable, value in mosaicDictionary.iteritems():
			config.set('mosaic', variable, value)
		surveyLogger.info("wrote mosaic config to "+configfilename+".cfg.")

	if(utilsDictionary):
		if(config.has_section('utils')):
			surveyLogger.info("utils config exists, overwriting...")        
		else:
			config.add_section('utils')            
		for variable, value in utilsDictionary.iteritems():
			config.set('utils', variable, value)
		surveyLogger.info("wrote utils config to "+configfilename+".cfg.")

	if(spectralDictionary):
		if(config.has_section('spectralSearch')):
			surveyLogger.info("spectralSearch config exists, overwriting...")        
		else:
			config.add_section('spectralSearch')            
		for variable, value in spectralDictionary.iteritems():
			config.set('spectralSearch', variable, value)
		surveyLogger.info("wrote spectralSearch config to "+configfilename+".cfg.")

	if(spatialDictionary):
		if(config.has_section('spatialSearch')):
			surveyLogger.info("spatialSearch config exists, overwriting...")        
		else:
			config.add_section('spatialSearch')            
		for variable, value in spatialDictionary.iteritems():
			config.set('spatialSearch', variable, value)
		surveyLogger.info("wrote spatialSearch config to "+configfilename+".cfg.")

	with open(configfilename+'.cfg', 'wb') as configfile:
		config.write(configfile)

def readConfig(surveyLogger,configfilename):
	"""
	Returns all of the needed information from the config file
	called <name>.cfg.  Also checks to make sure all of the 
	config parameters are set based on the configure dictionaries
	given in the configDictionaryList.
	"""
	surveyDictionary = {}
	mosaicDictionary = {}
	utilsDictionary = {}
	spectralDictionary = {}
	spatialDictionary = {}

	try:
		checkForFiles(surveyLogger,[configfilename+".cfg"])
		surveyLogger.info('Reading from config file ('+configfilename+'.cfg)')            
		config = ConfigParser.RawConfigParser()
		config.read(configfilename+'.cfg')
        
		if(config.has_section('survey')):
			surveyDictionary = dict(config.items('survey'))
		
		if(config.has_section('mosaic')):
			mosaicDictionary = dict(config.items('mosaic'))

		if(config.has_section('utils')):
			utilsDictionary = dict(config.items('utils'))

		if(config.has_section('spectralSearch')):
			spectralDictionary = dict(config.items('spectralSearch'))
	 
		if(config.has_section('spatialSearch')):
			spatialDictionary = dict(config.items('spatialSearch'))
	
		return surveyDictionary,mosaicDictionary,utilsDictionary,spectralDictionary,spatialDictionary

	except(FileNotFound):
		raise FileNotFound
		return

def checkConfig(surveyLogger, referenceDictionary,testDictionary):
	"""
	Checks a dictionary against a refernece to make sure that all
	of the parameters are there.  If all is good, it'll returen the 
	checked dictionary.  If not, it'll return the reference dictionary
	and raise an exception.
	"""
	try:
		for key in referenceDictionary:
			item = testDictionary[key]
		return testDictionary

	except KeyError as inst:
		surveyLogger.critical("Cannont find "+inst.args[0]+" in the config file.")
		raise KeyError
		return referenceDictionary

def Print(surveyLogger,configfile,label):
	"""
	Prints out information about the various objects to the terminal and to the log file.
	"""
	surveyLogger.info('Reading %s variables...'%label)
	i = 0
	for variable, value in configfile.iteritems():
		i += 1
		logString = "%s) %s = %s"%(i,variable,value)
		surveyLogger.info(logString)

def initLogger(name):
	"""
	Sets up and returns a properly configured logging object.
	"""
	surveyLogger = logging.getLogger(name)
	surveyLogger.setLevel(logging.DEBUG)
	# Prevents duplicate log entries after reinitialization.                                                        
	if(not surveyLogger.handlers):
		fh = logging.FileHandler('logdir/'+name+'.log')
		fh.setLevel(logging.DEBUG)
		ch = logging.StreamHandler()
		ch.setLevel(logging.DEBUG)
		FORMAT = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
		formatter = logging.Formatter(FORMAT,datefmt="%y-%m-%d %H:%M:%S")
		fh.setFormatter(formatter)
		ch.setFormatter(formatter)
		surveyLogger.addHandler(fh)
		surveyLogger.addHandler(ch)

	return surveyLogger
 
def _quotefn(filename):
	if filename is None:
		return None
	else:
		return "'" + filename + "'"
