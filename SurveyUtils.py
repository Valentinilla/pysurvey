#!/usr/bin/env python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

import re

import os
import sys
import logging
import ConfigParser

# lmfit info: http://newville.github.com/lmfit-py/index.html
from lmfit import minimize, Parameters
from lmfit.printfuncs import*

from numpy import *

from scipy.signal import fftconvolve
from scipy.optimize import fsolve
from scipy import ndimage

import pyfits

# Avoid math module because of conflicts with 
# functions in numpy (math module operates in 
# 1D arrays, numpy's in multidimensional arrays)
# If it's unavoidable then call functions as
# math.function or numpy.function

#################################################
# GLOBAL VARIABLES
#################################################
glob_Tb  = 'brightness temperature'
glob_ITb = 'integrated brightness temperature'
glob_N   = 'column density'
glob_ncpu = 14
glob_annuli = 'Galprop' # Ackermann2012:Galprop
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
	obj_lon_max = obs.x+((obs.nx-1.)*abs(obs.dx))/2.
	obj_lon_min = obs.x-((obs.nx-1.)*abs(obs.dx))/2.
	if survey == 'LAB':
		obj_lat_max = obs.y+((obs.ny-1.)*obs.dy)
		obj_lat_min = obs.y
	else:
		obj_lat_max = obs.y+((obs.ny-1.)*obs.dy)/2.
		obj_lat_min = obs.y-((obs.ny-1.)*obs.dy)/2.
	
	lon_max = max(lon1,lon2)
	lon_min = min(lon1,lon2)
	lat_max = max(lat1,lat2)
	lat_min = min(lat1,lat2)
	
	# It needs adjustments!!
	if 0:#not survey == 'Galprop':
		if lon_max>obj_lon_max or lon_min<obj_lon_min:
			surveyLogger.critical("The longitude within the .cfg file")
			surveyLogger.critical("doesn't match that of the mosaic!")
			sys.exit(0)
		if lat_max>obj_lat_max or lat_min<obj_lat_min:
			surveyLogger.critical("The latitude within the .cfg file")
			surveyLogger.critical("doesn't match that of the mosaic!")	
			sys.exit(0)
	
	l1 = int(round(obs.px+(lon1-obs.x)/obs.dx))
	l2 = int(round(obs.px+(lon2-obs.x)/obs.dx))
	b1 = int(round(obs.py+(lat1-obs.y)/obs.dy))
	b2 = int(round(obs.py+(lat2-obs.y)/obs.dy))
	
	l = [l1,l2]
	b = [b1,b2]

	# Get the sign of delta_l and delta_b
	lsign = int((l[1]-l[0])/fabs(l[1]-l[0]))
	bsign = int((b[1]-b[0])/fabs(b[1]-b[0]))
	sign = [lsign,bsign]
	
	# Order the indexes
	if l1>l2: l = [l2,l1]
	if b1>b2: b = [b2,b1]
	
	length = rint(side/fabs(obs.dy)*2)
	
	# Make an equal sides (Lx=Ly) mosaic
	if length%2 == 0: length = length+1
	if (l[1]-l[0]) < length: l[0] = l[0]-1
	if (b[1]-b[0]) < length: b[0] = b[0]-1
	
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

def get_nth_minvalue(a, nth):
	b = a.flatten()
	res_sort = sort(b)
	c = res_sort[nth:]
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

#################################################
#	END GENERAL FUNCTIONS
#################################################

#################################################
#	START MAP CORRECTIONS
#################################################
# CO corrections
def moment_mask(surveyLogger,T,zmax,dx,dv):
	'''
	CO correction: Moment Mask method (T.M.Dame)
	'''
	# Calculate the rms for raw data
	rms_t = getRMS(surveyLogger,T,zmax)
	
	# Degrading the resolution spatially and in velocity by a factor of 2
	fwhm_s = fabs(dx)*2 #px
	fwhm_v = fabs(dv)*2
	
	sig_s = fwhm_s/sqrt(8*log(2))
	sig_v = fwhm_v/sqrt(8*log(2))
	Tsmooth = ndimage.gaussian_filter(T,sigma=(sig_v,sig_s,sig_s),order=(0,0,0))
	
	# Calculate the rms for smoothed data
	rms_ts = getRMS(surveyLogger,Tsmooth,zmax)
	
	# Set the clipping level equal 5 times the rms noise in Tsmooth
	Tclipping = 5*rms_ts
	
	# Generate a masking cube initially filled with zeros with the same dimensions as Tb
	Mask = zeros(Tsmooth.shape)
	
	# Unmask the pixel with a value > Tclipping
	index = where(Tsmooth>Tclipping)
	Mask[index] = 1
	
	# Calculate the moment-masked cube
	T[:,:,:] = Mask*T
					
	return T

def getRMS(surveyLogger,T,zmax):
	'''
	CO correction: get root mean square
	'''
	# Compute the standard deviation along the specified axis:
	# std = sqrt(mean(abs(x - x.mean())**2))
	sigma_array = zeros(3,dtype=float)
	sigma_array[0] = sqrt(mean(T[zmax-15,:,:]**2))#std(T[235,:,:])
	sigma_array[1] = sqrt(mean(T[zmax-10,:,:]**2))#std(T[240,:,:])
	sigma_array[2] = sqrt(mean(T[zmax-5,:,:]**2))#std(T[245,:,:])
	rms = amin(sigma_array)
	
	surveyLogger.info("RMS = %s"%rms)
	return rms

# HI correction
def correct_continuum2( (T,vec) ):
	'''
	HI correction: remove continuum subtraction artifacts
	'''

	nz = T.shape[0]
	ny = T.shape[1]
	nx = T.shape[2]

	glon = vec[0]
	glat = vec[1]
	vel = vec[2]
	dv = vec[3]
	mosaic = vec[4]

	path = '/afs/ifh.de/group/that/work-sf/survey/batch/sgps/hi/regions/'
	
	def selection(file,region,zvec,drvec,location='up'):

		input = open(file,'r')
		regions = input.readlines()
		
		circle,x,y,r = [],[],[],[]
		for value in regions:
			if value[0:6] == 'circle':
				if value.split('circle(')[1].split(',')[2].split(')')[1][1:2] == '#':
					circle.append(value.split('text={')[1].split('}')[0])
					xx = round(float(value.split('circle(')[1].split(',')[0]))
					x.append(xx)
					yy = round(float(value.split('circle(')[1].split(',')[1]))
					y.append(yy)
					rr = round(float(value.split('circle(')[1].split(',')[2].split(')')[0]))
					r.append(1.)#rr)

		if len(zvec) != 2*len(drvec):
			print 'Error: zvec != 2*drvec in region %s'%region
			print len(zvec), zvec
			print 2*len(drvec), drvec
			sys.exit(0)
		
		for i,c in enumerate(circle):
			if c == region:
				j=0
				for k in xrange(0,len(zvec),2):
					r[i] = drvec[j]#r[i]+drvec[j]
					k1 = zvec[k]-2
					k2 = zvec[k+1]+1
					y1 = y[i]-r[i]
					y2 = y[i]+r[i]
					x1 = x[i]-r[i]
					x2 = x[i]+r[i]
					#print k,k+1,k1,k2,j,drvec[j],len(zvec)
					T[k1:k2,y1:y2,x1:x2]=patching(T,x1,x2,y1,y2,k1,k2,location=location)
					j += 1
	
	if mosaic == 'G258.0':
		# Region 1a		
		x11a,x12a = 609,624
		y11a,y12a = 180,195
		z11a,z12a = 136,147
		# Region 1b
		x11b,x12b = 611,621
		y11b,y12b = 182,192
		z11b,z12b = 152,173
		# Region 1c
		x11c,x12c = 612,620
		y11c,y12c = 183,191
		z11c,z12c = 195,201
		# Region 2a
		x21a,x22a = 777,795
		y21a,y22a = 141,159
		z21a,z22a = 135,166
		# Region 2b
		x21b,x22b = 777,795
		y21b,y22b = 141,159
		z21b,z22b = 168,234

		T[z11a:z12a,y11a:y12a,x11a:x12a] = patching(T,x11a,x12a,y11a,y12a,z11a,z12a,location='right')
		T[z11b:z12b,y11b:y12b,x11b:x12b] = patching(T,x11b,x12b,y11b,y12b,z11b,z12b,location='right')
		T[z11c:z12c,y11c:y12c,x11c:x12c] = patching(T,x11c,x12c,y11c,y12c,z11c,z12c,location='right')
		T[z21a:z22a,y21a:y22a,x21a:x22a] = patching(T,x21a,x22a,y21a,y22a,z21a,z22a)
		T[z21b:z22b,y21b:y22b,x21b:x22b] = patching(T,x21b,x22b,y21b,y22b,z21b,z22b)
		
	if mosaic == 'G268.0':
		# Region 1a
		x11a,x12a = 428,548
		y11a,y12a = 1,95
		z11a,z12a = 100,174
		# Region 1b
		x11b,x12b = 474,505
		y11b,y12b = 20,52
		z11b,z12b = 174,192
		# Region 2
		x21,x22 = 360,400
		y21,y22 = 10,50
		z21,z22 = 160,180
		
		T[z11a:z12a,y11a:y12a,x11a:x12a] = patching(T,x11a,x12a,y11a,y12a,z11a,z12a,location='right')
		T[z11b:z12b,y11b:y12b,x11b:x12b] = patching(T,x11b,x12b,y11b,y12b,z11b,z12b,location='right')
		T[z21:z22,y21:y22,x21:x22] = patching(T,x21,x22,y21,y22,z21,z22)

	if mosaic == 'G278.0':
		# Region 1
		x11,x12 = 110,140
		y11,y12 = 10,40
		z11,z12 = 98,151
		# Region 2a
		x21a,x22a = 835,850
		y21a,y22a = 20,35
		z21a,z22a = 107,115
		# Region 2b
		x21b,x22b = 819,872
		y21b,y22b = 1,55
		z21b,z22b = 115,132
		# Region 2c
		x21c,x22c = 838,849
		y21c,y22c = 24,30
		z21c,z22c = 132,135
		# Region 2d
		x21d,x22d = 827,864
		y21d,y22d = 10,44
		z21d,z22d = 138,150
		# Region 2e
		x21e,x22e = 815,880
		y21e,y22e = 2,68
		z21e,z22e = 150,171
		# Region 2f
		x21f,x22f = 827,864
		y21f,y22f = 10,44
		z21f,z22f = 171,177

		T[z11:z12,y11:y12,x11:x12] = patching(T,x11,x12,y11,y12,z11,z12,location='right')
		T[z21a:z22a,y21a:y22a,x21a:x22a] = patching(T,x21a,x22a,y21a,y22a,z21a,z22a)
		T[z21b:z22b,y21b:y22b,x21b:x22b] = patching(T,x21b,x22b,y21b,y22b,z21b,z22b)
		T[z21c:z22c,y21c:y22c,x21c:x22c] = patching(T,x21c,x22c,y21c,y22c,z21c,z22c)
		T[z21d:z22d,y21d:y22d,x21d:x22d] = patching(T,x21d,x22d,y21d,y22d,z21d,z22d)
		T[z21e:z22e,y21e:y22e,x21e:x22e] = patching(T,x21e,x22e,y21e,y22e,z21e,z22e)
		T[z21f:z22f,y21f:y22f,x21f:x22f] = patching(T,x21f,x22f,y21f,y22f,z21f,z22f)
	
	if mosaic == 'G288.0':
		file = path+'regions288.dat'
				
		z1a = [94,96,143,163]
		selection(file,'r1a',z1a,[10,50],location='right')
		z1b = [96,143]
		selection(file,'r1b',z1b,[40],location='up')
		z1c = [88,96,96,143]
		selection(file,'r1c',z1c,[5,60],location='up')
		z1d = [94,96,96,117,132,139]
		selection(file,'r1d',z1d,[5,20,20],location='left')
		z1e = [96,117]
		selection(file,'r1e',z1e,[28],location='right')
		
		z2a = [102,112]
		selection(file,'r2a',z2a,[5])
		z2b = [84,104,104,112,126,136]
		selection(file,'r2b',z2b,[15,8,15])
		z2c = [66,115,125,136]
		selection(file,'r2c',z2c,[15,15])
		z2d = [93,106,106,110]
		selection(file,'r2d',z2d,[20,5])
		z2e = [80,110,110,112]
		selection(file,'r2e',z2e,[25,5],location='right')
		z2f = [80,105,105,109]
		selection(file,'r2f',z2f,[15,5],location='left')
		z2g = [92,110]
		selection(file,'r2g',z2g,[18],location='left')

		z3a = [128,138,138,169,202,216,220,234,276,284]
		selection(file,'r3a',z3a,[5,18,5,5,3])

		z4a = [109,150]
		selection(file,'r4a',z4a,[35])
		
		z5a = [96,108,108,147,147,151]
		selection(file,'r5a',z5a,[10,50,5])

	if mosaic == 'G298.0':
		file = path+'regions298.dat'
		
		z1a = [104,133,155,165]
		selection(file,'r1a',z1a,[10,10],location='right')
		z2a = [106,203]
		selection(file,'r2a',z2a,[20],location='down')
		z3a = [104,140,140,207]
		selection(file,'r3a',z3a,[20,30],location='down')

	if mosaic == 'G308.0':
		file = path+'regions308.dat'
		
		z1a = [10,141,170,183]
		selection(file,'r1a',z1a,[10,18],location='right')

		z2a = [106,203]
		selection(file,'r2a',z2a,[20],location='down')

		z3a = [108,153,167,185,234,243]
		selection(file,'r3a',z3a,[10,10,10],location='down')

		z4a = [87,101,101,112,112,141,141,144,144,155,169,181]
		selection(file,'r4a',z4a,[10,20,38,30,26,40],location='down')		
		z4b = [170,180]
		selection(file,'r4b',z4b,[10],location='right')	

	if mosaic == 'G318.0':
		file = path+'regions318.dat'
		
		z1a = [119,135,144,182,184,203]
		selection(file,'r1a',z1a,[30,18,10],location='right')

		z2a = [136,161,162,184,189,203]
		selection(file,'r2a',z2a,[15,15,15],location='right')
		
		z3a = [92,140,165,174,186,254,297,306]
		selection(file,'r3a',z3a,[8,8,8,8],location='right')
		z3b = [110,136,163,177,182,201]
		selection(file,'r3b',z3b,[8,15,5],location='right')
		z3c = [97,203]
		selection(file,'r3c',z3c,[15],location='up')
		z3d = [125,137,148,176,179,188,189,198]
		selection(file,'r3d',z3d,[8,15,10,10],location='right')

		z4a = [194,198]
		selection(file,'r4a',z4a,[10],location='right')
		z4b = [103,183]
		selection(file,'r4b',z4b,[5],location='right')
		z4c = [124,174,186,198]
		selection(file,'r4c',z4c,[3,3],location='right')
		z4d = [116,206]
		selection(file,'r4d',z4d,[3],location='right')

		z5a = [126,129,129,160,160,168,168,176,176,186,186,199,199,203]
		selection(file,'r5a',z5a,[5,40,8,40,8,40,5],location='right')
		z5b = [118,154,189,202]
		selection(file,'r5b',z5b,[3,3],location='right')

		z6a = [112,179,192,200]
		selection(file,'r6a',z6a,[4,5],location='right')

	if mosaic == 'G328.0':
		file = path+'regions328.dat'
		
		z1a = [158,236]
		selection(file,'r1a',z1a,[5],location='right')
		z1b = [147,152,152,175,175,184,192,202,202,210,210,212,225,235]
		selection(file,'r1b',z1b,[13,20,15,15,24,13,13],location='up')
		z1c = [147,157]
		selection(file,'r1c',z1c,[5],location='up')
		z1d = [157,173,173,180,194,200,200,209,209,216,226,232]
		selection(file,'r1d',z1d,[30,5,15,30,5,5],location='up')
		z1e = [146,149,149,165,165,191,197,199,199,210,210,213,223,237]
		selection(file,'r1e',z1e,[5,40,5,5,40,5,5],location='up')

		z2a = [96,104,104,127,127,139,139,183,183,202,202,208,208,215,219,237]
		selection(file,'r2a',z2a,[5,40,5,40,5,40,5,10],location='right')
		z2b = [231,236]
		selection(file,'r2b',z2b,[10],location='down')
		z2c = [140,155,169,188,193,203,203,205,231,236]
		selection(file,'r2c',z2c,[5,4,5,14,5],location='down')
		z2d = [141,158,167,174,226,234]
		selection(file,'r2d',z2d,[5,10,4],location='up')
		z2e = [136,218,220,235]
		selection(file,'r2e',z2e,[10,5],location='up')
		z2f = [139,176,201,211,219,225]
		selection(file,'r2f',z2f,[4,4,4],location='down')

		z3a = [196,202]
		selection(file,'r3a',z3a,[5],location='left')
		z3b = [98,122,130,148,156,195,202,214]
		selection(file,'r3b',z3b,[4,5,6,6],location='left')
		z3c = [83,243,254,268]
		selection(file,'r3c',z3c,[7,7],location='left')
		z3d = [164,178,196,212,218,234]
		selection(file,'r3d',z3d,[5,10,10],location='left')
		z3e = [158,185,200,214]
		selection(file,'r3e',z3e,[5,6],location='right')
		z3f = [128,151,177,179]
		selection(file,'r3f',z3f,[8,8],location='left')
		z3g = [150,156,156,158,158,177,177,179,179,182,182,185,204,215,219,234]
		selection(file,'r3g',z3g,[8,40,50,40,50,40,50,50],location='right')
		z3gb = [185,188,196,203,235,239]
		selection(file,'r3g',z3gb,[20,10,20],location='up')
		z3h = [182,185]
		selection(file,'r3h',z3h,[12],location='up')

		z4a = [158,163,163,181,181,185,185,194,194,199,199,212,212,214,219,220,220,227,227,235]
		selection(file,'r4a',z4a,[14,40,5,40,20,40,8,15,35,15],location='left')
		z4b = [217,220,227,234]
		selection(file,'r4b',z4b,[8,8],location='right')

		z5a = [104,126,142,153,198,213,227,234]
		selection(file,'r5a',z5a,[8,8,8,8],location='right')
	
	if mosaic == 'G338.0':
		file = path+'regions338.dat'
		
		z1a = [150,165,175,207,214,220,229,247]
		selection(file,'r1a',z1a,[5,5,5,5],location='right')
		z1ab = [209,213,250,257]
		selection(file,'r1a',z1ab,[5,5],location='left')
		z1b = [220,223,251,254]
		selection(file,'r1b',z1b,[8,5],location='right')

		z2a = [186,212,222,234,248,253]
		selection(file,'r2a',z2a,[6,6,6],location='right')
		z2b = [212,216,219,226]
		selection(file,'r2b',z2b,[10,10],location='right')
		z2c = [204,206,206,217,217,220,221,237,237,251,251,258]
		selection(file,'r2c',z2c,[5,20,5,15,5,15],location='right')

		z3a = [165,178,178,184,221,229,243,254]
		selection(file,'r3a',z3a,[8,6,6,7],location='right')
		z3b = [87,121,124,162,174,190,215,224]
		selection(file,'r3b',z3b,[15,15,12,12],location='right')
		z3bb = [190,215,242,244,244,255]
		selection(file,'r3b',z3bb,[18,18,15],location='left')
		z3c = [95,115,122,147,150,166,179,260]
		selection(file,'r3c',z3c,[3,3,3,3],location='left')
		z3d = [92,155,156,165,174,215,219,221,221,228,234,238]
		selection(file,'r3d',z3d,[7,7,7,7,20,7],location='up')
		z3db = [240,249,251,254]
		selection(file,'r3d',z3db,[7,7],location='left')
		z3e = [172,199,199,206,206,239,245,257]
		selection(file,'r3e',z3e,[5,30,5,5],location='right')

		z4a = [92,147,147,163,185,190,190,210,212,215,215,230,230,246,249,254]
		selection(file,'r4a',z4a,[14,8,8,14,12,15,13,10],location='right')
		z4b = [90,97,97,117,117,139,139,164,191,207]
		selection(file,'r4b',z4b,[5,25,5,25,5],location='down')
		z4bb = [207,217,217,231,231,233,233,256]
		selection(file,'r4b',z4bb,[22,25,23,5],location='up')
		z4c = [87,93,98,107,107,115,127,140,140,147,196,207,207,212,214,228,228,232]
		selection(file,'r4c',z4c,[14,8,14,14,5,14,5,20,5],location='right')

		z5a = [187,195,198,200,222,227]
		selection(file,'r5a',z5a,[5,5,5],location='up')
		z5b = [182,242,246,255]
		selection(file,'r5b',z5b,[5,5],location='left')

		z5c = [150,154,154,173,173,205,205,214,214,237,237,238,243,246,246,249,249,257]
		selection(file,'r5c',z5c,[5,35,40,5,40,15,5,40,10],location='left')
		z5d = [205,214,242,246,249,255]
		selection(file,'r5d',z5d,[6,6,5],location='down')
		z5e = [164,168,168,173,173,205,207,216,216,233,243,245,245,254]
		selection(file,'r5e',z5e,[8,15,35,15,35,5,13],location='up')

	if mosaic == 'G348.0':
		file = path+'regions348.dat'
		
		z1a = [279,283,283,323]
		selection(file,'r1a',z1a,[8,29],location='right')

		z2a = [281,304]
		selection(file,'r2a',z2a,[5],location='right')
		z2b = [270,280,280,288,288,315,315,316]
		selection(file,'r2b',z2b,[5,25,45,25],location='left')
		z2c = [151,158,171,202,233,314]
		selection(file,'r2c',z2c,[8,8,8],location='right')
		z2d = [253,266,274,314]
		selection(file,'r2d',z2d,[10,10],location='right')

		z3a = [5,52,68,75,88,110,148,159,159,162,162,178,178,181,181,198,217,221,221,234,234,241,244,268,268,313]
		selection(file,'r3a',z3a,[5,5,5,5,35,40,30,5,5,36,5,5,36],location='left')
		z3b = [269,283,286,295]
		selection(file,'r3b',z3b,[5,5],location='left')
		z3c = [257,302]
		selection(file,'r3c',z3c,[5],location='left')
		z3d = [184,230,264,272,285,293]
		selection(file,'r3d',z3d,[4,4,5],location='down')
		z3db = [285,289]
		selection(file,'r3d',z3db,[21],location='up')
		z3e = [95,109,141,163,163,206,244,262,264,271,273,279,283,295,300,306]
		selection(file,'r3e',z3e,[5,5,6,5,10,10,10,6],location='up')
		z3f = [148,215,268,295,300,304]
		selection(file,'r3f',z3f,[8,8,8],location='left')
		z3g = [282,291]
		selection(file,'r3g',z3g,[5],location='left')
		z3h = [5,81,281,291,302,312]
		selection(file,'r3h',z3h,[5,10,10],location='right')
		z3i = [246,259,259,270,270,276,276,294,294,306,306,312,312,315]
		selection(file,'r3i',z3i,[5,35,10,35,10,35,8],location='up')
		z3l = [255,257,257,265,265,279,281,306,306,311,311,315]
		selection(file,'r3l',z3l,[5,35,10,5,25,5],location='right')

		z4a = [264,287,297,314]
		selection(file,'r4a',z4a,[6,6],location='down')		
		z4b = [150,164,269,287,311,313]
		selection(file,'r4b',z4b,[5,7,10],location='right')
		z4c = [242,255,255,317]
		selection(file,'r4c',z4c,[5,44],location='right')
		z4d = [205,228]
		selection(file,'r4d',z4d,[4],location='up')
		z4e = [253,276,306,313,]
		selection(file,'r4e',z4e,[5,5],location='up')

	return T

	
def patching(ar,x1,x2,y1,y2,z1,z2,location='up'):
	'''
	HI correction: remove continuum subtraction artifacts
	'''
	nx = ar.shape[2]
	ny = ar.shape[1]
	nz = ar.shape[0]
		
	region = ar[z1:z2,y1:y2,x1:x2]
	sample = array(region)#.shape)
	delx = int(fabs(x2-x1))
	dely = int(fabs(y2-y1))
	
	if location == 'up':
		sample = ar[z1:z2,y2:y2+dely,x1:x2]
		#sample  = sample[:,::-1,:]
	elif location == 'down':
		sample = ar[z1:z2,y1-dely:y1,x1:x2]
	elif location == 'left':
		sample = ar[z1:z2,y1:y2,x1-delx:x1]
	elif location == 'right':
		sample = ar[z1:z2,y1:y2,x2:x2+delx]
	
	delta = int(fabs(z1-z2))
	dx,dy = 1,1
	p1 = zeros((delta))
	p2 = zeros((delta))
	p3 = zeros((delta))
	for z in xrange(delta):
		if size(ar[z1+z,y2:y2+dy,x1:x2])==0.: 
			p1[z]=1.
		else:
			p1[z] = mean(ar[z1+z,y2:y2+dy,x1:x2])
		if size(ar[z1+z,y1:y2,x1-dx:x2])==0.:
			p2[z]=1.
		else:
			p2[z] = mean(ar[z1+z,y1:y2,x1-dx:x2])
		if size(ar[z1+z,y1:y2,x2:x2+dx])==0.:
			p3[z]=1.
		else:
			p3[z] = mean(ar[z1+z,y1:y2,x2:x2+dx])
	
	x0 = region.shape[-1]/2.
	y0 = region.shape[-2]/2.
	radius = max(x0,y0)
	output = array(region)
	
	wmaxp = zeros((delta))
	for z in xrange(delta):
		wmaxp[z] = max(p1[z],p1[z],p3[z])
	wmax = amax(wmaxp)
	
	#weights = (p1[:]+p2[:]+p3[:])/3.
	for x in xrange(0,region.shape[-1]):
		for y in xrange(0,region.shape[-2]):
			distance = sqrt((x-x0)**2 + (y-y0)**2)
			if distance <= radius:
				spec = sample[:,y,x]
				maxval = amax(spec)
				if maxval == 0.: maxval = 1.
				output[:,y,x] = (spec/maxval)*wmax
				for z in xrange(int(z2-z1)):
					if amax(output[z,y,x]) > wmaxp[z]:
						s = random.normal(wmaxp[z],10,100)
						#plotFunc(arange(100),[s])
						#exit(0)
						output[z,y,x] = s[50]
	
	for z in xrange(delta):
		output[z,:,:] = ndimage.gaussian_filter(output[z,:,:],sigma=(1,1),order=(0,0))
	for x in xrange(0,region.shape[-1]):
		for y in xrange(0,region.shape[-2]):
			distance = sqrt((x-x0)**2 + (y-y0)**2)
			if distance >= radius:
				output[:,y,x] = region[:,y,x]
	return output

# HI correction
def correct_continuum( (T,data) ):
	'''
	HI correction: remove continuum subtraction artifacts
	'''

	nz = T.shape[0]
	ny = T.shape[1]
	nx = T.shape[2]
	
	data[isnan(data)] = 0.
	
	from sklearn.mixture import GMM
	classif = GMM(n_components=2)
	classif.fit(data.reshape((data.size, 1)))
	threshold = mean(classif.means_)
	
	mask = data > threshold
	filled = ndimage.morphology.binary_fill_holes(mask)
	coded_regions, num_regions = ndimage.label(filled)
	
	data_slices = ndimage.find_objects(coded_regions)
	region = [coord for coord in data_slices]
	del data_slices
	#print classif.means_
	#print threshold
	#print num_regions
	#exit(0)
	
	for z in xrange(nz):

		slice = T[z,:,:]
		
		# loop over each region of the continuum
		for x in region:
			a,b = (max(x[0].start-1,0),min(x[0].stop+1,ny-1))
			c,d = (max(x[1].start-1,0),min(x[1].stop+1,nx-1))
			
			# Source in the continuum-data
			source = data[a:b,c:d]
			
			# Do not consider structure larger than 80 px
			if max(source.shape[0],source.shape[1]) > 130:#80:
				continue
			
			# Mean value of the line-data at the source location 
			test = mean(slice[a:b,c:d])
			
			lx = fabs(c-d)-1
			ly = fabs(a-b)-1

			x31 = c-lx
			x32 = d-lx
			x41 = c+lx
			x42 = d+lx
			
			y11 = a-ly
			y12 = b-ly
			y21 = a+ly
			y22 = b+ly

			# correct if outside the boundaries
			x31,x32 = check_boundaries(x31,x32,nx)
			x41,x42 = check_boundaries(x41,x42,nx)

			y11,y12 = check_boundaries(y11,y12,ny)
			y21,y22 = check_boundaries(y21,y22,ny)
			
			# Test-regions in the line-data around the source
			# match two arrays if zonei is smaller than source
			zone1 = slice[y11:y12,c:d]
			zone2 = slice[y21:y22,c:d]
			zone3 = slice[a:b,x31:x32]
			zone4 = slice[a:b,x41:x42]

			max_val = max(mean(zone1),mean(zone2),mean(zone3),mean(zone4))
			min_val = min(mean(zone1),mean(zone2),mean(zone3),mean(zone4))

			zone1 = match_arrays(source,zone1)
			zone2 = match_arrays(source,zone2)
			zone3 = match_arrays(source,zone3)
			zone4 = match_arrays(source,zone4)
			
			# weights depending on how many non-zero elements are in a region
			w1 = float(len(where(zone1 != 0.)[0]))/(zone1.shape[0]*zone1.shape[1])
			w2 = float(len(where(zone2 != 0.)[0]))/(zone2.shape[0]*zone2.shape[1])
			w3 = float(len(where(zone3 != 0.)[0]))/(zone3.shape[0]*zone3.shape[1])
			w4 = float(len(where(zone4 != 0.)[0]))/(zone4.shape[0]*zone4.shape[1])

			if test > (max_val) or test < (min_val):
				T[z,a:b,c:d] = (w1*zone1+w2*zone2+w3*zone3+w4*zone4)/(w1+w2+w3+w4)
			
	return T

# HI corrections
def correct_data(T,rms=4):
	'''
	HI correction: smooth negative pixels
	'''
	# Setting up thresholds
	threshold = -rms # root mean square
	
	nz = T.shape[0]
	ny = T.shape[1]
	nx = T.shape[2]
	#print nx,ny,nz

	for z in xrange(nz):	
		slice = T[z,:,:]
		
		mask = slice < threshold
		filled = ndimage.morphology.binary_fill_holes(mask)
		coded_regions, num_regions = ndimage.label(filled)
		del filled
		data_slices = ndimage.find_objects(coded_regions)
		region = [coord for coord in data_slices]
		del data_slices
		
		# loop over each region
		for x in region:
		
			a,b = (max(x[0].start-1,0),min(x[0].stop+1,ny-1))
			c,d = (max(x[1].start-1,0),min(x[1].stop+1,nx-1))
			
			area = slice[a:b,c:d]

			#T[z,a:b,c:d] = ndimage.median_filter(area,(3,3))
			fwhm_spat = 3 #px
			sig = fwhm_spat/sqrt(8*log(2))
			T[z,a:b,c:d] = ndimage.gaussian_filter(area,sigma=(sig,sig),order=(0,0))
	return T
	
def check_boundaries(x1,x2,n):
	
	# correct if outside the boundaries
	if x1 < 0: x1 = 0
	if x1 > n:
		x1 = n-1
		x2 = n
	
	if x2 > n: x2 = n
	if x2 < 0:
		x1 = 0
		x2 = 1
	
	return int(x1),int(x2)

def match_arrays(a,b):

	if b.shape != a.shape:
		dummy = zeros(a.shape)
		dummy[:b.shape[0],:b.shape[1]] += b
		return dummy
	else:
		return b


def nearest_neighbors(xarray,yarray):
	
	deltax = x - reshape(x,len(x),1) # compute all possible combination
	deltay = y - reshape(y,len(y),1) # 
	dist = sqrt(deltax**2+deltay**2)
	dist = dist + identity(len(x))*dist.max() # eliminate self matching
	# dist is the matrix of distances from one coordinate to any other
	return dist


def remove_region(cubemap,x1,x2,y1,y2,z1,z2,samples='l:r:u:d'):
	'''
	Needed to remove M31 from the LAB cubemap.
	It can also be used to remove different small regions.
	'''
	for z in xrange(z1,z2):
		
		source = cubemap[z,:,:]
		test = source[y1:y2,x1:x2]
		
		lx = fabs(x2-x1)-1
		ly = fabs(y2-y1)-1
		
		x31 = x1-lx
		x32 = x2-lx
		x41 = x1+lx
		x42 = x2+lx
				
		y11 = y1-ly
		y12 = y2-ly
		y21 = y1+ly
		y22 = y2+ly
		
		# correct if outside the boundaries
		x31,x32 = check_boundaries(x31,x32,cubemap.shape[2])
		x41,x42 = check_boundaries(x41,x42,cubemap.shape[2])
		
		y11,y12 = check_boundaries(y11,y12,cubemap.shape[1])
		y21,y22 = check_boundaries(y21,y22,cubemap.shape[1])
		
		# Test-regions in the line-data around the source
		# match two arrays if zonei is smaller than source
		# weights depending on how many non-zero elements are in a region
		if samples.find('d') >= 0:
			zone1 = source[y11:y12,x1:x2]
			zone1 = match_arrays(test,zone1)
			w1 = float(len(where(zone1 != 0.)[0]))/(zone1.shape[0]*zone1.shape[1])
		else:   zone1, w1 = 0.,0.
		if samples.find('u') >= 0:
			zone2 = source[y21:y22,x1:x2]
			zone2 = match_arrays(test,zone2)
			w2 = float(len(where(zone2 != 0.)[0]))/(zone2.shape[0]*zone2.shape[1])
		else:   zone2, w2 = 0.,0.
		if samples.find('l') >= 0:
			zone3 = source[y1:y2,x31:x32]
			zone3 = match_arrays(test,zone3)
			w3 = float(len(where(zone3 != 0.)[0]))/(zone3.shape[0]*zone3.shape[1])
		else:   zone3, w3 = 0.,0.
		if samples.find('r') >= 0:
			zone4 = source[y1:y2,x41:x42]
			zone4 = match_arrays(test,zone4)
			w4 = float(len(where(zone4 != 0.)[0]))/(zone4.shape[0]*zone4.shape[1])
		else:   zone4, w4 = 0.,0.
		
		cubemap[z,y1:y2,x1:x2] = (w1*zone1+w2*zone2+w3*zone3+w4*zone4)/(w1+w2+w3+w4)
		
	return cubemap[z1:z2,y1:y2,x1:x2]

def get_intervalEdges(xarray,lowerEdge,upperEdge):
	'''
	xarray = [180,-180]
	'''
	if lowerEdge > 180.: lowerEdge = lowerEdge-360
	if upperEdge > 180.: upperEdge = upperEdge-360
	subarray1 = where(xarray <= upperEdge)
	subarray2 = where(xarray >= lowerEdge)
	subarray = [x for x in subarray1[0] if x in subarray2[0]]
	#print 'Longitude:',subarray,xarray[subarray]
	#print min(subarray),max(subarray)

	return min(subarray),max(subarray),subarray # x1:x2

def moment_mask2(T,dx,dv):
	'''
	LAB Large/Small Magellanic Cloud correction: based on Moment Mask method (T.M.Dame)
	'''
	# Degrading the resolution spatially and in velocity by a factor of 2
	fwhm_s = fabs(dx)*2 #px
	fwhm_v = fabs(dv)*2
	
	sig_s = fwhm_s/sqrt(8*log(2))
	sig_v = fwhm_v/sqrt(8*log(2))
	Tsmooth = ndimage.gaussian_filter(T,sigma=(sig_v,sig_s,sig_s),order=(0,0,0))
	
	# Calculate the rms for smoothed data
	rms_ts = 0.02 # LAB
	
	# Set the clipping level equal 5 times the rms noise in Tsmooth
	Tclipping = 5*rms_ts
	
	# Generate a masking cube initially filled with zeros with the same dimensions as Tb
	Mask = zeros(Tsmooth.shape)
	
	# Unmask the pixel with a value < Tclipping
	index = where(Tsmooth<Tclipping)
	Mask[index] = 1
	
	# Calculate the moment-masked cube
	T[:,:,:] = Mask*T
					
	return T

#################################################
#	END MAP CORRECTIONS
#################################################

#################################################
#	START ROTATION CURVE
#################################################
def getAnnuli(ver):
	
	annuli = 0
	rmin = []
	rmax = []
	
	if ver == 'Ackermann2012':
		# Boundaries of Galactocentric Annuli - M.Ackermann et al
		# (The Astrophysical Journal, 750:3-38, 2012)
		annuli = 17
		rmin = [0.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.]
		rmax = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.5,7.,8.,10.,11.5,16.5,19.,50.]
		#ann_boundaries = [ [i,rmin[i],rmax[i]] for i in xrange(annuli)]

	if ver == 'Galprop':
		annuli = 9
		rmin = [0.,1.767249,3.524893,5.520308,7.505054,9.5,11.5,13.5,15.5]
		rmax = [1.767249,3.524893,5.520308,7.505054,9.5,11.5,13.5,15.5,50]

	return rmin,rmax,annuli

def getVelPec(glon,glat,paper="SBD2010"):
	''' 
	Calculate peculiar velocity of the sun based on one of the two papers
	- "DB1998":  Dehnen & Binney, 1998, MNRAS, 298, 387
	- "SBD2010": Schoenrich, Binney & Dehnen, 2010, MNRAS, 403, 1829
	'''
	U0 = 0. # radial component (velocity in/out from center)
	V0 = 0. # tangential component
	W0 = 0. # vertical component
	vpec = 0.

	if paper == "DB1998":
		U0 = 10.0
		V0 = 5.25
		W0 = 7.17
	elif paper == "SBD2010":
		U0 = 11.1
		V0 = 12.24
		W0 = 7.25

	vpec = U0*cos(glon)*cos(glat)+V0*sin(glon)*cos(glat)+W0*sin(glat)
	#print glon,glat,vpec
	#print U0*cos(glon)*cos(glat)
	#exit(0)

	return vpec*0.005

def RotCurveBissantz2003(parlist):
	'''
	Gas-flow simulation of the inner Galaxy using smoothed particle 
	hydrodynamics (SPH) and a realistic barred gravitational potential
	derived from the observed COBE DIRBE near-IR light distribution.
	(N.Bissantz et al, 2003)
	'''

	path = parlist[0]
	glon = parlist[1] 
	#glon = radians(145) 
	glat = parlist[2]
	#glat = radians(-65)
	r_sun = parlist[3]
	v_sun = parlist[4]
	r_proj = parlist[5]
	dv = parlist[6]
	dbin = parlist[7]
	N = parlist[8]

	#true_dis = dbin*(0.5+arange(N))
	#r_proj = true_dis*cos(glat)

	# Read in SPH model results
	# Get Bissantz's data
	file1 = path+'testvr1.dat'
	input = open(file1,'r')
	bissantz = input.read().split()
	n1 = 200
	vrs = array(bissantz).reshape(n1,n1)
	vrs = vrs.astype(float32)
	# free memory
	del bissantz
	#print vrs[100,98]  #-79.1474
	#print vrs[100,99]  #-56.3561
	#print vrs[100,100] #-25.6225
	
	R = 0.5
	xa = -10.+0.1*(R+arange(n1))
	ya = -10.+0.1*(R+arange(n1))
	rb = zeros((n1,n1),dtype=float32)
	rb = [sqrt(xa[i]**2+ya**2) for i in xrange(n1)]
	rb = array(rb)
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

	# Bissantz & Gerhard velocities
	xdir = ex_tan*cos(glon)+ex_norm*sin(glon)
	ydir = ey_tan*cos(glon)+ey_norm*sin(glon)
	xha = 100+floor(10*(x_sun_sph+r_proj*xdir))
	yha = 99-floor(10*(y_sun_sph+r_proj*ydir))
	ix = where((xha >= 0.) & (xha <= 199.) & (yha >= 0.) & (yha <= 199.))
	cx = size(ix[0])
	xha = xha.astype(int)
	yha = yha.astype(int)

	# Read in rotation curve
	file2 = path+'rotcurv4.dat'
	input = open(file2,'r')
	rotation = input.readlines() #4700 lines
	rotc = zeros((len(rotation)),dtype=float32)
	drot = zeros((len(rotation)),dtype=float32)
	i=0
	for value in rotation:
		ha = float(value.split('\n')[0].split()[0])
		rotc[i] = float(value.split('\n')[0].split()[1])
		drot[i] = float(value.split('\n')[0].split()[2])
		i=i+1

	# Calculate peculiar velocity of the sun: Equation (6)	
	vpec = getVelPec(glon,glat)	
	
	vbgr = zeros(N,dtype=float32)
	if(cx > 0):
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
			for k in xrange(cnt):
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
	#rad = 0.01*(0.5+arange(4700))
	#plotFunc(rad,[rotc])
	
	# Define distance from the GC (kpc)
	radi = sqrt(r_sun**2+r_proj**2-2*r_sun*r_proj*cos(glon)) # radi.shape = (760,)
	
	# Radial velocity (express rotc as a function of radi)
	x = floor(100*radi).astype(int)
	c = 100.*radi-x
	rot_curve = rotc[x]+c*(rotc[x+1]-rotc[x])
	#plotFunc(radi,[rot_curve])
	
	# Uncorrected radial velocity: Equation (5)
	v_lsr = (r_sun*rot_curve/radi-v_sun)*sin(glon)*cos(glat)
	#plotFunc(radi,[v_lsr])
	
	# Extrapolate vbgr
	pmean = r_sun*cos(glon)
	if(cos(glon) > 0.1):
		i = int(round(40.*pmean-5.))
		vmean = mean(vbgr[i-5:i+1])
		#vmean = vbgr[i]
		vrmean = mean(v_lsr[i-5:i+1])
		#vrmean = v_lsr[i]
		#vbgr[i:] = v_lsr[i:]+(vmean-vrmean)
		vbgr[i+1:559] = v_lsr[i+1:559]+(vmean-vrmean)
	
	# Merging - transition from inner and outer Galaxy
	Rt = 9. # kpc
	wradi = where(radi > Rt, 0., (Rt-radi)/2.)
	veff = v_lsr+(vbgr-v_lsr)*where(wradi>1.,1.,wradi)		
	# Corrected, effective velocity: Equation (7)
	fwhm = 3
	sigma = fwhm/sqrt(8*log(2))
	veff = ndimage.gaussian_filter(veff,sigma=sigma,order=0)-vpec
	veff[-3:] = veff[-4]

	#plotFunc(radi,[veff,v_lsr,vbgr],['Veff','Vlsr','Vbgr'],position='upper right')
	#exit(0)

	# Weights from delta veff
	dveff = array(veff)
	dveff[-1] = fabs(veff[-2]-veff[-1])
	dveff[0:N-1] = [fabs(veff[i+1]-veff[i]) for i in xrange(N-1)]
					
	# Equation (14)
	weight_veff = zeros(veff.shape)
	weight_veff = where((dveff+1.e-8)>dv,dv,dveff+1.e-8)

	return veff,dveff,weight_veff

def RotCurveClemens1985(parlist):
	'''
	Massachusetts-Stony Brook Galactic plane CO survey: The Galactic disk rotation curve
	(D.P.Clemens, The Astrophysical Journal, 295:422-436, 1985)
	'''

	path = parlist[0]
	glon = parlist[1] 
	#glon = radians(0) 
	glat = parlist[2]
	#glat = radians(0)
	r_sun = parlist[3]
	v_sun = parlist[4]
	r_proj = parlist[5]
	dv = parlist[6]
	dbin = parlist[7]
	N = parlist[8]
	
	#true_dis = dbin*(0.5+arange(N))
	#r_proj = true_dis*cos(glat)

	# Coefficients for Rsun = 8.5 kpc and Vsun = 220 kms-1
	A = [0.,3069.81,-15809.8,43980.1,-68287.3,54904.,-17731.]
	B = [325.0912,-248.1467,231.87099,-110.73531,25.073006,-2.110625]
	C = [-2342.6564,2507.60391,-1024.068760,224.562732,-28.4080026,2.0697271,-0.08050808,0.00129348]
	D = 234.88
	
	# Conditions
	r1 = 0.09*r_sun
	r2 = 0.45*r_sun
	r3 = 1.6*r_sun

	# Fill rotation curve
	rad = 0.01*(0.5+arange(4700))
	rotc = zeros(rad.shape,dtype=float32)
	for k in xrange(size(rad)):
		if rad[k] < r1:
			for i in xrange(len(A)):
				rotc[k] += A[i]*pow(rad[k],i)
		elif rad[k] < r2:
			for i in xrange(len(B)):
				rotc[k] += B[i]*pow(rad[k],i)
		elif rad[k] < r3:
			for i in xrange(len(C)):
				rotc[k] += C[i]*pow(rad[k],i)
		else:
			rotc[k] += D
	# Plot the Clemens curve
	#plotFunc(rad,[rotc])
	
	# Define distance from GC (kpc)
	radi = sqrt(r_sun**2+r_proj**2-2*r_sun*r_proj*cos(glon))
	#print amin(radi),amax(radi),amin(rad),amax(rad)
	
	# Radial velocity (express rotc as a function of radi)
	x = floor(100.*radi).astype(int)
	c = 100.*radi-x
	rot_curve = rotc[x]+c*(rotc[x+1]-rotc[x])
	#plotFunc(radi,[rot_curve])
	
	# Uncorrected radial velocity: Equation (5)
	v_lsr = (r_sun*rot_curve/radi-v_sun)*sin(glon)*cos(glat)
	
	# Calculate peculiar velocity of the sun: Equation (6)
	vpec = getVelPec(glon,glat)
		
	# Corrected, effective velocity: Equation (7)
	fwhm = 3
	sigma = fwhm/sqrt(8*log(2))
	veff = ndimage.gaussian_filter(v_lsr,sigma=sigma,order=0)-vpec
	veff[-3:] = veff[-4]
	
	#plotFunc(radi,[veff,v_lsr],['Veff','Vlsr'],position='upper right')
	#exit(0)
	
	# Weights from delta veff
	dveff = array(veff)
	dveff[-1] = fabs(veff[-2]-veff[-1])
	dveff[0:N-1] = [fabs(veff[i+1]-veff[i]) for i in xrange(N-1)]
					
	# Equation (14)
	weight_veff = zeros(veff.shape)
	weight_veff = where((dveff+1.e-8)>dv,dv,dveff+1.e-8)

	return veff,dveff,weight_veff


def Deconvolution( (Tb,Tcont,Tunab,coord,vec) ):
	"""
	Deconvolution technique - M.Pohl, P.Englmaier, and N.Bissantz
	(The Astrophysical Journal, 677:283-291, 2008)
	Limits for Bissantz2003 rotation curve: galactic longitude < |165 deg|
	"""
	#[path_curve,self.survey,self.mosaic,self.species,lat,vel,mosaic.dy,dv,utilsConf,rmin,rmax,rotcurve,maxis]
	path = vec[0]
	survey = vec[1]
	mosaic = vec[2]
	species = vec[3]
	vel = vec[5]
	dy = vec[6]
	dv = vec[7]
	C = float(vec[8]['c'])
	Ts = float(vec[8]['tspin'])
	Tcmb = float(vec[8]['tcmb'])
	rmin = vec[9]
	rmax = vec[10]
	rotcurve = vec[11]
	maxis = vec[12]

	annuli = len(rmin)		

	if maxis == 1:
		lon = vec[4]
		lat = coord
	elif maxis == 2:
		lon = coord
		lat = vec[4]

	nlon = Tb.shape[2]
	nlat = Tb.shape[1]
	nvel = Tb.shape[0]
	
	# free memory
	#del vec

	# Array to store results
	cubemap = zeros((annuli,nlat,nlon),dtype=float32)
	cosb = cos(radians(lat))

	# Line properties
	sigma = 0. 		# velocoty dispersion (standard deviation of the distribution)
	if species == 'HI' or species == 'HI_unabsorbed':
		sigma = 4.	#hi velocity dispersion (from LAB) [km s-1]
	elif species == 'CO':
		sigma = 3.	#co velocity dispersion [km s-1]
	elif species == 'HISA':
		sigma = 2.	#hisa velocity dispersion [km s-1]

	sigma_gas = sigma
	sigma_gas_inner = 5.	#velocity dispersion inner galaxy [km s-1]

	# Define dummy line profile and its weight
	L = 20 			    # half-lenght covered by the line (num of v-channels)
	ivzero = int(floor(L/dv))   # center of the line profile (index)
	vzero = ivzero*dv	    # center of the line profile (km s-1)
	iv_vec = arange(2*ivzero+1) # vector of velocity channels (indexes) of the line profile
	v_vec = iv_vec*dv	    # vector of velocity channels (km s-1) of the line profile

	#gaussian = p[0]*exp(-0.5*power( (x-p[1])/p[2],2))
	
	# Define dummy line profile and its weight
	line = gaussian(v_vec,[1,vzero,sigma_gas],normalized=False)
	sigma_line = sigma_gas*sqrt(2.*pi)
	lim = line/sigma_line

	# Line profile for GC region
	line_inner = gaussian(v_vec,[1,vzero,sigma_gas_inner],normalized=False)
	sigma_line_inner = sigma_gas_inner*sqrt(2.*pi)
	
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
	
	# Physical variables
	if rotcurve == 'Bissantz2003':
		r_sun = 8.	#kpc
		v_sun = 210. 	#km s-1
		gmax = 165. 	#+165, -165(=345) deg
		sol = 8 	#kinematically best-fitting location
	elif rotcurve == 'Clemens1985':
		r_sun = 8.5	#kpc
		v_sun = 220.	#km s-1
		gmax = 180. 	#+180, -180(=360) deg
		sol = 8		#kinematically best-fitting location
	#gmax = 1e3
	
	z_sun = 0.015		#kpc
	gal_radius = 20. 	#kpc
	gal_thick = 1. 		#kpc
	dbin = 1/gal_radius 	#0.05 kpc
	r_scale = 10. 		#radial scale [kpc]
	
	# Cuts
	v_offset = 10.		#velocity offset [10 km/s]
	lon_inner = 20.		#inner Galaxy longitude (|l|<=20) [deg]
	residual_line = 0.5	#threshold of residual line spectrum [K km/s]
	if species == 'HISA':
		residual_line = 0.01
	amp_frac = 0.2		#percentage of the peak value [x100 %]
	
	# Array definition
	N = 760
	true_dis = dbin*(0.5+arange(N))	# kpc
	
	report = open("report_"+species+".dat","w")
	
	for l in xrange(nlon):
		glo_deg = lon[l]
		#print "[%i/%i] longitude: %.3f deg"%(l,nlon,lon[l])
		if (abs(glo_deg) <= gmax):
			glon = radians(glo_deg)
			dismin = floor(r_sun*fabs(cos(glon))/dbin)
			radmin = floor(r_sun*fabs(sin(glon))/dbin)
			#r0 = 12+radmin
			#r1 = dismin-r0
			#r2 = dismin+r0
			r0 = radmin
			#r1 = r0-4
			#r2 = r0+4
			#print glon
			#print "px =",dismin,"py =",radmin
			#print r0,r1,r2
			#print r0*dbin,r1*dbin,r2*dbin
			#exit(0)
			for b in xrange(nlat):
				#print "  %i) latitude: %.3f"%(b,lat[b])
				gla_deg = lat[b]
				glat = radians(gla_deg)
			
				# If 0-pixel value in the continuum map then HISA cannot be extracted	
				if species == 'HISA' and Tcont[b,l] == 0.:
					continue
				
				# Line spectrum
				spec = array(nvel)
				spec = Tb[:,b,l]
				spec[0] = 0.
				spec[-1] = 0.
				
				# If the spectrum is empty skip this latitude (mostly for HISA)
				if not spec.any() != 0.:
					continue
				
				if species == 'HISA':
					# no need to smooth because the spectrum is quite clean (few channels)
					rspec = spec
				else:
					rspec = ndimage.gaussian_filter(spec,sigma=sigma_gas,order=0)
			
				# Find the zero-point of the spectrum
				if 1:
					zero_avg = 0.
					idx_zero_avg = where(rspec<0)
					
					if size(idx_zero_avg[0]) > 0:
						zero_avg = mean(rspec[idx_zero_avg])
					
					spec = spec-zero_avg		
					rspec = rspec-zero_avg
				
				# Define intervals and heights: Equations (4)   
				z = z_sun+true_dis*sin(glat)
				r_proj = true_dis*cos(glat)
				# distance from the GC
				radi = sqrt(r_sun**2+r_proj**2-2*r_sun*r_proj*cos(glon))
				
				# Get the rotation curve from a model (Bissantz, Clemens)
				parlist = [path,glon,glat,r_sun,v_sun,r_proj,dv,dbin,N]
				veff = zeros(radi.shape)
				dveff = zeros(veff.shape)
				weight_veff = zeros(veff.shape)
				if rotcurve == 'Bissantz2003':
					veff,dveff,weight_veff = RotCurveBissantz2003(parlist)
				elif rotcurve == 'Clemens1985':
					veff,dveff,weight_veff = RotCurveClemens1985(parlist)
				#plotFunc(r_proj,[veff,veff2],['Bissantz 2003','Clemens 1985'],position='lower right')
				#exit(0)

				wco = fabs(dv*sum(spec))
				wcb = wco/sigma_line
				wco_previous = 0
				cnt1 = 0
				
				# Start deconvolution
				while wco > residual_line:
					
					ivpeak = argmax(rspec)
					vpeak = vel[ivpeak]
					
					if species == 'HISA' and Tcont[b,l]>Tunab[ivpeak,b,l]:
						break
					
					ivlow = ivpeak-ivzero
					ivhigh = ivpeak+ivzero+1
					
					const = 1#sqrt(2.*pi)
					if ivlow >= 0:
						iv1 = 0
						if ivhigh < nvel:
							iv2 = size(iv_vec)
							sigma_line = sigma_gas*sqrt(2.*pi)
							sigma_line_inner = sigma_gas_inner*sqrt(2.*pi)
						else:
							iv2 = size(iv_vec)+(nvel-ivhigh)
							ivhigh = nvel+1
							sigma_line = fabs(dv*sum(line[iv1:iv2]))*const
							sigma_line_inner = fabs(dv*sum(line_inner[iv1:iv2]))*const
					else:
						iv1 = -ivlow
						iv2 = size(iv_vec)
						ivlow = 0
						sigma_line = fabs(dv*sum(line[iv1:iv2]))*const
						sigma_line_inner = fabs(dv*sum(line_inner[iv1:iv2]))*const
					
					# Find a match between gas velocity and rotation curve
					ivgood = where((vpeak >= veff) & (vpeak <= (veff+dveff)))
					cnt_ivgood = size(ivgood[0])
					
					#linevpeak = ones(veff.shape)*vpeak
					#plotFunc(r_proj[0:30],[veff[0:30],veff[0:30]+dveff[0:30],linevpeak[0:30]])
					#exit(0)
					
					# Standard selection of locations
					# -------------------------------
					# Gas with forbidden velocity is placed in distance bins with the best matching velocity
					# except toward the inner Galaxy (|l| < 20 deg) where for a velocity offset of more than
					# 10 km/s to the nearest allowed velocity we accept only distance bins in the GC region.
					#delta_v_list = fabs(veff-vpeak)
					#if((amin(delta_v_list) > v_offset) and (fabs(glo_deg) < lon_inner)):
					#	igalactic_rad = r1+argsort(delta_v_list[r1:r2+1])
					
					delta_v_list = fabs(veff-vpeak)
					if((amin(delta_v_list) > v_offset)):# and (fabs(glo_deg) < lon_inner)):
						#igalactic_rad = r1+argsort(delta_v_list[r1:r2+1])
						#igalactic_rad = r0+argsort(delta_v_list[r1:r2+1])
						#if r1 > r2: r1 = r2-10
						#if r1 < 0.: r1 = 0
						if r0 < 4: r0 = 4
						#igalactic_rad = argsort(delta_v_list[r0-4:r0+5])#r2+10])
						igalactic_rad = arange(r0-4,r0+5)
						#velo1 = delta_v_list[r0-4:r0+5]
						#velo2 = delta_v_list[index2[0]:index2[-1]+1]

						#print igalactic_rad
						#print index2
						#print velo1
						#print velo2
						#exit(0)
						#if len(igalactic_rad) == 0:
						#	print r0,r1,r2
					else:
						igalactic_rad = argsort(delta_v_list)

				
					#delta_v_list = veff-vpeak
					#cnt_delta = size(where(delta_v_list<0)[0])
					
					#veff_sign = veff/fabs(veff)
					#vpeak_sign = vpeak/fabs(vpeak)
					#diff_sign = where(veff_sign != vpeak_sign)
					#cnt_diff = size(diff_sign[0])

					#if cnt_diff > 0:
						#print 'here2'
						#print veff_sign,vpeak_sign
					#	ir_sun = floor(r_sun/dbin)
					#	igalactic_rad = argsort(delta_v_list[ir_sun-4:ir_sun+4])
					#elif cnt_delta > 0:
					#	igalactic_rad = argsort(delta_v_list[])					


					# The line signal is distributed among n=sol solutions with weights
					if(cnt_ivgood == 0):
						roots = sol
						ilocation = zeros(roots,dtype=float32)
						ilocation[0:roots] = igalactic_rad[0:roots]
					else:
						roots = cnt_ivgood+sol
						ilocation = zeros(roots,dtype=float32)
						ilocation[0:cnt_ivgood] = ivgood[0][0:cnt_ivgood]
						ilocation[cnt_ivgood:roots] = igalactic_rad[0:sol]
					
					# Product of three weights (velocity,kinematic,height)
					wa = zeros(roots,dtype=float32)
					for i,j in enumerate(ilocation):
						# Thickness of the gas layer - equation (15)
						#sigma_z = 1.204*((0.06-0.04*radi[j]/r_sun)+0.095*(radi[j]/r_sun)**2) # pc
						if species == 'HI' or species == 'HI_unabsorbed' or species == 'HISA':
							sigma_z = (100.+30.*radi[j]/r_sun+90.*(radi[j]/r_sun)**2)*1e-3 # kpc
						elif species == 'CO':
							sigma_z = (50-50*radi[j]/r_sun+90*(radi[j]/r_sun)**2)*1e-3 # kpc
						
						zc = 0.							
						# Warp in the outer region of the Galaxy (r >= 11 kpc)
						if(radi[j] >= 11.):
							sphi = r_proj[j]*sin(glon)/radi[j]
							cphi = (r_proj[j]*cos(glon)-r_sun)/radi[j]
							nrx = floor(10.*(radi[j]-10.))
							sphia = sphi*cos(phia[nrx])-cphi*sin(phia[nrx])
							sphib = 2.*sphi*cphi*cos(phia[nrx])+(sphi**2-cphi**2)*sin(phia[nrx])
							# equation (16)
							zc = (warpa[nrx]+warpb[nrx]*sphia+warpc[nrx]*sphib)*1e-3  # kpc
						
						# Weights from height above plane
						weight_z = gaussian(z[j],[1,zc,sigma_z],normalized=False)
						weight_z = where(weight_z>exp(-20),weight_z,exp(-20))
						# Minimize the kinematically allowed but physically unlikely placing of gas
						weight_k = gaussian(radi[j],[1,0,r_scale],normalized=False)
						
						wa[i] = weight_veff[j]*weight_k*weight_z
						
					wtot = sum(wa)
					
					amp = amp_frac*rspec[ivpeak]
					amp = where(wcb>amp,amp,wcb)
					
					# Add intensity line (= sigma_line*amp) to cubemap
					w2 = 0.
					for i,j in enumerate(ilocation):
						sigma = sigma_line
						if(radi[j] < 1.):
							w2 += wa[i]/wtot
							sigma = sigma_line_inner
						wga = wa[i]/wtot
						
						wamp = 0.
						deltaV = sigma#(sigma/sqrt(2*pi))*sqrt(8*log(2))# sigma*1.2
						if species == 'HI' or species == 'HI_unabsorbed':
							if amp > Ts-5: amp = Ts-5
							NHI = log((Ts-Tcmb)/(Ts-Tcmb-amp))*Ts*deltaV*C
							wamp = wga*NHI
							if fabs(z[j]) > 1.: radi[j] = r_sun
						elif species == 'HISA':
							Tc = Tcont[b,l]
							Tu = Tunab[ivpeak,b,l]
							NHISA = get_ampHISA(rspec[ivpeak],Tc,Tu,dy,dv,r_proj[j],sigma,vec[8])*deltaV*C
							if NHISA == None: break
							wamp =  wga*amp_frac*NHISA
							if fabs(z[j]) > 1.: radi[j] = r_sun
						elif species == 'CO':
							wamp = wga*amp*sigma
							if fabs(z[j]) > 0.2: radi[j] = r_sun
						
						for a in xrange(annuli):
							if(radi[j] >= rmin[a]) and (radi[j] < rmax[a]):
								cubemap[a,b,l] += wamp
							if (radi[j] > 50.):
								print "Distance > 50. kpc! (%.2f)"%radi[j]
					
					# Subtract the results from the original spectrum
					w1 = 1.-w2
					line1 = line[iv1:iv2]/sigma_line
					line2 = line_inner[iv1:iv2]/sigma_line_inner
					spec[ivlow:ivhigh] -= (w1*line1+w2*line2)*amp*sigma_line

					if species == 'HISA':
						rspec = spec
					else:
						rspec = ndimage.gaussian_filter(spec,sigma=sigma_gas,order=0)
					
					wco_previous = wco
					wco = fabs(dv*sum(spec))
					wcb = wco/sigma_line	
					if wco > wco_previous:
						wco = residual_line
						wcb = wco/sigma_line

					cnt1 += 1
					if cnt1 > 600:
						string = "\ncnt1 = %i\n"%(cnt1)
						string += "[glo,glat] = [%.4f,%.4f] - [l,b] = [%i,%i]\n"%(glo_deg,gla_deg,l,b)
						string += "1) wco_previous = %.3f\n"%(wco_previous)
						string += "2) wco = %.3f\n"%(wco)
						report.write(string)
					
					#print "wco = %.3f, wco_previous = %.3f"%(wco,wco_previous)
					#if fabs(vpeak) > 250.:
					#plotFunc(vel,[spec,rspec],['observed','gaussian filter'], position='upper right')
					#exit(0)
				#zcor=total(densi(370:379,gl,gb))
				#znorm=total(densi(0:378,gl,gb))
				#densi(370:379,gl,gb)=0.
				#densi(0:378,gl,gb)=densi(0:378,gl,gb)*(1.+zcor/znorm)

				#if fabs(lat[b]) < 20.:
				#	pa = sum(cubemap[-1,b,l])
				#	pb = sum(cubemap[:-1,b,l])
				#	cubemap[-1,b,l] = 0.
				#	cubemap[:-1,b,l] *= (1.+pa/pb)
				
				#print cnt1,wco
				#plotFunc(vel,[spec,rspec],['observed','gaussian filter'], position='upper right')
				#exit(0)
	report.close()

	#radius_max = r_sun*sin(radians(10.))
	#for l in xrange(nlon):
	#	if fabs(lon[l]) < 10.:
	#		glon = radians(lon[l])
	#		for b in xrange(nlat):
	#			glat = radians(lat[b])
	#			r_proj = true_dis*cos(glat)
	#			radi = sqrt(r_sun**2+r_proj**2-2*r_sun*r_proj*cos(glon))			
				#cubemap[j,:,l] += cubemap[j,:,l]
	#			j = where(radi < radius_max)			
	#			print j,size(j)
	#print radi[j]
	#exit(0)

	# Corrections for latitude
	for b in xrange(nlat):
		cubemap[:,b,:] *= cosb[b]
	
	return cubemap

#################################################
# END ROTATION CURVE
#################################################

def get_ampHISA(Tpeak,Tc,Tu,dy,dv,r,sigma,utilsConf):
	
	debug = 0	

	deltaT = -Tpeak # Tpeak = |deltaT|, since deltaT < 0 ==> deltaT = -Tpeak
	Tcmb = float(utilsConf['tcmb'])
	C = float(utilsConf['c']) # Costant (cm-2)		
	if debug: print "> deltaT = %f, Tc = %f, Tu = %f"%(deltaT,Tc,Tu)	

	# Define params
	amplitude = 0.
	theta = radians(dy) #rad
	ds = r*tan(theta)*1e3 #pc
	
	A1 = float(utilsConf['pc2cm'])*float(utilsConf['poverk'])
	A2 = float(utilsConf['fn'])*ds/(sigma*1.2*C)#(sigma*sqrt(8*log(2))*C)
	A = A1*A2
	B = Tc+float(utilsConf['p'])*Tu
	
	# Starting estimate
	par_tau = fabs(1.05*A/(0.9*B)**2)
	par_tspin = fabs(0.9*B + deltaT/(1-exp(-1.05*A)))
	if debug: print "- A = %f, B = %f"%(A,B)
	if debug: print "- partau = %f, partspin = %f"%(par_tau,par_tspin)
	
	init = [par_tau,par_tspin]
	def equations(var,A,B,deltaT):
		'''
		Ts and tau functions - S.J.Gibson et al
		(The Astrophysical Journal, 540:852-862, 2000)
		'''
		tau,temp = var
		
		# Equation (6)
		Ftemp = 1-(deltaT/(temp-B))
		if Ftemp <= 0.: Ftemp = 1.#exp(-tau)
		
		# Equation (9)
		if tau <= 0.: tau = 0.01
		Ftau = A/tau

		f1 = -log(Ftemp)-tau
		f2 = sqrt(Ftau)-temp
		
		#print tau,temp
		return array([f1, f2], dtype=float32)

	#from scipy.optimize import fsolve
	(tau,Ts),infodict,ier,mesg = fsolve(equations,init,args=(A,B,deltaT),full_output=1)
	
	# Boundary conditions
	if Ts < Tcmb+5.: # 7-10 K due to UV radiation and CR heating
		Ts = Tcmb+5.
		Ftemp = 1-deltaT/(Ts-B)
		if Ftemp <= 0.:
			tau = par_tau
		else:
			tau = -log(Ftemp)
	TsMax = Tc+(deltaT+Tu)+(float(utilsConf['p'])-1.)*Tu
	if Ts > TsMax: # For Ts = TsMax, tau --> +oo
		Ts = TsMax
		tau = par_tau #100.#-log(1-deltaT/(Ts-B))
	if tau < 0.:
		Ts = par_tspin #0.
		tau = par_tau #0.#-log(1-deltaT/B)	
	
	if debug: print "- tau = %f, Ts = %s"%(tau,Ts)
	solution = False
	if ier == 1:
		solution = True
		amplitude = tau*Ts
	if not solution:
		#print "Could not find a valid solution:"
		#print mesg
		amplitude = None#par_tau*par_tspin
	if debug: print "- amplitude = %f"%amplitude
	return fabs(amplitude)

#################################################
# START COMBINE 
#################################################
def concatenateMosaics( (list,vec) ):
	# Concatenate mosaics along the longitude axis
	dim = vec[0]
	dx = vec[1]
	msc_x = vec[2]
	
	index = 0
	c = []
	#for i in xrange(10):
	#	print i+1,list[i].header['NAXIS']
	#exit(0)
	for current, next in zip(list, list[1:]):
		d1 = current.data.shape
		d2 = next.data.shape
		if d1 != d2:
			print "Mosaic dimensions don't agree: current = %s, next = %s"%(str(d1),str(d2))
			sys.exit(0)
		if index == 0:
			if dim == '2D':
				c = concatenate((current.data[:,0:-dx], next.data[:,dx:]), axis=1)
			if dim == '3D':
				c = concatenate((current.data[:,:,0:-dx], next.data[:,:,dx:]), axis=2)
		elif index > 0 and index < len(list):
			if dim == '2D':
				c = concatenate((c[:,:-dx], next.data[:,dx:]), axis=1)
			if dim == '3D':
				c = concatenate((c[:,:,:-dx], next.data[:,:,dx:]), axis=2)
		elif index == len(list):
			if dim == '2D':
				c = concatenate((c, next.data[:,dx:msc_x]), axis=1)
			if dim == '3D':
				c = concatenate((c, next.data[:,:,dx:msc_x]), axis=2)
		index += 1
	
	return array(c)

#################################################
# END COMBINE 
#################################################
#################################################
# START SPATIAL SEARCH FUNCTIONS
#################################################
def movingaverage2D(array, w):
	window = ones((w,w))/float(w**2)
	return fftconvolve(array, window, 'same')

def spatialAverage2D(array,w):
	wmin = int(w/2)
	wmax = array.shape[0]-wmin

	a2 = movingaverage2D(array, w)

	a3 = zeros(array.shape)
	a3[wmin:wmax,wmin:wmax]=a2[wmin:wmax,wmin:wmax]

	a3[:wmin,:]=array[:wmin,:]
	a3[:,:wmin]=array[:,:wmin]
	a3[wmax:,:]=array[wmax:,:]
	a3[:,wmax:]=array[:,wmax:]

	return a3

def rms_estimation2D(signal,noise,wsize):

	rms_noise = signal - noise
	sigma = zeros(signal.shape,dtype=float32)
	
	noise_avg = ndimage.median_filter(rms_noise,wsize)
	noise_sqr_avg = ndimage.median_filter(rms_noise**2,wsize)
	sigma = sqrt(noise_sqr_avg-noise_avg**2)

	return sigma

	
	species = vec[0]
	vel = vec[2]
	dv = vec[3]
	path = vec[4]
	C = vec[5]
	Ts = vec[6]
	rmin = vec[7]
	rmax = vec[8]
	rotcurve = vec[9]
	maxis = vec[10]
	
	annuli = len(rmin)

	nlon = T.shape[2]
	nlat = T.shape[1]
	nvel = T.shape[0]

	if maxis == 1:
		lon = vec[1]
		lat = coord
	elif maxis == 2:
		lon = coord
		lat = vec[1]


def spatialSearch( (T,vec) ):
	#[obs.zmin,obs.zmax,obs.dx,spatialConf,maxis]
	kmin = vec[0]
	kmax = vec[1]
	dx = vec[2] # deg
	params = vec[3]
		
	nx = T.shape[2]
	ny = T.shape[1]
	nz = T.shape[0]
	
	N_SPATIAL         = int(params['n_spatial'])
	MAX_LOOPS         = int(params['max_loops'])
	HIGH              = int(params['high'])
	RESIDUAL_FRAC     = float(params['residual_frac'])
	CLIP_SPATIAL      = float(params['clip_spatial'])
	GAIN_SPATIAL      = float(params['gain_spatial'])
	FWHM_SPATIAL      = float(params['fwhm_spatial']) # 20 arcmin = 20'
	NOISE_RESOLVE     = float(params['noise_resolve'])
	HISA_F_SPATIAL    = float(params['hisa_f_spatial'])
	TEMP_CRITICAL     = float(params['temp_critical'])
	AMP_MIN_FIRST     = float(params['amp_min_first'])
	FWHM_SPATIAL_HISA = float(params['fwhm_spatial_hisa'])
	MIN_HISA          = float(params['min_hisa'])
	
	result = zeros((nz,ny,nx),dtype=float32)

	# spatial gaussian
	fwhm_spat = FWHM_SPATIAL/(fabs(dx)*60.) # [fwhm] = px (1 arcmin = 1/60 deg)
	sigma_spat = fwhm_spat/sqrt(8*log(2)) # [sigma_spat] = px
	L = 1.5*fwhm_spat
	xx,yy = mgrid[-L:L,-L:L]
	gauss_spat = 1/(2*pi*sigma_spat**2)*exp(-0.5*(xx**2+yy**2)/(sigma_spat**2) )
	gauss_spat *= 1/(gauss_spat).sum()
	#print FWHM_SPATIAL,fwhm_spat,FWHM_SPATIAL_HISA,dx
	#print gauss_spat.sum()
	#plotFunc2D(xx,yy,gauss_spat)
	#exit(0)

	# HISA gaussian
	sigma_spat_HISA = FWHM_SPATIAL_HISA/sqrt(8*log(2)) # FWHM = 5px	
	L = 1.5*FWHM_SPATIAL_HISA
	xx,yy = mgrid[-L:L,-L:L]
	gauss_spat_HISA = 1/(2*pi*sigma_spat_HISA**2)*exp(-0.5*(xx**2+yy**2)/(sigma_spat_HISA**2) )
	gauss_spat_HISA *= 1/(gauss_spat_HISA).sum()
	#plotFunc2D(xx,yy,gauss_spat_HISA)
	#exit(0)

	# Spatial search algorithm
	for k in xrange(kmin,kmax):

		#print "k = %s"%k
			
		# O(i,j) is the observed spectrum
		observed = T[k,:,:]
		# S(i,j) is a smoothed version of O(i,j)
		smoothed = ndimage.median_filter(observed,N_SPATIAL) #15px
		# U(i,j) is the unabsorbed spectrum (initially set to zero)
		unabsorbed = zeros(observed.shape,dtype=float32)
		# R(k) is the residual spectrum (initially set equal to S(k))
		residual = smoothed
						
		# N(i,j) is the smoothed noise on a large angular scale [20' = 66.7 px]
		noise_box = NOISE_RESOLVE/(60.*fabs(dx)) #67px
		noise = ndimage.gaussian_filter(observed,sigma=(noise_box,noise_box),order=(0,0))
			
		# Estimatation of the rms noise in O(i,j)
		sigma_obs = rms_estimation2D(observed,noise,noise_box)
			
		# Estimation of the average rms noise in S(i,j)	
		sigma_sm = zeros(observed.shape,dtype=float32)
		sigma_sm1 = rms_estimation2D(smoothed,noise,noise_box)
			
		imin = int(noise_box/2.)
		imax_lat = int(ny-(imin+1))
		imax_lon = int(nx-(imin+1))
			
		for y in xrange(imin,imax_lat):
			for x in range(imin,imax_lon):
				min1 = min(sigma_sm1[y-imin,x-imin],sigma_sm1[y-imin,x+imin])
				min2 = min(sigma_sm1[y+imin,x-imin],sigma_sm1[y+imin,x+imin])
				sigma_sm[y,x] = min(min1,min2)
			
		den = (nx-noise_box)**2
		sigma_sm_avg = (sigma_sm[imin:ny-(imin+1),imin:nx-(imin+1)]).sum()/den
			
		# Clean loop
		s_max = amax(smoothed)
		for loop in xrange(MAX_LOOPS):

			# Find the 10th highest residual
			rmax = get_nth_maxvalue(residual,HIGH)
					
			# 1. 
			if(rmax < RESIDUAL_FRAC*s_max or rmax == 0):
				break
				
			# 2.
			correction = zeros(observed.shape,dtype=float32)
			index = where(residual > CLIP_SPATIAL*rmax) 
			correction[index] = GAIN_SPATIAL*residual[index]
				
			# 3.
			#unabsorbed += fftconvolve(correction,gauss_spat,"same")
			unabsorbed += ndimage.gaussian_filter(correction,sigma=(sigma_spat,sigma_spat),order=(0,0))
			# 4.
			# calculate the rms of positive values of residual
			residual = smoothed - unabsorbed
				
			x_pos = residual[residual > 0.]
			if( len(x_pos) > 0):
				sigma_pos = sqrt(mean( power(x_pos,2) ))
			else:
				sigma_pos = 0.
				
			if(sigma_pos < sigma_sm_avg):
				break
									
			if loop == -2:
				exit(0)

		# Process results to look for potential HISA
		observed_dif = observed-unabsorbed

		# eliminate pixels with no valid noise calculations from consideration
		omit = [arange(noise_box),nx-arange(noise_box)-1]
		omit2 = [arange(noise_box),ny-arange(noise_box)-1]
		residual[:,omit] = 0.
		residual[omit2,:] = 0.
		observed_dif[:,omit] = 0.
		observed_dif[omit2,:] = 0.

		res_check = where(residual < HISA_F_SPATIAL*sigma_sm)
		obs_check = where(observed_dif < HISA_F_SPATIAL*sigma_obs)
			
		hisa_detected = zeros(observed.shape,dtype=float32)
		if(len(res_check) > 0):
			hisa_detected[res_check] = unabsorbed[res_check]-observed[res_check]
		if(len(obs_check) > 0):
			hisa_detected[obs_check] = unabsorbed[obs_check]-observed[obs_check]
			
		# perform initial filtering of detected HISA
		filter1 = where(unabsorbed < TEMP_CRITICAL)
		filter2 = where(hisa_detected < AMP_MIN_FIRST)
		hisa_detected[filter1] = 0.
		hisa_detected[filter2] = 0.
			
		# Do spatial smoothing of HISA
		#result[k,:,:] = fftconvolve(hisa_detected,gauss_spat_HISA,"same")
		result[k,:,:] = ndimage.gaussian_filter(hisa_detected,sigma=(sigma_spat_HISA,sigma_spat_HISA),order=(0,0))

	#self.logger.info("Only take smoothed HISA > 2K...")
	## Turn results into a map with 1 where smoothed HISA > 2K and 0 everywhere else
	result[result < MIN_HISA] = 0.
	
	return result

#################################################
# END SPATIAL SEARCH FUNCTIONS
#################################################
#################################################
# START SPECTRAL SEARCH FUNCTIONS
#################################################
def gaussian(x,p,normalized=True):
	# Returns a 1D gauss for convolution
	# p[0] = norm, p[1] = mu, p[2] = sigma
	if(normalized is True):
		A = 1/(p[2]*sqrt(2*pi))
	if(normalized is not True):
		A = p[0]
	g = A*exp(-0.5*power( (x-p[1])/p[2],2))
	return g

def residualG(pars, x, data=None):
	amp = pars['amp'].value
	mu = pars['mu'].value
	sigma = pars['sigma'].value
	#p0 = pars['p0'].value
	#p1 = pars['p1'].value

	if abs(sigma) < 1.e-10:
		sigma = sign(sigma)*1.e-10

	model = amp*exp(-0.5*power((x-mu)/sigma,2))#+p0+p1*x

	if data is None:
		return model

	return (data-model)

def dipFilter(spectrum,x0,dx):
	# Check for "dip" in O(k) at kmin (=x0)
	hdx = int(dx/2)

	i11 = x0-dx#+1
	i12 = x0-hdx
	i21 = x0+hdx
	i22 = x0+dx#+1

	y0 = spectrum[x0]
	y1 = sum(spectrum[i11:i12])/hdx
	y2 = sum(spectrum[i21:i22])/hdx

	Diff = 2*max([abs(y0-y1),abs(y0-y2)])
	Frac = 1 - abs(y1-y2)/Diff
							
	if(y1<y0 and y2<y0):
		Frac = 0

	if(False):
		print "Left:   y1 = %s - O[%s:%s]"%(y1,i11,i12)
		print "Center: y0 = %s - O[%s]"%(y0,x0)
		print "Right:  y2 = %s - O[%s:%s]\n"%(y2,i21,i22)
		print "Diff:  %s"%Diff
		print "Frac:  %s"%Frac

	return Frac

def spectralSearch( (T,vec) ):
	#[obs.zmin,obs.zmax,obs.dx,dv,vel,spectralConf]
	kmin = vec[0]
	kmax = vec[1]
	dx = vec[2] # deg
	dz = vec[3]
	zarray = vec[4]
	params = vec[5]
		
	nx = T.shape[2]
	ny = T.shape[1]
	nz = T.shape[0]

	N_SPECTRAL        = int(params['n_spectral'])
	MAX_LOOPS         = int(params['max_loops'])
	RESIDUAL_FRAC     = float(params['residual_frac'])
	CLIP_SPECTRAL     = float(params['clip_spectral'])
	GAIN_SPECTRAL     = float(params['gain_spectral'])
	FWHM_SPECTRAL     = float(params['fwhm_spectral'])
	HISA_F_SPECTRAL   = float(params['hisa_f_spectral'])
	TEMP_CRITICAL     = float(params['temp_critical'])
	FIT_NARROW        = float(params['fit_narrow'])
	FIT_BROAD         = float(params['fit_broad'])
	FIT_QUAL          = float(params['fit_qual'])
	DIP_CUT           = float(params['dip_cut'])
	FWHM_SPATIAL_HISA = float(params['fwhm_spatial_hisa'])
	MIN_HISA          = float(params['min_hisa'])
	
	# Array to store results
	result = zeros((nz,ny,nx),dtype=float32)

	# Set up gaussians to be used for convolution
	# spectral convolution
	sigma_spec = FWHM_SPECTRAL/sqrt(8*log(2))
	L = 3*FWHM_SPECTRAL
	deltax = fabs(dz)
	n_half = floor( (L/2)/deltax )
	xi = linspace(-n_half*deltax,n_half*deltax,num=1+2*n_half)
	gauss_spec = gaussian(xi,[0,0,sigma_spec],normalized=True)
	gauss_spec *= 1/(gauss_spec).sum()
	#plotFunc(xi,[gauss_spec])
	
	# spatial convolution (HISA gaussian)
	sigma_HISA = FWHM_SPATIAL_HISA/sqrt(8*log(2))
	L = 3*FWHM_SPATIAL_HISA #3*5
	xx,yy = mgrid[-L:L,-L:L]
	gauss_HISA = 1/(2*pi*sigma_HISA**2)*exp(-0.5*(xx**2+yy**2)/(sigma_HISA**2) )
	gauss_HISA *= 1/(gauss_HISA).sum()
	#plotFunc2D(xx,yy,gauss_HISA)
	
	# Find sigma values for gaussian fits, in units of pixels
	sigma_narrow = FIT_NARROW/(sqrt(8*log(2))*fabs(dz))
	sigma_broad = FIT_BROAD/(sqrt(8*log(2))*fabs(dz))

	# Start spectral search algorithm
	# needed for rms (see below)
	a=round((nz-kmin)/3)+kmin-1
	b=round((nz-kmin)/3)+kmin
	c=round(2*(nz-kmin)/3)+kmin-1
	d=round(2*(nz-kmin)/3)+kmin
	#print a,b,c,d #a=101, b=102, c=186, d=187
	
	for i in xrange(7,nx-9):#7
		for j in xrange(7,ny-9):
			
			# O(k) is the observed spectrum
			observed = T[:,j,i]
			# S(k) is a smoothed version of O(k)
			smoothed = ndimage.median_filter(observed,N_SPECTRAL)
			# U(k) is the unabsorbed spectrum (initially set to zero)
			unabsorbed = zeros(observed.shape,dtype=float32)
			# R(k) is the residual spectrum (initially set equal to S(k))
			residual = smoothed
			# Spectral rms noise in O(k)
			noise = observed - smoothed
			sigma_array = zeros(3,dtype=float32)
			# Compute the standard deviation along the specified axis:
			# std = sqrt(mean(abs(x - x.mean())**2))
			sigma_array[0] = std(noise[kmin:a])
			sigma_array[1] = std(noise[b:c])
			sigma_array[2] = std(noise[d:])
			sigma_obs = amin(sigma_array)
			
			# Initialize the HISA container
			HISA_merged = zeros(nz,dtype=float32)
			
			# Initialize the rms of positive values of residual
			sigma_pos = 0
			
			# Clean loop
			smax = amax(smoothed)
			for loop in xrange(MAX_LOOPS):
				rmax = amax(residual)
				# 1. 
				if(rmax < RESIDUAL_FRAC*smax or rmax == 0.):
					break
				# 2.
				correction = zeros(nz,dtype=float32)
				correction = where(residual>CLIP_SPECTRAL*rmax,residual*GAIN_SPECTRAL,correction)
				# 3.
				#unabsorbed += fftconvolve(correction,gauss_spec,"same")
				unabsorbed += ndimage.gaussian_filter(correction,sigma=sigma_spec,order=0)
				# 4.
				# calculate the rms of positive values of residual
				residual = smoothed-unabsorbed
				x_pos = residual[residual > 0.]
				
				#plotFunc(zarray,[observed,smoothed,unabsorbed,residual])
				if(len(x_pos) > 0):
					sigma_pos = sqrt(mean( power(x_pos,2) ))
				else:
					sigma_pos = 0.
				
				if(sigma_pos < sigma_obs):
					break
				
			# Process results to look for potential HISA
			if(rmax > 0):
				# output (array([190, 191, 207]),)
				suspected_hisa_index = where(residual < sigma_pos*HISA_F_SPECTRAL)
				# counting starts from right in the graph (18 initial 0-events)
				cnt1 = size(suspected_hisa_index)
			else: 
				cnt1 = 0
			
			# Skip HISA search if HICA present
			if(amin(observed) < -6.*sigma_obs):
				cnt1 = 0	
			
			counter = 0
			segments=[]
			
			# Build consecutive HISA candidates into segments
			#while(counter < cnt1):
			#print suspected_hisa_index
			#for iseg,ires in enumerate(suspected_hisa_index[0]):
			#	segments.append(residual[ires]) 
				#print iseg,ires
				
			while(counter < cnt1):
				segments = [residual[suspected_hisa_index[0][counter]]]
				ikmin = [suspected_hisa_index[0][counter],0]
				#print "> counter %i"%counter,ikmin,segments
				
				count_start = counter
				counter += 1
				while(counter < cnt1 and (suspected_hisa_index[0][counter]-1) == suspected_hisa_index[0][counter-1]):
					segments.append(residual[suspected_hisa_index[0][counter]])
					if(segments[counter-count_start] < segments[ikmin[1]] ):
						ikmin = [suspected_hisa_index[0][counter],counter-count_start]
					counter += 1

				#exit(0)
				# ikmin[0][0] = index in residual, ikmin[0][1] = index in segments, ikmin[1] = size of segment
				ikmin = [ikmin, counter-count_start]
				#print "counter %i"%counter,ikmin,segments
	
				# Fit Gaussians to candidate segments
				#perform initial filtering on candidate segments
				if(smoothed[ikmin[0][0]] > 2.*HISA_F_SPECTRAL*sigma_pos and unabsorbed[ikmin[0][0]] > TEMP_CRITICAL):
					# Pad array (segments) with 0s to make sure fit possible
					num_of_zeros = 6
					data = zeros(ikmin[1]+num_of_zeros)
					i_inf = 0.5*num_of_zeros
					i_sup = ikmin[1] + 0.5*num_of_zeros
					data[i_inf:i_sup] = fabs(array(segments))
				
					amp = amax(data)
					mu = ikmin[0][0]
					par = [amp,mu,sigma_narrow,sigma_broad]

					x = xrange( ikmin[1] + num_of_zeros ) + ( ikmin[0][0] - (ikmin[0][1]+ 0.5*num_of_zeros) )

					parN = Parameters()
					parB = Parameters()

					if((data.size-num_of_zeros) > 3): #3 (at least 3 points to fit)
						parN.add('amp', value=par[0], vary=True)
						parN.add('mu', value=par[1], vary=False)
						parN.add('sigma', value=par[2], vary=False)
						parB.add('amp', value=par[0], vary=True)
						parB.add('mu', value=par[1], vary=False)
						parB.add('sigma', value=par[3], vary=False)
					else:
						parN.add('amp', value=par[0], vary=False)
						parN.add('mu', value=par[1], vary=True)
						parN.add('sigma', value=par[2], vary=True)
						parB.add('amp', value=par[0], vary=False)
						parB.add('mu', value=par[1], vary=True)
						parB.add('sigma', value=par[3], vary=True)

					fitN = minimize(residualG, parN, args = (x, data), engine = "leastsq")
					fitB = minimize(residualG, parB, args = (x, data), engine = "leastsq")

					bestfit_parN = [float(parN['amp'].value),float(parN['mu'].value),float(parN['sigma'].value)]
					bestfit_parB = [float(parB['amp'].value),float(parB['mu'].value),float(parB['sigma'].value)]
						
					# Check quality of fit
					sigma_fitN = sqrt((fitN.residual**2).sum()/data.size)
					sigma_fitB = sqrt((fitB.residual**2).sum()/data.size)

					condition = min([bestfit_parN[0]/max([sigma_fitN,sigma_pos]),bestfit_parB[0]/max([sigma_fitB,sigma_pos])])
					if( (condition > FIT_QUAL) and fitN.success and fitB.success):
						fit_quality = True
					else:
						fit_quality = False
						
					# Process and combine fitted Gaussians
					if(fit_quality):
						HISA_narrow = gaussian(arange(nz),bestfit_parN,normalized=False)
						HISA_broad  = gaussian(arange(nz),bestfit_parB,normalized=False)

						# Only take values > 5% of peak (cut the gaussian tails)
						filterN = 0.05 * HISA_narrow[ikmin[0][0]]
						filterB = 0.05 * HISA_broad[ikmin[0][0]]

						HISA_narrow[HISA_narrow < filterN] = 0
						HISA_broad[HISA_broad < filterB] = 0

						HISA_merged_temp = where(HISA_narrow > HISA_broad, HISA_narrow, HISA_broad)

						# Check for "dip" in O(k) at ikmin
						dip_frac = dipFilter(observed,ikmin[0][0],num_of_zeros)
						if(dip_frac > DIP_CUT):
							HISA_merged += HISA_merged_temp
						
						#ar = -1*ones(observed.shape,dtype=float32)
						#trh1 = ones(nz)*(sigma_pos*HISA_F_SPECTRAL-10)
						#trh2 = 30*ones(nz)
						#f1 = [observed,smoothed,unabsorbed,10*ar+residual]
						#f2 = [60*ar+HISA_narrow,90*ar+HISA_broad,120*ar+HISA_merged]
						#f3 = [trh1,trh2]
						# U(k),S(k),R(k),HISA_N,HISA_B,HISA_M
						#plotFunc(zarray,f1[1:]+f2+f3,lbl='no label',position='lower right')
						#label=['smoothed','unabsorbed','residual']
						#plotFunc(zarray,f1[1:]+f3,lbl=label,position='upper right')
						#exit(0)
		
			# Store HISA in result_array
			result[:,j,i] = HISA_merged
			
		#print "Done with (i,j) = (%i,%i)"%(i,j)
		
	# Do spatial smoothing of HISA
	print "Spatial smoothing (convolution)..."
	for k in xrange(nz):
		result[k,:,:] = fftconvolve(result[k,:,:],gauss_HISA,"same")

	#self.logger.info("Only take smoothed HISA > 2K...")
	# Turn results into a map with 1 where smoothed HISA > 2K and 0 everywhere else
	result[result < MIN_HISA] = 0.

	return result

#################################################
# END SPECTRAL SEARCH FUNCTIONS
#################################################
#################################################
# START PLOTTING FUNCTIONS
#################################################
def plotFunc(x,func,lbl=None,position=None):
	'''
	Usage: plotFunc(x, [func1,...,func8])
	       plotFunc(x, [func1,...,func8], [lbl1,...lbl8])
	       plotFunc(x, [func1,...,func8], [lbl1,...lbl8], position='lower right')
	       plotFunc(x, [func1,...,func8], lbl='no label')
	'''
	import matplotlib.pyplot as plt
	n = len(func)
	if lbl==None:
		lbl = []
		for i in xrange(n):
			lbl.append('func%i'%i)
	if n==1:
		plt.plot(x,func[0])
	if n==2:
		plt.plot(x,func[0],x,func[1])
	if n==3:
		plt.plot(x,func[0],x,func[1],x,func[2])
	if n==4:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3])
	if n==5:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4])
	if n==6:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4],x,func[5])
	if n==7:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4],x,func[5],x,func[6])
	if n==8:
		plt.plot(x,func[0],x,func[1],x,func[2],x,func[3],x,func[4],x,func[5],x,func[6],'--',x,func[7],'--')
	if n>8:
		print "plotFunc: Max number of functions allowed is 8, %i found!"%n
		exit(0)
	if lbl!='no label':
		if position==None: position='upper left'
		plt.legend( (lbl),loc='%s'%position,shadow=False,fancybox=True)#color='black',lw=1
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

def histT():
	data  = Tb[0,50,:,:].flatten()
	data1 = Tb[0,100,:,:].flatten()
	data2 = Tb[0,150,:,:].flatten()
	data3 = Tb[0,200,:,:].flatten()
				
	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	x = data2#[data,data1,data2,data3]
	colors = ['black']#,'blue','green','red']
	# the histogram of the data
	n, bins, patches = plt.hist(x,500,normed=1,color=colors,alpha=0.75)
				
	# add a 'best fit' line
	#y = data
	#l = plt.plot(bins, y, 'r--', linewidth=1)
				
	plt.xlabel('T$_{b}$')
	plt.ylabel('Counts')
	plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
	#plt.axis([amin(data), amax(data), 0, 0.1])
	plt.grid(True)
	plt.show()


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

def setNaN2Zero(array):
	array[isnan(array)] = 0.
	return array

def getPath2(surveyLogger, key=list, mode="DESY"):

	survey = list[0]
	spec = list[1]
	flag = list[2]

	if mode == "DESY":
		disk1 = "/afs/ifh.de/group/that"
		disk2 = "/lustre/fs4/group/that/sf/Surveys"

		data1 = disk1+"/data/Radio/skymaps"
		result1 = disk1+"/data/HISA"

		workdir = disk1+"/work-sf/survey"

	elif mode == "HOME":
		disk1 = "/Volumes/My Book/DESY"
		disk2 = disk1+"/analysis/survey/results"

		data1 = disk1+"/data/Radio"
		result1 = data1+"/hisa"

		workdir = disk1+"/analysis/survey"

	elif mode == "BATCH":
		disk1 = "."#"$TMPDIR"
		disk2 = disk1+"/results"

		data1 = disk1+"/data"
		result1 = disk1+"/results"

		workdir = disk1

	path = False
	# Rot.Curve
	if key == 'rotcurve_mpohl':
		path = True
		return workdir+'/rotcurve/'

	# Mosaic list
	if key == 'list_mosaic':
		path = True
		return workdir+'/lstmosaic/'
	
	# TODO
	#
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


def getPath(surveyLogger, key="cgps_hi"):

	mode = "DESY"
	#mode = "HOME"
	#mode = "BATCH"

	if mode == "DESY":
		disk1 = "/afs/ifh.de/group/that"
		disk2 = "/lustre/fs4/group/that/sf/Surveys"

		data1 = disk1+"/data/Radio/skymaps"
		result1 = disk1+"/data/HISA"

		workdir = disk1+"/work-sf/survey"

	elif mode == "HOME":
		disk1 = "/Volumes/My Book/DESY"
		disk2 = disk1+"/analysis/survey/results"

		data1 = disk1+"/data/Radio"
		result1 = data1+"/hisa"

		workdir = disk1+"/analysis/survey"

	elif mode == "BATCH":
		disk1 = "."#"$TMPDIR"
		disk2 = disk1+"/results"

		data1 = disk1+"/data"
		result1 = disk1+"/results"

		workdir = disk1
		

	path = False
	# Rot.Curve
	if key == 'rotcurve_mpohl':
		path = True
		return workdir+'/rotcurve/'

	# Mosaic list
	if key == 'list_mosaic':
		path = True
		return workdir+'/lstmosaic/'
	
	# TODO
	#
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
		return disk1+'/soft/opt/Galprop/FITS/'
	if key=="galprop_co":
		path = True
		return disk1+'/soft/opt/Galprop/FITS/'
	
	if key=="lustre_galprop":
		path = True
		return disk2+'/Galprop/'
	if key=="lustre_galprop_hi":
		path = True
		return disk2+'/Galprop/HI/'
	if key=="lustre_galprop_hi_column_density":
		path = True
		return disk2+'/Galprop/HI/col_den/'
	if key=="lustre_galprop_hisa":
		path = True
		return disk2+'/Galprop/HISA/'
	if key=="lustre_galprop_hisa_column_density":
		path = True
		return disk2+'/Galprop/HISA/col_den/'
	if key=="lustre_galprop_co":
		path = True
		return disk2+'/Galprop/CO/'
	if key=="lustre_galprop_co_column_density":
		path = True
		return disk2+'/Galprop/CO/col_den/'

	# Dame
	if key=="dame_co":
		path = True
		return disk1+'/data/DiffuseEmission/co/'
	if key=="lustre_dame":
		path = True
		return disk2+'/Dame/'
	if key=="lustre_dame_co_column_density":
		path = True
		return disk2+'/Dame/col_den/'
	
	# LAB
	if key=="lab_hi":
		path = True
		return data1+'/hi/LAB/'
	
	if key=="lustre_lab":
		path = True
		return disk2+'/LAB/'
	if key=="lustre_lab_hi":
		path = True
		return disk2+'/LAB/HI/'
	if key=="lustre_lab_hi_split":
		path = True
		return disk2+'/LAB/HI/split/'
	if key=="lustre_lab_hi_column_density":
		path = True
		return disk2+'/LAB/HI/col_den/'

	# CGPS
	if key=="cgps_hi":
		path = True
		return data1+'/hi/CGPS/'
	if key=="cgps_hi_continuum":
		path = True
		return data1+'/continuum/cgps/'
	if key=="cgps_hisa_dat":
		path = True
		return result1+'/cgps/results/'
	if key=="cgps_co":
		path = True
		return data1+'/co/cgps/'

	if key=="lustre_cgps":
		path = True
		return disk2+'/CGPS/'
	if key=="lustre_cgps_hi":
		path = True
		return disk2+'/CGPS/HI/'
	if key=="lustre_cgps_hi_split":
		path = True
		return disk2+'/CGPS/HI/split/'
	if key=="lustre_cgps_hi_column_density":
		path = True
		return disk2+'/CGPS/HI/col_den/'
	if key=="lustre_cgps_hi_unabsorbed":
		path = True
		return disk2+'/CGPS/HI_unabsorbed/'
	if key=="lustre_cgps_hi_unabsorbed_split":
		path = True
		return disk2+'/CGPS/HI_unabsorbed/split/'
	if key=="lustre_cgps_hi_unabsorbed_column_density":
		path = True
		return disk2+'/CGPS/HI_unabsorbed/col_den/'
	if key=="lustre_cgps_hisa":
		path = True
		return disk2+'/CGPS/HISA/'
	if key=="lustre_cgps_hisa_split":
		path = True
		return disk2+'/CGPS/HISA/split/'
	if key=="lustre_cgps_hisa_column_density":
		path = True
		return disk2+'/CGPS/HISA/col_den/'
	if key=="lustre_cgps_co":
		path = True
		return disk2+'/CGPS/CO/'
	if key=="lustre_cgps_co_split":
		path = True
		return disk2+'/CGPS/CO/split/'
	if key=="lustre_cgps_co_column_density":
		path = True
		return disk2+'/CGPS/CO/col_den/'

	if key=="lustre_cgps_hisa_spectral":
		path = True
		return disk2+'/CGPS/HISA/spectral/'
	if key=="lustre_cgps_hisa_spatial":
		path = True
		return disk2+'/CGPS/HISA/spatial/'
	
	# SGPS
	if key=="sgps_hi":
		path = True
		return data1+'/hi/SGPS/'
	if key=="sgps_hi_continuum":
		path = True
		return data1+'/continuum/sgps/'
	if key=="sgps_hisa_dat":
		path = True
		return result1+'/sgps/results/'

	if key=="lustre_sgps":
		path = True
		return disk2+'/SGPS/'
	if key=="lustre_sgps_hi":
		path = True
		return disk2+'/SGPS/HI/'
        if key=="lustre_sgps_hi_split":
                path = True
                return disk2+'/SGPS/HI/split/'
	if key=="lustre_sgps_hi_column_density":
		path = True
		return disk2+'/SGPS/HI/col_den/'
	if key=="lustre_sgps_hi_unabsorbed":
		path = True
		return disk2+'/SGPS/HI_unabsorbed/'
	if key=="lustre_sgps_hi_unabsorbed_split":
		path = True
		return disk2+'/SGPS/HI_unabsorbed/split/'
	if key=="lustre_sgps_hi_unabsorbed_column_density":
		path = True
		return disk2+'/SGPS/HI_unabsorbed/col_den/'
	if key=="lustre_sgps_hisa":
		path = True
		return disk2+'/SGPS/HISA/'
        if key=="lustre_sgps_hisa_split":
                path = True
                return disk2+'/SGPS/HISA/split/'
	if key=="lustre_sgps_hisa_column_density":
		path = True
		return disk2+'/SGPS/HISA/col_den/'

	# VGPS
	if key=="vgps_hi":
		path = True
		return data1+'/hi/VGPS/'
	if key=="vgps_hi_continuum":
		path = True
		return data1+'/continuum/vgps/'
	if key=="vgps_hisa_dat":
		path = True
		return result1+'/vgps/results/'

	if key=="lustre_vgps":
		path = True
		return disk2+'/VGPS/'
	if key=="lustre_vgps_hi":
		path = True
		return disk2+'/VGPS/HI/'
        if key=="lustre_vgps_hi_split":
                path = True
                return disk2+'/VGPS/HI/split/'
	if key=="lustre_vgps_hi_column_density":
		path = True
		return disk2+'/VGPS/HI/col_den/'
	if key=="lustre_vgps_hi_unabsorbed":
		path = True
		return disk2+'/VGPS/HI_unabsorbed/'
	if key=="lustre_vgps_hi_unabsorbed_split":
		path = True
		return disk2+'/VGPS/HI_unabsorbed/split/'
	if key=="lustre_vgps_hi_unabsorbed_column_density":
		path = True
		return disk2+'/VGPS/HI_unabsorbed/col_den/'
	if key=="lustre_vgps_hisa":
		path = True
		return disk2+'/VGPS/HISA/'
        if key=="lustre_vgps_hisa_split":
                path = True
                return disk2+'/VGPS/HISA/split/'
	if key=="lustre_vgps_hisa_column_density":
		path = True
		return disk2+'/VGPS/HISA/col_den/'
		
	if(not path):
		surveyLogger.critical("Path '%s' doesn't exist."%key)
		raise FileNotFound

def getFile(surveyLogger,survey,mosaic,species,type,datatype,nmsc,totmsc,mypath):

	# Select the file according to Survey and Mosaic
	path = ''
	flag = ''
	sur = survey.lower()
	spec = species.lower()
	nmsc = int(nmsc)
	totmsc = int(totmsc)

	#if datatype == 'original':
	#	path = getPath(surveyLogger, key=sur+'_'+spec)
	#elif datatype == 'clean':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec)
	#elif datatype == '2D_col_density':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_column_density')
	#elif datatype == '3D_col_density':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_column_density')
	#elif datatype == '3D_integrated_line':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
	#elif datatype == 'processed':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec)
	#elif datatype == 'split':
	#	path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_split')


	# GALPROP
	if survey == 'Galprop':
		if species == 'HI':
			if datatype == 'original':
				path = getPath(surveyLogger, key=sur+'_hi')
				mosaic = 'TOT'
				flag = species+'_column_density_rbands_image'
			elif datatype == '2D_col_density':
				path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
				flag = species+'_column_density'
			else:
				datatypeErrorMsg(surveyLogger,datatype,surveyEntry=survey)
		elif species == 'WCO':
			if datatype == 'original':
				path = getPath(surveyLogger, key=sur+'_co')
				mosaic = 'TOT'
				flag = species+'_rbands_image'
			elif datatype == '3D_integrated_line':
				path = getPath(surveyLogger, key='lustre_'+sur+'_co')
				flag = species+'_line_rings'
			elif datatype == '2D_col_density':
				path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
				flag = 'H2_column_density'
			else:
				datatypeErrorMsg(surveyLogger,datatype,surveyEntry=survey)
		else:
			surveyLogger.critical("Only HI and WCO species available for "+survey+" survey.")
	
	# LAB
	elif survey == 'LAB':
		if species == 'HI':
			if datatype == 'original':
				path = getPath(surveyLogger, key='lab_hi')
				mosaic = 'TOT'
				flag = species+'_line_image'
			elif datatype == 'new': #glob_Tb:
				path = getPath(surveyLogger, key='lustre_lab_hi')
				flag = species+'_line_image'
			elif datatype == 'clean':
				if not mosaic == 'TOT':
					surveyLogger.warning("Only mosaic TOT has datatype = clean,")
					surveyLogger.warning("i.e. without M31,SMC, and LMC.")
					sys.exit(0)
				path = getPath(surveyLogger, key='lustre_lab_hi')
				flag = species+'_line_image_clean'
			elif datatype == '2D_col_density':
				path = getPath(surveyLogger, key='lustre_lab_hi_column_density')
				flag = species+'_column_density'
			elif datatype == 'split':
				if nmsc==0: nmsc=1
				path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_split')
				flag = '%s_line_part_%i-%i'%(species,nmsc,totmsc)
			elif not datatype == 'generic':
				datatypeErrorMsg(surveyLogger,datatype,surveyEntry=survey)
		else:
				surveyLogger.critical("Only HI species available for "+survey+" survey.")
	# DAME
	elif survey == 'Dame':
		if species == 'WCO':
			if datatype == 'original':
				path = getPath(surveyLogger, key='dame_co')
				mosaic = 'TOT'
				flag = species+'_line_image'
			elif datatype == 'new':
				path = getPath(surveyLogger, key='lustre_dame')
				flag = species+'_line_image'
			elif datatype == '2D_col_density':
				path = getPath(surveyLogger, key='lustre_dame_co_column_density')
				flag = 'H2_column_density'
			else:
				datatypeErrorMsg(surveyLogger,datatype,surveyEntry=survey)
		else:
			surveyLogger.critical("Only WCO species available for "+survey+" survey.")
	
	else:
		if datatype == 'original':
			if not (species == 'HI' or species == 'CO'):
				surveyLogger.warning("Only HI and CO have datatype = original.")
				sys.exit(0)
			path = getPath(surveyLogger, key=sur+'_'+spec)
			flag = species+'_line_image'

		elif datatype == 'clean':
			if not (species == 'HI' or species == 'HI_unabsorbed' or species == 'HISA' or species == 'CO'):
				surveyLogger.warning("Only HI and CO have datatype = clean and")
				surveyLogger.warning("only HI_unabsorbed and HISA can load it.")
				sys.exit(0)
			if species == 'HI_unabsorbed' or species == 'HISA':
				species = 'HI'
				spec = species.lower()
			path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec)
			flag = species+'_line_image_clean'

		elif datatype == '2D_col_density':
			if species == 'HI' or species == 'HISA':
				path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_column_density')
				flag = species+'_column_density'
			elif species == 'CO':
				path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
				flag = 'H2_column_density'
			elif species == 'HI+HISA':
				path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
				flag = species+'_column_density'
			elif species == 'HI+CO':
				path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
				flag = 'HI+H2_column_density'
			elif species == 'HI+HISA+CO':
				path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
				flag = 'HI+HISA+H2_column_density'

		elif datatype == '3D_col_density':
			if not (species == 'HISA' or species == 'HI_unabsorbed'):
				surveyLogger.warning("Only HISA and HI_unabsorbed have datatype = 3D_col_density.")
				sys.exit(0)
			path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_column_density')
			flag = species+'_column_density_rings'

		elif datatype == '3D_integrated_line':
			if not species == 'CO':
				surveyLogger.warning("Only CO has datatype = 3D_integrated_line.")
				sys.exit(0)
			path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
			flag = 'WCO_line_rings'

		elif datatype == 'processed':
			if not (species == 'HISA' or species == 'HI_unabsorbed'):
				surveyLogger.warning("Only HISA and HI_unabsorbed have datatype = processed.")
				sys.exit(0)
			path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec)
			flag = species+'_line'

		elif datatype == 'lowres':
			path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec)
			flag = species+'_line_lowres'

		elif datatype == 'split':
			if nmsc==0: nmsc=1
			path = getPath(surveyLogger, key='lustre_'+sur+'_'+spec+'_split')
			flag = '%s_line_part_%i-%i'%(species,nmsc,totmsc)

	filename = path+survey+'_'+mosaic+'_'+flag+'.fits'

	if datatype == 'generic':
		filename = mypath

	return	filename,mosaic

def getFile2(surveyLogger,survey,mosaic,species,type,load):
	# Select the file according to Survey and Mosaic
	path = ''
	flag = ''
	sur = (survey).lower()

	if survey == 'Galprop':
	# If Galprop
		if species == 'HI':
			if not load:
				path = getPath(surveyLogger, key=sur+'_hi')
				mosaic = 'TOT'
				flag = species+'_column_density_rbands_image'
			else:
				if type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type)
		elif species == 'WCO':
			if not load:
				if type == glob_ITb:
					path = getPath(surveyLogger, key=sur+'_co')
					mosaic = 'TOT'
					flag = species+'_rbands_image'
				else:
					typeErrorMsg(surveyLogger,type,species)
			else:
				if type == glob_ITb:
					path = getPath(surveyLogger, key='lustre_'+sur+'_co')
					flag = species+'_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type)
		else:
			surveyLogger.critical("Only HI and WCO species available for "+survey+" survey.")

	elif survey == 'LAB':
		# If LAB
		if species == 'HI':
			if not load:
				path = getPath(surveyLogger, key='lab_hi')
				mosaic = 'TOT'
				flag = species+'_line_image'
			else:
				if type == glob_Tb:
					path = getPath(surveyLogger, key='lustre_lab')
					flag = species+'_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_lab_hi_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type)
		else:
			surveyLogger.critical("Only HI species available for "+survey+" survey.")

	elif survey == 'Dame':
		# If Dame
		if species == 'WCO':
			if not load:
				path = getPath(surveyLogger, key='dame_co')
				if type == glob_ITb:
					mosaic = 'TOT'
					flag = species+'_line_image'
			else:
				if type == glob_ITb:
					path = getPath(surveyLogger, key='lustre_dame')
					flag = species+'_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_dame_co_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type)
		else:
			surveyLogger.critical("Only WCO species available for "+survey+" survey.")

	else:
		if not load:
			if species == 'HI':
				path = getPath(surveyLogger, key=sur+'_hi')
				if type == glob_Tb:
					flag = species+'_line_image'
			elif species == 'CO':
				path = getPath(surveyLogger, key=sur+'_co')
				if type == glob_Tb:
					flag = species+'_line_image'
		else:
			if species == 'HI':
				if type == glob_Tb:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi')
					flag = species+'_unabsorbed_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
					flag = species+'_unabsorbed_column_density'
				else:
					typeErrorMsg(surveyLogger,type)

			elif species == 'HISA':                            
				if type == glob_Tb:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hisa')
					flag = species+'_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hisa_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type)

			elif species == 'HI+HISA':
				if type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
					flag = species+'_column_density'
				else:
					typeErrorMsg(surveyLogger,type,typeEntry=species)

			elif species == 'CO':                              
				if type == glob_Tb:
					path = getPath(surveyLogger, key='lustre_'+sur+'_co')
					flag = 'CO_line'
				elif type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_co_column_density')
					flag = 'H2_column_density'
				else:
					typeErrorMsg(surveyLogger,type,typeEntry=species)

			elif species == 'HI+CO':                          
				if type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
					flag = 'HI+H2_column_density'
 				else:
					typeErrorMsg(surveyLogger,type,typeEntry=species)

			elif species == 'HI+HISA+CO':                              
				if type == glob_N:
					path = getPath(surveyLogger, key='lustre_'+sur+'_hi_column_density')
					flag = 'HI+HISA+H2_column_density'
				else:
					typeErrorMsg(surveyLogger,type,typeEntry=species)
	
	return	path,flag,mosaic


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

def datatypeErrorMsg(surveyLogger,datatype,surveyEntry='CGPS'):

	list = []
	if surveyEntry == 'CGPS':
		list = ['original','clean','2D_col_density','3D_col_density','3D_integrated_line','processed','split','generic']
	elif surveyEntry == 'SGPS' or surveyEntry=='VGPS':
		list = ['original','clean','2D_col_density','3D_col_density','processed','split','generic']
	elif surveyEntry == 'LAB':
		list = ['original','2D_col_density','3D_col_density','processed','split','generic']
	else:
		list = ['original']

	surveyLogger.critical("Allowed datatypes are:")
	for i,item in enumerate(list):
		surveyLogger.critical("%i. %s"%(i+1,item))
	surveyLogger.critical("Your entry is: %s."%datatype)

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

