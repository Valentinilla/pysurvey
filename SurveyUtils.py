#!/usr/bin/env python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

import os
import sys
import logging
import math
import ConfigParser

# lmfit info: http://newville.github.com/lmfit-py/index.html
from lmfit import minimize, Parameters
from lmfit.printfuncs import*

from numpy import *
from numpy.lib.stride_tricks import as_strided
from scipy.signal import fftconvolve


#################################################
#		LAB Survey
#################################################
def getMosaicCoordinate(obs,lon,lat,side):
	"""
	LAB Survey: this function allows to generate 'mosaics' like those
	defined in CGPS, SGPS, and VGPS.
	Input parameters:
	- longitude: galactic longitude of the mosaic's center
	- latitude:  galactic latitude of the mosaic's center
	- side:      length in degrees of the mosaic's side
	Outout parameters:
	- rl1, rl2:  longitude in pixels
	- rb1, rb2:  latitude in pixels
	"""
	
	lon1 = lon+side
	lon2 = lon-side
	lat1 = lat-side
	lat2 = lat+side

	#l1 = (obs.msc_lon/2. + lon1/obs.msc_del_lon)
	#l2 = (obs.msc_lon/2. + lon2/obs.msc_del_lon)
	#b1 = (obs.msc_lat/2. + lat1/abs(obs.msc_del_lat))
	#b2 = (obs.msc_lat/2. + lat2/abs(obs.msc_del_lat))

	l1 = int(round(obs.msc_ind_lon-1.+(lon1-obs.msc_ref_lon)/obs.msc_del_lon))
	l2 = int(round(obs.msc_ind_lon-1.+(lon2-obs.msc_ref_lon)/obs.msc_del_lon))
	b1 = int(round(obs.msc_ind_lat-1.+(lat1-obs.msc_ref_lat)/obs.msc_del_lat))
	b2 = int(round(obs.msc_ind_lat-1.+(lat2-obs.msc_ref_lat)/obs.msc_del_lat))

	rl1 = rint(l1)
	rl2 = rint(l2)
	rb1 = rint(b1)
	rb2 = rint(b2)

	length = rint(side/abs(obs.msc_del_lat)*2)

	if length%2 == 0:
		length = length+1

	if (rl2-rl1) < length:
		rl1 = rl1-1

	if (rb2-rb1) < length:
		rb1 = rb1-1

	return int(rl1),int(rl2),int(rb1),int(rb2)


#################################################
#		END LAB Survey
#################################################


def get_nth_maxvalue(a, nth):
	b = a.flatten()
	res_sort = sort(b)
	c = res_sort[-nth:]
	return c[0]

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
	
	#glo=lon(gl)+1.d-6
	#gla=lat(gbo)
	
	if (abs(glo_deg) < 165.):
		glon = glo_deg*pi/180.
		glat = gla_deg*pi/180.
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
		#plotFunc(proj_dis,vr)
		#print amin(vr),amax(vr)
		
		try:
			i = 0
			diff=[]
			for vel in vr:
				#print round(vel)/round(vlsr)
				if round(vel)/round(vlsr) == 1.:
					diff.append([abs(vel-vlsr),i])
					#print vel,proj_dis[i],i
				i=i+1
			#print diff
			radius = proj_dis[min(diff)[1]]
			return radius # kpc
		except ValueError:
			surveyLogger.critical("The rotation curve doesn't contain the")
			surveyLogger.critical("vlsr value of the mosaic!!")
			surveyLogger.critical("...[vrot_min, vrot_max] = [%.2f,%.2f]"%(amin(vr),amax(vr)))
			surveyLogger.critical("...vlsr_msc = %.2f"%vlsr)
			sys.exit(0)	

#--------------------------------
# START PLOTTING FUNCTIONS
#--------------------------------
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

#--------------------------------
# END ANALYSIS AREA
#--------------------------------

class FileNotFound: pass
class CommandNotFound: pass

def getPath(surveyLogger, key="cgps_hi"):

	path = False
	# Rot.Curve
	if key=="rotcurve_mpohl":
		path = True
		return '/afs/ifh.de/group/that/work-sf/survey/rotcurve/'

	# LAB
	if key=="lab_hi":
		path = True
		return '/afs/ifh.de/group/that/data/Radio/skymaps/hi/LAB/'
	if key=="lustre_lab":
		path = True
		return '/lustre/fs4/group/that/sf/LAB/'
	if key=="lustre_lab_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/LAB/col_den/'

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
		return '/lustre/fs4/group/that/sf/CGPS/'
	if key=="lustre_cgps_hi":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/HI/'
	if key=="lustre_cgps_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/HI/col_den/'
	if key=="lustre_cgps_hisa":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/HISA/'
	if key=="lustre_cgps_hisa_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/HISA/col_den/'
	if key=="lustre_cgps_co":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/CO/'
	if key=="lustre_cgps_co_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/CGPS/CO/col_den/'

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
		return '/lustre/fs4/group/that/sf/SGPS/'
	if key=="lustre_sgps_hi":
		path = True
		return '/lustre/fs4/group/that/sf/SGPS/HI/'
	if key=="lustre_sgps_hi_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/SGPS/HI/col_den/'
	if key=="lustre_sgps_hisa":
		path = True
		return '/lustre/fs4/group/that/sf/SGPS/HISA/'
	if key=="lustre_sgps_hisa_column_density":
		path = True
		return '/lustre/fs4/group/that/sf/SGPS/HISA/col_den/'
		
	if(not path):
		surveyLogger.critical("Path '%s' doesn't exist."%key)
		raise FileNotFound


def checkForFiles(surveyLogger, fileList, existence=False):
	"""
	Checks for the existence of needed files in the list.
	"""
	for filename in fileList:
		if not os.path.exists(filename) and not existence:
			surveyLogger.critical(filename+" doesn't exist.")
			raise FileNotFound
		elif os.path.exists(filename) and existence:
			surveyLogger.critical(filename+" already exist.")
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

def initLogger(mosaic, name):
	"""
	Sets up and returns a properly configured logging object.
	"""
	surveyLogger = logging.getLogger(name)
	surveyLogger.setLevel(logging.DEBUG)
	#Prevents duplicate log entries after reinitialization.                                                        
	if(not surveyLogger.handlers):
		fh = logging.FileHandler('logdir/'+mosaic+'_'+name+'.log')
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
