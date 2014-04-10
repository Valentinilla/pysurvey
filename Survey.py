#!/usr/bin/python

__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

from SurveyUtils import *
from Mosaic import *
from makeMosaic import *
from makeCorrection import *
from combineMosaics import *
from dsampleMosaic import *
from testModule import *
from spectralAnalysis import *
from spatialAnalysis import *
from cleanMosaic import *
from splitMosaic import *

class Survey:
	
	def __init__(self, survey='MySurvey', species='HI', mosaic='skymap', configFile = False, 
		surveyConfig = {'survey':'MySurvey','species':'HI'},
		mosaicConfig = {'mosaic':'skymap','lon':'INDEF','lat':'INDEF','z1':'INDEF','z2':'INDEF','side':'INDEF'},
		utilsConfig = {'tcmb':2.7, # Cosmic Microwave Background temperature (K)
		'tspin':150., # Excitation or Spin temperature (K)
		'xfactor':1.9e20, # CO Factor - Strong & Mattox (1996): X=NH2/Wco (K-1 cm-2 km-1 s)
		'c':1.823e18, # Costant (cm-2)
		'pc2cm':3.08567758e18, # Conversion factor from pc to cm (cm)
		'poverk':4000.,
		'p':1.0, # Fraction of HI emission originating behind the HISA cloud
		'fn':1.0}, # Fraction of particle density contributed by the HISA gas, fn = n_hisa/n_tot
		spectralConfig = {"n_spectral" : 7, #size of box for spectral smoothing
		"max_loops" : 1000,#maximum number of CLEAN iterations
		"residual_frac" : 0.03,#residual fraction of smoothed for CLEAN cutoff
		"clip_spectral" : 0.8, #fraction of r_max for adding to correction
		"gain_spectral" : 0.25, #fraction of residual height added to correction
		"fwhm_spectral" : 8, #FWHM of the Gaussian for CLEAN loop, in km/s
		"hisa_f_spectral" : -2.0, #residual amplitude factor for potential HISA in spectral search
		"temp_critical" : 30., #brightness temperature threshold
		"fit_narrow" : 2.0, #FWHM of Gaussian for narrow fit (km/s)
		"fit_broad" : 4.0,  #FWHM of Gaussian for broad fit (km/s)
		"fit_qual" : 2.0, #Gaussian fit reliability cutoff
		"dip_cut" : 0.6, #Cutoff for min morphological "dip" for spectral HISA
		"fwhm_spatial_hisa" : 5, #FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
		"min_hisa" : 2.0}, #cutoff for min HISA amplitude after spatial smoothing
		spatialConfig = {"n_spatial" : 15, #size of box for spatial smoothing
		"max_loops" : 1000, #max number of loops for CLEAN algorithm
		"high" : 10, #10th (or Mth) highest peak in residual used as rmax
		"residual_frac" : 0.03, #fraction of max smoothed height for CLEAN loop cutoff
		"clip_spatial" : 0.5, #fraction of r_max for adding to correction
		"gain_spatial" : 0.25, #fraction of residual added to correction
		"fwhm_spatial" : 20, #fwhm of gaussian for CLEAN, in arcmin
		"noise_resolve" : 20, #angular resolution for calculation of sigma_obs, in minutes
		"hisa_f_spatial" : -2., #residual amplitude factor for potential HISA in spatial search
		"temp_critical" : 30., #min unabsorbed temperature for candidate HISA
		"amp_min_first" : 4., #cutoff amplitude for first HISA filter (pre-smoothing)
		"fwhm_spatial_hisa" : 5, #FWHM of Gaussian for spatial smoothing of HISA, in units of pixels
		"min_hisa" : 2.}): #cutoff for min HISA amplitude after spatial smoothing

		surveyConfig['survey'] = survey
		surveyConfig['species'] = species
		mosaicConfig['mosaic'] = mosaic
		
		self.logger = initLogger(survey+'_'+mosaic+'_'+species+'_Analysis')
		self.configfilename = 'config/'+survey+'_'+mosaic

		if(configFile):
			try:
				surveyConfigRead,mosaicConfigRead,utilsConfigRead,spectralConfigRead,\
				spatialConfigRead = readConfig(self.logger,self.configfilename)
			except(FileNotFound):
				self.logger.critical("One or more needed files do not exist")
				return
			try:
				surveyConfig = checkConfig(self.logger,surveyConfig,surveyConfigRead)
			except(KeyError):
				return
			try:
				mosaicConfig = checkConfig(self.logger,mosaicConfig,mosaicConfigRead)
			except(KeyError):
				return
			try:
				utilsConfig = checkConfig(self.logger,utilsConfig,utilsConfigRead)
			except(KeyError):
				return
			try:
				spectralConfig = checkConfig(self.logger,spectralConfig,spectralConfigRead)
			except(KeyError):
				return
			try:
				spatialConfig = checkConfig(self.logger,spatialConfig,spatialConfigRead)
			except(KeyError):
				return

		self.surveyConf   = surveyConfig
		self.mosaicConf   = mosaicConfig
		self.utilsConf    = utilsConfig
		self.spectralConf = spectralConfig
		self.spatialConf  = spatialConfig

		self.flag_existance = False

		self.ret = re.compile('\n')
		Print(self.logger,self.surveyConf,'survey')

	def writeConfig(self):
		"""
		Writes all of the initialization variables to the config file called <surveyname>.cfg.
		"""
		writeConfig(self.logger, surveyDictionary=self.surveyConf,mosaicDictionary=self.mosaicConf,\
		utilsDictionary=self.utilsConf,spectralDictionary=self.spectralConf,spatialDictionary=self.spatialConf)


	def makeObs(self,type='brightness temperature'):
		"""
		Reads the header and gets the data
		"""
		try:				
			self.obs = Mosaic(self.surveyConf,self.mosaicConf,type)
			self.logger.info(self.ret.subn(', ',str(self.obs))[0])
				
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return

	def cleanMosaic(self,scale_data=False):
		"""
		Access mosaic's attributes through self.clean
		Input parameters: 
			species = 'HI'(default), 'CO'
			scale   = True (write scaled data in FITS file), False (do not scale data)
		"""
		try:
			self.mosaic
 		except AttributeError:
			#self.logger.critical("Obs object does not exist. Create it first with the 'makeObs' function.")
			self.logger.critical("Mosaic object does not exist. Create it first with the 'loadMosaic' function.")
			return

		try:
			self.clean = cleanMosaic(self.mosaic,scale_data)
			self.logger.info(self.ret.subn(', ',str(self.clean))[0])

		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return

	def generateMosaic(self,species='HI'):
		"""
		Generate CGPS-like mosaic. 
		Input parameters: species = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS),
		'CO' (Wco) (only for CGPS)
		Access mosaic's attributes through self.msc
		"""
		try:
			self.mosaic
			self.mosaic.species = species
 		except AttributeError:
			self.logger.critical("Obs object does not exist. Create it first with the 'makeObs' function.")
			return

		try:
			self.msc = makeMosaic(self.mosaic,self.mosaicConf)
			self.logger.info(self.ret.subn(', ',str(self.msc))[0])
			self.flag_existance = True

		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return


	#def loadMosaic(self,species='HI',type='brightness temperature'):
	def loadMosaic(self,species='HI',type='brightness temperature',datatype='original',nmsc=1,totmsc=1):
		"""
		Load a mosaic.
		Input parameters: 
			- species  = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS column density),
				    'CO' (Wco) (only for CGPS)
			- type     = 'brightness temperature'(defualt),'column density'
			- datatype = 'original' (default),'clean' (after applying data-clean methods, see cleanMosaic)
			- load     = True, False (original survey mosaic)
		Access mosaic's attributes through self.mosaic
		"""
		try:
			#if self.flag_existance:
			#	del self.msc
			#	self.logger.info("Free memory")
			
			self.mosaic = Mosaic(self.surveyConf,self.mosaicConf,type,species,datatype,nmsc,totmsc)
			#self.mosaic = Mosaic(self.surveyConf,self.mosaicConf,type,species,load=True)
			self.logger.info(self.ret.subn(', ',str(self.mosaic))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return

	def downSampleMosaic(self,scale=1):
		"""
		Access mosaic's attributes through self.lowres
		Input parameters: 
			species = 'HI'(default), 'CO'
			scale   = True (write scaled data in FITS file), False (do not scale data)
		"""
		try:
			self.mosaic
 		except AttributeError:
			#self.logger.critical("Obs object does not exist. Create it first with the 'makeObs' function.")
			self.logger.critical("Mosaic object does not exist. Create it first with the 'loadMosaic' function.")
			return

		try:
			self.lowres = dsampleMosaic(self.mosaic,scale)
			self.logger.info(self.ret.subn(', ',str(self.lowres))[0])

		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return


	# It needs to be implemented
	def spectralSearch(self):
		"""
		Access mosaic's attributes through self.spec
		"""
		try:
			self.mosaic
 		except AttributeError:
			self.logger.critical("Mosaic object "+species+" does not exist. Create it first with the 'loadMap' function.")
			return

		try:
			self.spec = spectralAnalysis(self.mosaic,self.spectralConf)
			self.logger.info(self.ret.subn(', ',str(self.spec))[0])

		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return


	def getColumnDensity(self,species='HI'):
		"""
		Calculate the column density of a mosaic.
		Input parameters: species = 'HI'(default),'HISA'(only for CGPS),'HI+HISA'(only for CGPS)
		Access mosaic's attributes through self.coldens
		"""
		try:
			self.mosaic
			self.mosaic.newspec = species
 		except AttributeError:
			self.logger.critical("Mosaic object "+species+" does not exist. Create it first with the 'loadMap' function.")
			return
		try:
			self.coldens = makeCorrection(self.mosaic,self.mosaicConf,self.utilsConf)
			self.logger.info(self.ret.subn(', ',str(self.coldens))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return
			
	def getGalpropMap(self,species='HI',type='column density'):
		try:
			self.mosaic
			self.mosaic.newspec = species
 		except AttributeError:
			self.logger.critical("Mosaic object "+species+" does not exist. Create it first with the 'loadMap' function.")
			return
		try:
			self.galprop = makeAnnuli(self.mosaic,self.mosaicConf,self.utilsConf)
			self.logger.info(self.ret.subn(', ',str(self.galprop))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return

	def combineMosaics(self,mosaic='skymap',species='HI',type='column density',dim='2D'):
		
		try:
			self.skyregion = combineMosaics(self.surveyConf,mosaic,species,type,dim)
			self.logger.info(self.ret.subn(', ',str(self.skyregion))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return
	
	def makePlot(self,plot='NH vs Ts',l=0.,b=0.):
		
		try:
			self.data = plotNvsTs(self.logger,self.obs,self.utilsConf,plot,l,b)
			self.logger.info(self.ret.subn(', ',str(self.data))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
			return
					

	def deleteMosaic(self):
		"""
		Delete the last object loaded with the loadMosaic function.
		Usage: 
			survey.loadMosaic(species='some',type='some')
			deleteMosaic()
		"""
		try:
			self.mosaic
 		except AttributeError:
			self.logger.critical("Cannot delete Mosaic object because does not exist.")
			return

		filename = self.mosaic.mosaicFile
		checkForFiles(self.logger,[filename])
		self.logger.info('Removing file: '+filename)
		os.remove(filename)

	def splitMosaic(self,num=1):
		"""
		Access mosaic's attributes through self.test
		"""
		try:
			self.mosaic
 		except AttributeError:
			self.logger.critical("Mosaic object "+species+" does not exist. Create it first with the 'loadMap' function.")
			return
		try:
			self.split = splitMosaic(self.mosaic,num)
			self.logger.info(self.ret.subn(', ',str(self.split))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")		


	def testFunction(self,species='HI'):
		"""
		Access mosaic's attributes through self.test
		"""
		try:
			self.mosaic
			self.mosaic.newspec = species
 		except AttributeError:
			self.logger.critical("Mosaic object "+species+" does not exist. Create it first with the 'loadMap' function.")
			return
		try:
			self.test = testModule(self.mosaic,self.mosaicConf,self.utilsConf)
			self.logger.info(self.ret.subn(', ',str(self.test))[0])
			
		except(FileNotFound):
			self.logger.critical("One or more needed files do not exist")
		

def printCLIHelp():
    	"""
	This function prints out the help for the CLI.
	"""
	cmd = os.path.basename(sys.argv[0])
	
	print """
				- Survey - 
	
	Perform Survey analysis for different line emissions. 
	You can use the command line functions listed below or run this module
	from within python.
	
	%s (-h|--help) ... This help text.
	
	%s (-a|--analyze) (-n |--surveyname=)<surveyname> ...  Perform an analysis
	on <surveyname>.  <surveyname> is the prefix used for this analysis.
	You must already have a configuration file if using the command
	line interface.
	
	%s (-i|--initialize) ... Generate a default config file called
	example.cfg.  Edit this file and rename it <surveyname>.cfg for use
	in the Survey module.
	
	""" %(cmd,cmd,cmd)
	
def main():
	"""
	Command-line interface.  Call this without any options for usage notes.
	"""
	import getopt
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hiamxcb:n:', ['help',
									'analyze',
									'initialize',
									'surveyname',])
		# Loop through first and check for the surveyname
		haveMosaic = False
		surveyname = 'example'
		for opt,val in opts:
			if opt in ('-n','--surveyname'):
				haveMosaic = True
				surveyname = val

		for opt, val in opts:
			if opt in ('-h', '--help'):
				printCLIHelp()
				return
			elif opt in ('-a', '--analyze'):
				if not haveMosaic: raise getopt.GetoptError("Must specify surveyname, printing help.")
				mysurvey = Survey(surveyname, True)
				print "Analysis start here!!"        
				return
			elif opt in ('-i', '--initialize'):
				print "Creating example configuration file called example.cfg"
				mysurvey = Survey(surveyname)
				mysurvey.writeConfig()
				return
                
		if not opts: raise getopt.GetoptError("Must specify an option, printing help.")
		
	except getopt.error as e:
		print "Command Line Error: " + e.msg
		printCLIHelp()


if __name__ == '__main__': 
	main()

