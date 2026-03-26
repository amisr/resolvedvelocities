"""
LoggerInit.py

Written by: Steven Chen

This class will set up the logging module

CHANGELOG:

    11/21/2013 - Initial implementation for creation of logging
    02/04/2014 - Added change to parse different extensions

"""
import logging
import logging.config
import logging.handlers
import sys
import glob

"""
Logger Class:
    Class will attempt to initialize the logging module
    
"""
class LoggerInit:
    
    def __init__(self,config='config/log.ini'):
        """
            Initializes the logging module
            
            Inputs:
                config : passable string input for where to load the log config file
            
            Outputs:
                None
        """
        # Read in logging configuration file
        try:
        	
            # Read in logging configuration file
            config = glob.glob(config)
            logging.config.fileConfig(config, disable_existing_loggers=False)
            self.logger = logging.getLogger(__name__)
            self.logger.info("Logger created")
            self.logger.info("**** Using Python %s.%s.%s ****" % sys.version_info[0:3]) 
            
        except:
            # On exception will create an error log file and quit the program
            format = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
            datefmt = '%Y-%m-%d %H:%M:%S'
            logging.basicConfig(filename="log/"+self.__class__.__name__+".log",
                format=format,
                datefmt=datefmt,)
            self.logger = logging.getLogger(__name__)
            consoleHandler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter(fmt=format,datefmt=datefmt)
            consoleHandler.setFormatter(formatter)
            self.logger.addHandler(consoleHandler)
            self.logger.critical("Could not load logging... quitting now")
            sys.exit(0)