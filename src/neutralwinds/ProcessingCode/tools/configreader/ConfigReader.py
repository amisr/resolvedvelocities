"""
ConfigReader.py

Written by: Steven Chen

This class will read in the configuration files.

CHANGELOG:

    11/21/2013 - Modified code, now works with both Python 2.6+ and 3.x
    12/11/2013 - Heavily modified ConfigReader code. Now acts in more of a 
                 class-like structure
    01/21/2014 - Added functionality for parsing types (dict,list,tuple)
    02/04/2014 (MJN) - Revised so that interpolation works    
    02/04/2014 - Added functionality for lists separated by commas and newline
                 Added literal string parsing with quotes as well for overriding commas and newline
                 Switched to using ExtendedInterpolation
    03/05/2014 - Fixed logger not having handlers when called as a singular instance
    05/08/2014 - Modified code, now it should work with Python 2.4
    05/12/2014 - Added functionality for return keys to be lowercase
    05/13/2014 - Fixed logging parameter for Python3.4
    10/17/2014 - Added a quick patch to enable list input for files.
    
    BUGS : 
        Sometimes the dictionaries come out as OrderedDict... although this is not really
        a big issue it will cause issues when trying to use functions to operate of the 
        config dictionary
"""

# Check for python version
import sys
if sys.version_info < (3,0):
    import ConfigParser as configparser
else:
    import configparser

# Import system-wide packages **NECESSARY**
import logging
import glob
import os
import re
    

# Import User packages
from .Struct import *

# Implemented for backward compatability
try:
    import ast
except:
    import SafeEval



"""
Functions used for backwards compatibility
"""
def all(iterable):
    """
    This function imitates the 'all' function that is found in later
    versions of python.
    """
    for element in iterable:
        if not element:
            return False
    return True


"""
ConfigReader Class:
    Class handles reading in all the configuration options into a dictionary
    for later access
"""
class ConfigReader:
    
    def __init__(self,logger=None,logger_level='DEBUG'):
        """
            Initializes the ConfigReader module
            
            Inputs:
                logger = passed logger object to use
                logger_level = level of logging to use if logger is passed
                				(ignored otherwise)
            
            Outputs:
                None
        """

        # Set up logging module
        self.logger = logger or logging.getLogger(__name__)

        # If imported logger
        if logger:
            level = logging.getLevelName(logger_level)
            try:
                self.logger.setLevel(level)
            except:
                # Fix added for Python3.4
                self.logger.setLevel(logger_level)
               
        # If logger handlers do not exist, create one
        if not self.logger.root.handlers:
            # Apply handlers
            format = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
            datefmt = '%Y-%m-%d %H:%M:%S'
            consoleHandler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter(fmt=format,datefmt=datefmt)
            consoleHandler.setFormatter(formatter)
            self.logger.addHandler(consoleHandler)
           
            level = logging.getLevelName(logger_level)
            try:
                self.logger.setLevel(level)
            except:
                # Fix added for Python3.4
                self.logger.setLevel(logger_level)
        
        self.logger.debug("ConfigReader class initialized")
        
        # Try import ast library for formating
        try:
            ast
        except:
            self.logger.warning("Could not import ast library...")
           
        # Set up 
        if sys.version_info < (3,0):
            self.Config = configparser.ConfigParser()
            self.RawConfig = configparser.RawConfigParser()
        else:
            self.Config = configparser.ConfigParser(strict=False)
            self.RawConfig = configparser.RawConfigParser(strict=False)

        # Keeps case sensitivity
        self.Config.optionxform=str
            
        # Regex pattern match
        self.pattern = re.compile('\s*[,\n]\s*')
        
        
    
    def read(self,*args,**kwargs):
        """
        This function mimics the configParser read functionality and returns
        the passed arguments back as a dictionary.
        
        If length of args is 1 it will just return a dictionary in the 
        following format:
        
        {SECTIONS: {VARIABLES:} }
        
        Otherwise:
        
        {FILE: {SECTIONS: {VARIABLES:} } }
        
        This function also accepts the asterisk keyword for accepting all files
        in a given directory.

		Inputs:
                args = files to be parsed
                kwargs = type of processing to do on the files
                		Currently only 'type' key is supported with values of dict or
						object
            
        Outputs:
                cfg = returned variable with all config file values
            
        """
        
        # Parse kwargs options - Will be extended in future revisions
        if 'key' in kwargs:
            func = {}
            func['key'] = kwargs['key']
        else:
            func = None


            
        self.cfg = {}

        # Formatted file list
        files = self.format_args(args)
        
        
        for file in files:
            
            # Get the name of the file
            name = os.path.basename(file).split('.')[0]
            
            # Read the config file
            self.Config.read(file)
            self.RawConfig.read(file)
    
            # Grab all the sections in the file
                        
            sections = {}
            
            sections.update(dict(self.Config._sections))

            
            # Loop through sections and parse using Config.get, which will handle 
            # interpolation of values such as %(path)/more_path

            for section, content in sections.items():
                content.pop("__name__", None)
                for key,value in content.items():
                    try:
                        sections[section][key] = self.format_values(self.Config.get(section, key))
                    except:
                        sections[section][key] = self.format_values(self.RawConfig.get(section, key))

            # add default (from [DEFAULT] heading) which configparser handles as a special case
            if self.Config._defaults:
                sections["DEFAULT"] = {}
                for key, value in self.Config._defaults.items():
                    sections["DEFAULT"][key] = self.format_values(self.Config._defaults[key])
            else:
                self.logger.debug("No DEFAULT section found in %s, ignoring..." % (file))
                
                
                
                
            if len(files) > 1:
                self.cfg[name] = sections
                
                # Needs to be reset
                self.Config = configparser.ConfigParser()
                # Keeps case sensitivity
                self.Config.optionxform=str
                
            else:
                self.cfg = sections
                        
        
        # Parse optional parameters
        if func:
            self.cfg = self.parse_options(self.cfg, func)

        
        # DEBUG statements
        if len(files) >1:
            for filename, contents in self.cfg.items():
                self.logger.debug("%s:", filename)
                for sections, keys in contents.items():
                    self.logger.debug("\t[%s]" %(sections))
                    for key, value in keys.items():
                        self.logger.debug("\t\t%s : %s" %(key,value))
        else:
            for sections, keys in self.cfg.items():
                self.logger.debug("\t[%s]" %(sections))
                for key, value in keys.items():
                    self.logger.debug("\t\t%s : %s" %(key,value))
        
        if 'type' not in kwargs:
            # Default will return dict
            self.logger.debug("***Returned as nested dictionary***")
            return self.cfg
        elif 'type' in kwargs:
            if kwargs['type'] == dict:
                self.logger.debug("***Returned as nested dictionary***")
                return self.cfg
            elif kwargs['type'] == object:
                self.logger.debug("***Returned as object***")
                return Struct(**self.cfg)
            else:
                self.logger.error("Illegal type specified for return")
                return None
        

        
    def format_args(self, args):
        """
        This is a helper function to parse for the correct files from a
        ConfigReader.read() function call

		Inputs:
                args = files to be parsed
		
		Outputs:
				files = set of files
        """
            
        files = []
        for arg in args:
            if type(arg) is list:
                for elem in arg:
                    if '*' in elem:
                        for file in glob.glob(elem):
                            print(file)
                            if "log" not in file:
                                files.append(file)
                    else:
                        files.append(elem)
            elif '*' in arg:
                for file in glob.glob(arg):
                    if "log" not in file:
                        files.append(file)
            else:
                files.append(arg)
        
        
        files = set(files)
        return files
        
        
        
    def format_values(self, value):
        """
        This function adds functionality for parsing special string cases.
        e.g. comma separated strings with newlines -> lists

		Inputs:
				value = variable to parse
		Outputs:
				fvalue = formatted variable
        """

        try:
            fvalue = ast.literal_eval(value)
        except:
            try:
                fvalue = SafeEval.safe_eval(value)
            except:
                fvalue = value
            
        if type(fvalue) is str:
            if ("'" in value or '"' in value) and ',' in value:
                pass
            else:
                if ',' in fvalue:
                    fvalue = self.pattern.split(fvalue)
            
        elif type(fvalue) is tuple or type(fvalue) is list:
            # Not yet implemented
            pass
        elif type(fvalue) is dict:
            # Not yet implemented
            pass
            
        return fvalue
        
        
    def parse_options(self, cfg, *args):
        """
        Parse options given from user
		
		Inputs:
				cfg = dictionary of config variables
				args = list of functions to apply to each value 

		Outputs:
				fcfg = formatted dictionary of config variables
        """
        fcfg = self.formatDict(cfg, args)
        
        return fcfg

       
    def formatDict(self, cfg, args, cfg_copy={}):
        """
        Format the output based on user options/functions
		Recursively goes through the list

		Inputs:
				cfg = dictionary of config variables
				args = options to parse by
				cfg_copy = recursively fill up dictionary
		Outputs:
				cfg_copy = returned formatted dictionary
		
        """
        
        for key, value in cfg.items():
            for option in args:
                try:
                    # Test line to see if it is a dict and/or ordereddict
                    value.keys
                    if 'key' in option:
                        cfg_copy[self.key_wrapper(key,option['key'])] = {}
                        self.formatDict(cfg[key],args,cfg_copy[self.key_wrapper(key,option['key'])])
                    else:
                        cfg_copy[key] = {}
                except:
                    if 'key' in option:
                        cfg_copy[self.key_wrapper(key,option['key'])] = value
                    else:
                        cfg_copy[key] = value
        

        return cfg_copy

        
    def key_wrapper(self, key, func):
        """
        This is a wrapper function used for various
        string function calls

		Inputs:
				key = key from config file to operate on
				func = function to use on the key
						Acceptable functions: lower, upper, ...
		
		Outputs:
				key = parsed key
		
        """
        if func:
            return getattr(key,func)()
        else:
            return key
                
        
