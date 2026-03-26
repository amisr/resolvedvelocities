"""
Example.py

Written by: Steven Chen

This class tests the functionality of the ConfigReader package

CHANGELOG:

    5/6/14 - First implementation
    
"""
# Add to PYTHONPATH

# Import user module
from lib.configreader.ConfigReader import *

if __name__ == "__main__":
    
    # Initializes ConfigReader Class
    config = ConfigReader(logger_level='DEBUG')
    
    cfg = config.read("testconfig/*.ini", key='lower')
    
    cfg_obj = config.read(["testconfig/*.ini","testconfig/smtp.ini"],type=object, key='lower')
    