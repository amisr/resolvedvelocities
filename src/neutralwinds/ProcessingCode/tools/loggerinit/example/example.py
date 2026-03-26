"""
Example.py

Written by: Steven Chen

This class tests the functionality of the LoggerInit package

CHANGELOG:

    5/6/14 - First implementation
    
"""
# Add to PYTHONPATH
import sys
sys.path.append("../")

# Import user module
import LoggerInit as LoggerInit

if __name__ == "__main__":
    
    # Initializes LoggerInit Class
    LoggerInit.LoggerInit(config="testconfig/log.ini")
    