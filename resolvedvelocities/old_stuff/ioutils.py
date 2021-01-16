#!/usr/bin/env python

"""

"""

import tables

class h5file():

    def __init__(self,fname):
        """ initialization function """
        self.fname = fname
        self.fhandle = None        
        return
    
    def openFile(self):
        """ open file self.fname """
        self.fhandle = tables.open_file(self.fname, mode = "a")
        return
    
    def readWholeh5file(self):
        
        h5file=tables.open_file(self.fname)
        output={}
        for group in h5file.walk_groups("/"):
            output[group._v_pathname]={}
            for array in h5file.list_nodes(group, classname = 'Array'):
                output[group._v_pathname][array.name]=array.read()
        h5file.close()
         
        return output
        
def ini_tool(config,secName,parmName,required=0,defaultParm=''):

    try: 
        if config.has_option(secName,parmName):
            parm=config.get(secName,parmName)
        elif required==1:
            raise IOError('%s must have parameter %s' % (secName,parmName))
        else:
            parm=defaultParm
    except:
        raise IOError('Error reading %s from %s' % (parmName,secName))
        
    return parm
