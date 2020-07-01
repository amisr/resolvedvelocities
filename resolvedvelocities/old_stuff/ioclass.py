#!/usr/bin/env python

"""

"""

import os
import tables
import numpy as np
import distutils.dir_util

class OutputFileClass():
    
    def __init__(self):
        """ initialization function """
        self.h5Paths={}
        self.h5Attribs = {}
        self.fname = ''
        self.title = ''
        self.fhandle = None
        return
    
    def createFile(self, fname):
        """ create file fname, set self.fname and self.fhandle """
        self.fname = fname
        try:
            distutils.dir_util.mkpath(os.path.dirname(self.fname))
        except:
            raise IOError('Unable to create output path %s' % os.path.dirname(self.fname))
        self.fhandle=tables.open_file(self.fname, mode = "w", title = self.title)
        return
    
    def openFile(self):
        """ open file self.fname """
        self.fhandle = tables.open_file(self.fname, mode = "a")
        return
        
    def closeFile(self):
        """ close self.fhandle """
        self.fhandle.close()
        
    def createh5groups(self):
        """ creates groups """
        tvals = sorted(self.h5Paths.values())
        for v0,v1 in tvals:
            gp,gn = os.path.split(v0)      
            self.fhandle.create_group(gp,gn,v1)
        return
    
    def createStaticArray(self,path,data,keys2do=[]):  
        """ creates a static array """
        if len(keys2do)==0:
            dp,dn = os.path.split(path)
            self.fhandle.create_array(dp,dn,data,'Static array')
        else:
            for key in keys2do:
                self.fhandle.create_array(path,key,np.array(data[key]),'Static array')
        return
    
    def createDynamicArray(self,fhandle,path,rec,keys2do=[]):  
        """ creates a dynamic array """
        if len(keys2do)==0:
            dp,dn = os.path.split(path)
            data = rec.copy()
            data.shape = (1,)+data.shape  ## add integration dimension to data array
            if not fhandle.__contains__(path):
                shape = list(data.shape)
                shape[0] = 0
                atom = tables.Atom.from_dtype(data.dtype)
                arr = fhandle.create_earray(dp,dn,atom,shape)
                arr.flavor='numpy'
            arr = fhandle.getNode(path)
            if (len(arr.shape)>2) and (data.shape[2] != arr.shape[2]):
                if data.shape[2] > arr.shape[2]:
                    # read array 
                    tarr = arr.read() 
                    # remove old node                 
                    arr.remove() 
                    tshape=list(tarr.shape); tshape[2]=data.shape[2]-tarr.shape[2]
                    tarr=np.append(tarr,np.zeros(tshape)*np.nan,axis=2)   
                    # create new node                 
                    shape = list(tarr.shape)
                    shape[0] = 0
                    atom = tables.Atom.from_dtype(tarr.dtype)
                    arr = fhandle.create_earray(dp,dn,atom,shape)
                    arr.flavor='numpy'
                    arr = fhandle.getNode(path)
                    # dump data
                    arr.append(tarr)
                else:
                    tshape=list(data.shape); tshape[2]=arr.shape[2]-data.shape[2]
                    data=np.append(data,np.zeros(tshape)*np.nan,axis=2)
            arr.append(data)
        else:
            for key in keys2do:
                data = np.array(rec[key])
                data.shape = (1,)+data.shape  ## add integration dimension to data array
                if not fhandle.__contains__(path+'/'+key):
                    shape = list(data.shape)
                    shape[0] = 0
                    atom = tables.Atom.from_dtype(data.dtype)
                    arr = fhandle.create_earray(path,key,atom,shape)
                    arr.flavor='numpy'
                arr = fhandle.getNode(path+'/'+key)
                if (len(arr.shape)>2) and (data.shape[2] != arr.shape[2]):
                    if data.shape[2] > arr.shape[2]:
                        # read array  
                        tarr = arr.read() 
                        # remove old node                 
                        arr.remove() 
                        tshape=list(tarr.shape); tshape[2]=data.shape[2]-tarr.shape[2]
                        tarr=np.append(tarr,np.zeros(tshape)*np.nan,axis=2)   
                        # create new node                 
                        shape = list(tarr.shape)
                        shape[0] = 0
                        atom = tables.Atom.from_dtype(tarr.dtype)
                        arr = fhandle.create_earray(path,key,atom,shape)
                        arr.flavor='numpy'
                        arr = fhandle.getNode(path+'/'+key)
                        # dump data
                        arr.append(tarr)
                    else:
                        tshape=list(data.shape); tshape[2]=arr.shape[2]-data.shape[2]
                        data=np.append(data,np.zeros(tshape)*np.nan,axis=2)
                arr.append(data)        
        return
    
    def setAtrributes(self):
        for key in self.h5Attribs.keys():
            for attr in self.h5Attribs[key]:
                try:  self.fhandle.set_node_attr(key,attr[0],attr[1])
                except: ''
        return
