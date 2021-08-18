#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test_nenotr.py

Written by: Ashton Reimer

A script to unit test the CalcNeNoTr class.

CHANGELOG:

        19/07/2018 - Initial implementation

"""

import os
import sys
import tables
import unittest
import numpy as np
import shutil


# ensure we are running from the tests directory
path = os.getcwd()
if not os.path.basename(path) == 'tests':
    print('Need to run this in the tests directory.')
    sys.exit(1)

exclude_groups = ['/ProcessingParams/ComputerInfo','/ProcessingParams/ConfigFiles/File1']
exclude_nodes = ['/ProcessingParams/ProcessingTimeStamp','/ProcessingParams/SoftwareVersion',
                 '/Calibration/CalDate','/ProcessingParams/OutputPath','/ProcessingParams/RawFiles']

def validate(expected,newfile):
    # compare data in one file with data in another
    with tables.open_file(expected,'r') as f1:
        with tables.open_file(newfile,'r') as f2:

            for group in f1.walk_groups():
                group_name = group._v_pathname

                if group_name in exclude_groups:
                    continue                

                # First make sure the group in the input file is in the output file.
                f2_groups = [g._v_pathname for g in f2.walk_groups()]
                if group_name not in f2_groups:
                    print("Group: %s is missing from file: %s" % (group_name,self.file2))
                    print("Files are NOT the same...")
                    return False

                # Now list the nodes in the group. For any node that isn't a group, write it
                # to the output file.

                f2_nodes = [n for n in f2.list_nodes(group_name)]
                f2_node_names = np.array([x._v_pathname for x in f2_nodes])

                for f1_node in f1.list_nodes(group_name):
                    if not hasattr(f1_node,'_v_nchildren'):        #detect if it's a group or not

                        if f1_node._v_pathname in exclude_nodes:
                            continue

                        ind = np.where(f2_node_names == f1_node._v_pathname)[0]

                        if len(ind) != 1:
                            print("Node: %s NOT found in file: %s" % (f1_node._v_pathname,self.file2))
                            print("Files are NOT the same...")
                            return False
                        else:
                            ind = ind[0]

                        # Read the array from file1 and file2 and compare them!
                        data1 = f1_node.read()
                        data2 = f2_nodes[ind].read()

                        try:
                            np.testing.assert_array_equal(data1,data2)
                        except AssertionError:
                            print("Arrays: %s NOT equal!" % (f1_node._v_pathname))
                            return False

    print("Files are equivalent. Data in them is exactly the same.")
    return True

from resolvedalts import ResolveVectorsAlt
class SingleFileResolvedAlts(unittest.TestCase):
    def test_ResolvedAlts(self):

        # Run the nenotr program
        config = os.path.join(path,'config_test.ini')
        print(config)

        vvelsalt = ResolveVectorsAlt(config)
        vvelsalt.run()

        output_file = 'output/20150126.001_vvelsalt_5min.h5'
        expected_file = 'expected/20150126.001_vvelsalt_5min.h5'
        
        status = validate(expected_file,output_file)
        self.assertTrue(status)

        print("Cleaning up test files...")
        shutil.rmtree('output')


if __name__ == '__main__':
    unittest.main()

