#!/usr/bin/env python
# Copyright (C) 2012 Ion Torrent Systems, Inc. All Rights Reserved

import unittest
import sys
import os
import subprocess

class Analysis_BasicTests(unittest.TestCase):


    def test_SanityCheck(self):
	# a very simple test to see if Analysis responds when called at the command line
        self.assertTrue(subprocess.call("Analysis | grep 'Command line = Analysis'", shell=True)==0,"Analysis failed to return info to stdout... compilation error?")


if __name__ == "__main__":
    try:
        import xmlrunner
        test_runner = xmlrunner.XMLTestRunner(stream=sys.stdout,output='test-reports')

    except ImportError:
        test_runner = unittest.TextTestRunner(stream=sys.stdout)

    unittest.main(testRunner=test_runner)

