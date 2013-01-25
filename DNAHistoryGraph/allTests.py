#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt
import unittest

from pyAVG.DNAHistoryGraph.graphTest import DNAHistoryGraphTest
from pyAVG.DNAHistoryGraph.segmentTest import SegmentTest

from cactus.shared.test import parseCactusSuiteTestOptions

def allSuites(): 
    allTests = unittest.TestSuite((unittest.makeSuite(DNAHistoryGraphTest, 'test'),
                                   unittest.makeSuite(SegmentTest, 'test')))
    return allTests
        
def main():
    suite = allSuites()
    runner = unittest.TextTestRunner()
    i = runner.run(suite)
    return len(i.failures) + len(i.errors)
        
if __name__ == '__main__':
    import sys
    sys.exit(main())
                
