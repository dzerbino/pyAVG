import unittest
from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph

class DNAHistoryGraphTest(unittest.TestCase):
    """Tests the DNA history graph functions, particularly the acyclicity functions
    """
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.g = DNAHistoryGraph()
        self.s1 = self.g.newSegment()
        self.s2 = self.g.newSegment()
        self.s3 = self.g.newSegment()
        self.s4 = self.g.newSegment()
        self.g.createBranch(self.s1, self.s2)
        self.g.createBranch(self.s1, self.s3)
        self.g.createBranch(self.s1, self.s4)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testInterpolateSegment(self):
        #Test case where segment is root
        self.assertEqual(self.s1.parent, None)
        s5 = self.g.interpolateSegment(self.s1)
        self.assertEqual(self.s1.parent, s5)
        self.assertEqual(s5.children, set([ self.s1 ]))
        #Test case where segment is not root
        self.assertEqual(self.s2.parent, self.s1)
        self.assertEqual(self.s1.children, set([ self.s2, self.s3, self.s4 ]))
        s6 = self.g.interpolateSegment(self.s2)
        self.assertEqual(self.s2.parent, s6)
        self.assertEqual(s6.children, set([ self.s2 ]))
        self.assertEqual(self.s1.children, set([ s6, self.s3, self.s4 ]))
    
    def testPullDown(self):
        s5 = self.g.pullDown(self.s1, [ self.s2, self.s3 ])
        self.assertEqual(s5.parent, self.s1)
        self.assertEqual(self.s1.children, set([ s5, self.s4 ]))
        self.assertEqual(s5.children, set([self.s2, self.s3 ]))
        self.assertEqual(s5, self.s2.parent)
        self.assertEqual(s5, self.s3.parent)
        
if __name__ == '__main__':
    unittest.main()