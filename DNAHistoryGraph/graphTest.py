import unittest
from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph

class DNAHistoryGraphTest(unittest.TestCase):
    """Tests the DNA history graph functions, particularly the acyclicity functions
    """
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.g = DNAHistoryGraph()
        self.s1 = self.g.newSegment("A")
        self.s2 = self.g.newSegment("T")
        self.s3 = self.g.newSegment("T")
        self.s4 = self.g.newSegment("A")
        self.s5 = self.g.newSegment("A")
        self.g.createBranch(self.s1, self.s2)
        self.g.createBranch(self.s1, self.s3)
        self.g.createBranch(self.s3, self.s4)
        self.g.createBranch(self.s3, self.s5)
        
        self.s1b = self.g.newSegment("A")
        self.s2b = self.g.newSegment("A")
        self.s3b = self.g.newSegment()
        self.s4b = self.g.newSegment("T")
        self.s5b = self.g.newSegment("T")
        self.g.createBranch(self.s1b, self.s2b)
        self.g.createBranch(self.s1b, self.s3b)
        self.g.createBranch(self.s3b, self.s4b)
        self.g.createBranch(self.s3b, self.s5b)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testInterpolateSegment(self):
        #Test case where segment is root
        self.assertEqual(self.s1.parent, None)
        s6 = self.g.interpolateSegment(self.s1)
        self.assertEqual(self.s1.parent, s6)
        self.assertEqual(s6.children, set([ self.s1 ]))
        #Test case where segment is not root
        self.assertEqual(self.s2.parent, self.s1)
        self.assertEqual(self.s1.children, set([ self.s2, self.s3 ]))
        s6 = self.g.interpolateSegment(self.s2)
        self.assertEqual(self.s2.parent, s6)
        self.assertEqual(s6.children, set([ self.s2 ]))
        self.assertEqual(self.s1.children, set([ s6, self.s3 ]))
    
    def testPullDown(self):
        s6 = self.g.pullDown(self.s1, [ self.s2, self.s3 ])
        self.assertEqual(s6.parent, self.s1)
        self.assertEqual(self.s1.children, set([ s6 ]))
        self.assertEqual(s6.children, set([self.s2, self.s3 ]))
        self.assertEqual(s6, self.s2.parent)
        self.assertEqual(s6, self.s3.parent)
        
    def testCreateBond(self):
        self.g.createBond(self.s1.left, self.s1b.left)
        self.assertEquals(self.s1.left.bond, self.s1b.left)
        self.assertEquals(self.s1b.left.bond, self.s1.left)
        
        self.g.createBond(self.s3.left, self.s4b.left) 
        self.assertEquals(self.s3.left.bond, self.s4b.left)
        self.assertEquals(self.s4b.left.bond, self.s3.left)
        
        try: #Try to make bond that violates acyclicity 
            self.g.createBond(self.s4.left, self.s3b.left) 
            self.assertTrue(False)
        except RuntimeError:
            pass
    
    def testDeleteBond(self):
        self.g.createBond(self.s1.left, self.s1b.left)
        self.g.createBond(self.s2.left, self.s2b.left)
        self.g.createBond(self.s3.left, self.s3b.left)
        self.g.createBond(self.s4.left, self.s4b.left)
        self.g.createBond(self.s5.left, self.s5b.left)
        
        self.assertEquals(self.s3.left.bond, self.s3b.left)
        self.assertEquals(self.s3b.left.bond, self.s3.left)
        self.g.deleteBond(self.s3.left)
        self.assertEquals(self.s3.left.bond, None)
        self.assertEquals(self.s3b.left.bond, None)
        
        self.g.deleteBond(self.s4.left)
        self.g.deleteBond(self.s5.left)
        
        self.g.createBond(self.s4.left, self.s3b.left)
        self.assertEquals(self.s4.left.bond, self.s3b.left)
        self.assertEquals(self.s3b.left.bond, self.s4.left)
        
        self.g.deleteBond(self.s4.left)
        self.g.createBond(self.s3.left, self.s4b.left) 
        self.assertEquals(self.s3.left.bond, self.s4b.left)
        self.assertEquals(self.s4b.left.bond, self.s3.left)
    
    def testAreSiblings(self):
        pass
    
    def testSubstitutionAmbiguity(self):
        self.assertEquals(self.g.substitutionAmbiguity(), 3)
    
    def testRearrangementAmbiguity(self):
        pass
    
    def testAmbiguity(self):
        pass
    
    def testLowerBoundSubstitutionCost(self):
        self.assertEquals(self.g.lowerBoundSubstitutionCost(), 3)
    
    def testUpperBoundSubstitutionCost(self):
        self.assertEquals(self.g.upperBoundSubstitutionCost(), 6)
        
if __name__ == '__main__':
    unittest.main()