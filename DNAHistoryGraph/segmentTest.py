import unittest
from pyAVG.DNAHistoryGraph.segment import Segment
from pyAVG.DNAHistoryGraph.label import Label

class SegmentTest(unittest.TestCase):
    """Tests the segment and associated label classes, including testing label functionality, but ignoring sides, which are 
    tested in sideTest.py"""
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        #Creates a very basic branch-tree
        self.s1 = Segment("A")
        self.s2 = Segment("A")
        self.s3 = Segment()
        self.s4 = Segment("T")
        self.s1.createBranch(self.s2)
        self.s1.createBranch(self.s3)
        self.s3.createBranch(self.s4)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testChildren(self):
        self.assertEqual(self.s1.children, set([ self.s2, self.s3 ]))
        self.assertEqual(self.s2.children, set())
        self.assertEqual(self.s3.children, set([ self.s4 ]))
        self.assertEqual(self.s4.children, set())
        
    def testParent(self):
        self.assertEqual(self.s1.parent, None)
        self.assertEqual(self.s2.parent, self.s1)
        self.assertEqual(self.s3.parent, self.s1)
        self.assertEqual(self.s4.parent, self.s3)
    
    def testDeleteBranch(self):
        self.s1.deleteBranch(self.s2)
        self.assertEqual(self.s2.parent, None)
        self.assertEqual(self.s1.children, set([ self.s3 ]))
        self.s1.deleteBranch(self.s3)
        self.assertEqual(self.s3.parent, None)
        self.assertEqual(self.s1.children, set())
        
    def testDisconnect(self):
        self.s3.disconnect()
        self.assertEqual(self.s3.parent, None)
        self.assertEqual(self.s2.parent, self.s1)
        self.assertEqual(self.s4.parent, self.s1)
        self.assertEqual(self.s1.children, set([ self.s2, self.s4 ]))
        self.assertEqual(self.s3.children, set())
    
    def testStr(self): 
        self.assertEqual(str(self.s1), "A")
        self.assertEqual(str(self.s2), "A")
        self.assertEqual(str(self.s3), str(None))
        self.assertEqual(str(self.s4), "T")
    
    def testLabel(self): 
        self.assertEqual(self.s1.label, Label("A"))
        self.assertEqual(self.s2.label, Label("A"))
        self.assertEqual(self.s3.label, None)
        self.assertEqual(self.s4.label, Label("T"))
    
    def testDeleteLabel(self): 
        self.s1.deleteLabel()
        self.assertEqual(self.s1.label, None)
        
    def testAncestor(self):
        self.assertEqual(self.s1.ancestor(), self.s1)
        self.assertEqual(self.s2.ancestor(), self.s1)
        self.assertEqual(self.s3.ancestor(), self.s1)
        self.assertEqual(self.s4.ancestor(), self.s1)
        self.s3.label = Label("A")
        self.assertEqual(self.s4.ancestor(), self.s3)
        
    def testLiftedLabels(self):
        """The lifted label function returns the segments that contain the lifting labels
        """
        self.assertEqual(self.s1.liftedLabels(), set([ self.s2, self.s4 ]))
        self.assertEqual(self.s2.liftedLabels(), set())
        self.assertEqual(self.s3.liftedLabels(), set([ self.s4 ]))
        self.assertEqual(self.s4.liftedLabels(), set())
    
    def testNonTrivialLiftedLabels(self):
        self.assertEqual(self.s1.nonTrivialLiftedLabels(), set([ self.s4 ]))
        s5 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.nonTrivialLiftedLabels(), set([ self.s4, s5 ]))
    
    def testSubstitutionAmbiguity(self):
        self.assertEqual(self.s1.substitutionAmbiguity(), 0)
        self.assertEqual(self.s2.substitutionAmbiguity(), 0)
        self.assertEqual(self.s3.substitutionAmbiguity(), 0)
        self.assertEqual(self.s4.substitutionAmbiguity(), 0)
        s5 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.substitutionAmbiguity(), 1)
        
    def testLowerBoundSubstitutionCost(self):
        s5 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.lowerBoundSubstitutionCost(), 1)
        self.assertEqual(self.s2.lowerBoundSubstitutionCost(), 0)
        self.assertEqual(self.s3.lowerBoundSubstitutionCost(), 0)
        self.assertEqual(self.s4.lowerBoundSubstitutionCost(), 0)
        s6 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.lowerBoundSubstitutionCost(), 1)
        s7 = Segment("G", parent=self.s1)
        self.assertEqual(self.s1.lowerBoundSubstitutionCost(), 2)
    
    def testUpperBoundSubstitutionCost(self):
        s5 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.upperBoundSubstitutionCost(), 2)
        self.assertEqual(self.s2.upperBoundSubstitutionCost(), 0)
        self.assertEqual(self.s3.upperBoundSubstitutionCost(), 0)
        self.assertEqual(self.s4.upperBoundSubstitutionCost(), 0)
        s6 = Segment("T", parent=self.s1)
        self.assertEqual(self.s1.upperBoundSubstitutionCost(), 3)
        s7 = Segment("G", parent=self.s1)
        self.assertEqual(self.s1.upperBoundSubstitutionCost(), 4)
    
if __name__ == '__main__':
    unittest.main()