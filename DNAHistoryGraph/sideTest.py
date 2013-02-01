import unittest
from pyAVG.DNAHistoryGraph.segment import Segment
from pyAVG.DNAHistoryGraph.label import Label

class SideTest(unittest.TestCase):
    """Tests the segment and associated label classes, including testing label functionality, but ignoring sides, which are 
    tested in sideTest.py"""
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        #Creates a very basic branch-tree
        self.s1 = Segment().left
        self.s2 = Segment().right
        self.s3 = Segment().left
        self.s4 = Segment().right
        self.s1.segment.createBranch(self.s2.segment)
        self.s1.segment.createBranch(self.s3.segment)
        self.s3.segment.createBranch(self.s4.segment)
        
        self.s1B = Segment().left
        self.s2B = Segment().right
        self.s4B = Segment().right
        self.s1B.segment.createBranch(self.s2.segment)
        
        self.s1.createBond(self.s1B)
        self.s2.createBond(self.s2B)
        self.s4.createBond(self.s4B)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testChildren(self):
        self.assertEqual(self.s1.children(), set([ self.s2, self.s3 ]))
        self.assertEqual(self.s2.children(), set())
        self.assertEqual(self.s3.children(), set([ self.s4 ]))
        self.assertEqual(self.s4.children(), set())
        
    def testParent(self):
        self.assertEqual(self.s1.parent(), None)
        self.assertEqual(self.s2.parent(), self.s1)
        self.assertEqual(self.s3.parent(), self.s1)
        self.assertEqual(self.s4.parent(), self.s3)
    
    def testCreateAndDeleteBond(self):
        self.assertEqual(self.s1.bond, self.s1B)
        self.assertEqual(self.s1B.bond, self.s1)
        self.s1.deleteBond()
        self.assertEqual(self.s1.bond, None)
        self.assertEqual(self.s1B.bond, None)
        
    def testAncestor(self):
        self.assertEqual(self.s1.ancestor(), self.s1)
        self.assertEqual(self.s2.ancestor(), self.s1)
        self.assertEqual(self.s3.ancestor(), self.s1)
        self.assertEqual(self.s4.ancestor(), self.s1)
        self.s3.createBond(Segment().left)
        self.assertEqual(self.s4.ancestor(), self.s3)
        
    def testLiftedBonds(self):
        """The lifted bond function returns the segments that contain the lifting Bonds
        """
        self.assertEqual(self.s1.liftedBonds(), set([ self.s2, self.s4 ]))
        self.assertEqual(self.s2.liftedBonds(), set())
        self.assertEqual(self.s3.liftedBonds(), set([ self.s4 ]))
        self.assertEqual(self.s4.liftedBonds(), set())
    
    def testNonTrivialLiftedBonds(self):
        self.assertEqual(self.s1.nonTrivialLiftedBonds(), set([ self.s4 ]))
        s5 = Segment(parent=self.s1.segment).left
        s5.createBond(Segment())
        self.assertEqual(self.s1.nonTrivialLiftedBonds(), set([ self.s4, s5 ]))
    
    def testRearrangementAmbiguity(self):
        self.assertEqual(self.s1.rearrangementAmbiguity(), 0)
        self.assertEqual(self.s2.rearrangementAmbiguity(), 0)
        self.assertEqual(self.s3.rearrangementAmbiguity(), 0)
        self.assertEqual(self.s4.rearrangementAmbiguity(), 0)
        s5 = Segment(parent=self.s1.segment).left
        s5.createBond(Segment())
        self.assertEqual(self.s1.rearrangementAmbiguity(), 1)

if __name__ == '__main__':
    unittest.main()