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
        self.s2 = Segment().left
        self.s3 = Segment().left
        self.s4 = Segment().left
        self.s5 = Segment().left
        self.s1.segment.createBranch(self.s2.segment)
        self.s1.segment.createBranch(self.s3.segment)
        self.s3.segment.createBranch(self.s4.segment)
        self.s3.segment.createBranch(self.s5.segment)
        
        self.s1B = Segment().left
        self.s2B = Segment().left
        self.s3B = Segment().left
        self.s4B = Segment().left
        self.s5B = Segment().left
        self.s1B.segment.createBranch(self.s2B.segment)
        self.s1B.segment.createBranch(self.s3B.segment)
        self.s3B.segment.createBranch(self.s4B.segment)
        self.s5B.segment.createBranch(self.s5B.segment)
        
        self.s1.createBond(self.s1B)
        self.s2.createBond(self.s2B)
        self.s4.createBond(self.s4B)
        self.s5.createBond(self.s5B.opposite)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testOpposite(self):
        self.assertEquals(self.s1.opposite, self.s1.segment.right)
        self.assertEquals(self.s2.opposite, self.s2.segment.right)
        
    def testChildren(self):
        self.assertEqual(set(self.s1.children()), set([ self.s2, self.s3 ]))
        self.assertEqual(set(self.s2.children()), set())
        self.assertEqual(set(self.s3.children()), set([ self.s4, self.s5 ]))
        self.assertEqual(set(self.s4.children()), set())
        
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
        self.assertEqual(self.s5.ancestor(), self.s1)
        self.s3.createBond(Segment().left)
        self.assertEqual(self.s4.ancestor(), self.s3)
        self.assertEqual(self.s5.ancestor(), self.s3)
        
    def testIsJunction(self):
        self.assertTrue(self.s1.isJunction())
        self.assertTrue(not self.s2.isJunction())
        self.assertTrue(self.s3.isJunction())
        self.assertTrue(not self.s4.isJunction())
        self.assertTrue(not self.s5.isJunction())
        
    def testLiftedBonds(self):
        """The lifted bond function returns the segments that contain the lifting Bonds
        """
        self.assertEqual(self.s1.liftedBonds(), set([ self.s2, self.s4, self.s5 ]))
        self.assertEqual(self.s2.liftedBonds(), set())
        self.assertEqual(self.s3.liftedBonds(), set([ self.s4, self.s5 ]))
        self.assertEqual(self.s4.liftedBonds(), set())
        self.assertEqual(self.s5.liftedBonds(), set())
    
    def testNonTrivialLiftedBonds(self):
        self.assertEqual(self.s1.nonTrivialLiftedBonds(), set([ self.s4, self.s5 ]))
        self.assertEqual(self.s3.nonTrivialLiftedBonds(), set([ self.s4, self.s5 ]))
        self.assertEqual(self.s4.nonTrivialLiftedBonds(), set())
        
    def testRearrangementAmbiguity(self):
        self.assertEqual(self.s1.rearrangementAmbiguity(), 1)
        self.assertEqual(self.s2.rearrangementAmbiguity(), 0)
        self.assertEqual(self.s3.rearrangementAmbiguity(), 0)
        self.s3.createBond(Segment().left)
        self.assertEqual(self.s3.rearrangementAmbiguity(), 1)
        self.assertEqual(self.s4.rearrangementAmbiguity(), 0)

if __name__ == '__main__':
    unittest.main()