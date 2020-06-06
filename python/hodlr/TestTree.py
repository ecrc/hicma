import unittest
from Tree import Tree

class TestTree(unittest.TestCase):
    superms=[0,1,2,3,4,5,6]
    lowersubms=[1,3,5,7,9,11,13]
    uppersubms=[2,4,6,8,10,12,14]
    def super_sub(self):
        for i, s in enumerate(self.superms):
            ls = Tree.lower_subm(s)
            self.assertEqual(ls, self.lowersubms[i])
            us = Tree.upper_subm(s)
            self.assertEqual(us, self.uppersubms[i])
            superml = Tree.super_matrix(ls)
            self.assertEqual(superml, s)
            supermu = Tree.super_matrix(us)
            self.assertEqual(supermu, s)
    numtotalrow = 16 # number of rows of matrix
    tss = [ 16, 8, 4, 2, 1]
    def tile_size(self):
        for l in range(0,5):
            self.assertEqual(Tree.max_tile_size(l, self.numtotalrow), self.tss[l])
    depths = [0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]
    def depth(self):
        for i, n in enumerate(range(len(self.depths))):
            self.assertEqual(Tree.depth(i), self.depths[i])
    numrows = [16, 8, 8, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2]
    def number_of_rows(self):
        for i, n in enumerate(range(len(self.numrows))):
            self.assertEqual(Tree.numrow(i, self.numtotalrow), self.numrows[i])
    tiles_at_levels = [[0],[1,2],[3,4,5,6],[7,8,9,10,11,12,13,14]]
    def tiles_at_depth(self):
        for level in range(0,4):
            #print(level, Tree.tiles_at_depth(level))
            tiles = Tree.tiles_at_depth(level)
            self.assertEqual(tiles, self.tiles_at_levels[level])
        tiles = Tree.tiles_at_depth_before(3, 11)
        self.assertEqual(tiles, [7,8,9,10])
        tiles = Tree.tiles_at_depth_before(3, 14)
        self.assertEqual(tiles, [7,8,9,10,11,12,13])
        tiles = Tree.tiles_at_depth_before(3, 15)
        self.assertEqual(tiles, [7,8,9,10,11,12,13,14])
        tiles = Tree.tiles_at_depth_before(3, 16)
        self.assertEqual(tiles, [7,8,9,10,11,12,13,14])
    def tile_storage_size(self):
        tree = Tree(16, 3)
        self.assertEqual(tree.stsize(10), 6)
    ststarts = [None, 0, 8, 16, 20, 24, 28, 32, 34, 36, 38, 40, 42, 44, 46]
    def tile_storage_start_end(self):
        tree = Tree(16, 1)
        for i in range(1, 15):
            self.assertEqual(tree.ststart(i), self.ststarts[i])
            self.assertEqual(tree.ststart(i)+tree.stsize(i), tree.stend(i))
    rowstarts = [None, 8, 0, 4, 0, 12, 8, 2, 0, 6, 4, 10, 8, 14, 12]
    def tile_row_start_end(self):
        tree = Tree(16, 1)
        for i in range(1, 15):
            self.assertEqual(tree.rowstart(i), self.rowstarts[i])
            self.assertEqual(tree.rowstart(i)+Tree.numrow(i, 16), tree.rowend(i))
    colstarts = [None, 0, 8, 0, 4, 8, 12, 0, 2, 4, 6, 8, 10, 12, 14]
    def tile_col_start_end(self):
        tree = Tree(16, 1)
        for i in range(1, 15):
            self.assertEqual(tree.colstart(i), self.colstarts[i])
            self.assertEqual(tree.colstart(i)+Tree.numcol(i, 16), tree.colend(i))
if __name__ == '__main__':
    test = TestTree()
    test.super_sub()
    test.tile_size()
    test.depth()
    test.number_of_rows()
    test.tiles_at_depth()
    test.tile_storage_size()
    test.tile_storage_start_end()
    test.tile_row_start_end()
    test.tile_col_start_end()
