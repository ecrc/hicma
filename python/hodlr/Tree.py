import math
class Tree(object):
    # Tree class provides utilities for a full binary tree.
    # Warnings:
    # 1. Be aware that divisions like /2 should result in natural numbers
    #    if tree is full, o.w., wrong result is produced
    #
    elmsize = 1 # size of double float etc
    def __init__(self, numtotalrows, stmax, elmsize = 1):
        self.numrows = numtotalrows
        self.numcols = self.numrows #ASSUMPTION: square tiles
        self.stmax = stmax
        self.elmsize = elmsize
    # static methods
    def lower_subm(i):
        return 2*i+1
    def upper_subm(i):
        return 2*i+2
    def super_matrix(i):
        fr = math.ceil(float(i)/2.0)-1
        ir = int(fr)
        assert(ir == fr)
        return ir
    def depth(i):
        return math.floor(math.log(i+1, 2))
    def max_tile_size(level, numrows): #ts
        fr = float(numrows)/math.pow(2, level)
        ir = int(fr)
        assert(ir == fr)
        return ir
    def numrow(i, numrows): #ASSUMPTION: uniform tile size
        return Tree.max_tile_size(Tree.depth(i), numrows)
    def numcol(i, numcols):
        return Tree.numrow(i, numcols) #ASSUMPTION
    def tiles_at_depth(d):
        start = int(math.pow(2, d)-1)
        end   = int(math.pow(2, d+1)-2)
        return list(range(start, end+1))
    def tiles_at_depth_before(d, i):
        start = int(math.pow(2, d)-1)
        end   = int(math.pow(2, d+1)-2)
        if i < end+1:
            end = i-1
        return list(range(start, end+1))
    # class methods
    def stsize(self, i):
        return Tree.numrow(i, self.numrows) * self.stmax * self.elmsize
    def ststart(self, i):
        if i == 0:  # tile 0
            raise ValueError("Tile 0 is not valid for ststart()")
        depth = Tree.depth(i)
        # skip previous levels
        start = (depth-1) * self.numrows * self.stmax * self.elmsize
        tiles = Tree.tiles_at_depth_before(depth, i)
        for t in tiles:
            start += self.stsize(t)
        return start
    def stend(self, i):
        if i == 0:  # tile 0
            raise ValueError("Tile 0 is not valid for ststart()")
        return self.ststart(i) + self.stsize(i)
    def rowstart(self, i):
        if i%2 == 0:
            before = i - 1
            start = 0
        else:
            before = i
            start = Tree.numrow(i+1, self.numrows)
        depth = Tree.depth(i)
        tiles = Tree.tiles_at_depth_before(depth, before)
        for t in tiles:
            start += Tree.numrow(t, self.numrows)
        return start
    def rowend(self, i):
        return self.rowstart(i) + Tree.numrow(i, self.numrows)
    def colstart(self, i):
        depth = Tree.depth(i)
        tiles = Tree.tiles_at_depth_before(depth, i)
        start = 0
        for t in tiles:
            start += Tree.numcol(t, self.numcols)
        return start
    def colend(self, i):
        return self.colstart(i) + Tree.numcol(i, self.numcols)


