
"""
Multialign module

Perform multiple sequence alignment using tree as guide

Written by Marshall Beddoe <mbeddoe@baselineresearch.net>
Copyright (c) 2004 Baseline Research

Licensed under the LGPL
"""

import align
from Numeric import *

class Multialign:

    """Implementation of multialign base class"""

    def __init__(self, tree):
        self.tree = tree
        self.aligned = []
        self.index = 0

        self._go()

    def _go(self):
        pass

    def __len__(self):
        return len(self.aligned)

    def __getitem__(self, index):
        return self.aligned[index]

    def __iter__(self):
        self.index = 0
        return self

    def next(self):
        if self.index == len(self.aligned):
            raise StopIteration

        self.index += 1

        return self.aligned[self.index - 1]

class NeedlemanWunsch(Multialign):

    """Perform global multiple sequence alignment"""
    def __init__(self, tree, smx=None):
        self.smx = smx

        # If similarity matrix is None, make a quick identity matrix
        if self.smx == None:
            self.smx = zeros((257, 257), Float)

            for i in range(257):
                for j in range(257):
                    if i == j:
                        self.smx[i][j] = 1.0
                    else:
                        self.smx[i][j] = 0.0

        Multialign.__init__(self, tree)

    def _go(self):

        self._assign(self.tree)
        self._alignSum(self.tree, [])

    def _assign(self, root):
        """Traverse tree and align sequences"""

        if root.getValue()[2] == None:
            if root.getLeft().getValue()[2] == None:
                self._assign(root.getLeft())

            if root.getRight().getValue()[2] == None:
                self._assign(root.getRight())

            # Get sequences
            seq1 = root.getLeft().getValue()[2][1]
            seq2 = root.getRight().getValue()[2][1]

            (a1, a2, e1, e2, score, gaps) = \
                align.NeedlemanWunsch(seq1, seq2, self.smx, 0, 0)

            v1 = root.getLeft().getValue()
            v2 = root.getRight().getValue()

            nv1 = (v1[0], v1[1], v1[2], e1)
            nv2 = (v2[0], v2[1], v2[2], e2)

            root.getLeft().setValue(nv1)
            root.getRight().setValue(nv2)

            # Choose the sequence with least gaps
            if e1 < e2:
                nseq = a1
            else:
                nseq = a2

            v1 = root.getValue()
            nseq = (v1[0], nseq)
            nv1 = (v1[0], v1[1], nseq, v1[3])
            root.setValue(nv1)

    def _alignSum(self, root, edits):

        if root.getLeft() == None and root.getRight() == None:
            seq1 = root.getValue()[2][1]
            id = root.getValue()[2][0]
            new = seq1

            for i in range(len(edits)):
                e = edits[i]
                self._applyEdits(new, e)

            self.aligned.append((id, new))
        else:
            e = root.getLeft().getValue()[3]
            edits.insert(0, e)
            self._alignSum(root.getLeft(), edits)
            k = edits.pop(0)

            e = root.getRight().getValue()[3]
            edits.insert(0, e)
            self._alignSum(root.getRight(), edits)
            k = edits.pop(0)

    def _applyEdits(self, seq, edits):
        i = 0
        gap = 256

        edits.sort()

        for e in edits:
            seq.insert(e, gap)

        return seq
