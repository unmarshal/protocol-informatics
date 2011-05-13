
"""
Phylogenetic tree module

Implementation of multiple tree building algorithms:
    - UPGMA
    - Maximum Parsimony
    - Maximum Likelihood

Also, contains class to perform clustering based on results of tree generation

Written by Marshall Beddoe <mbeddoe@baselineresearch.net>
Copyright (c) 2004 Baseline Research

Licensed under the LGPL
"""

from sets import *
from tree import *
from pydot import *

class Phylogeny:

    """Implementation of base phylogenetic class"""

    def __init__(self, sequences, dmx, minval=1.0):
        self.dmx = dmx
        self.index = 0
        self.tree = None
        self.clusters = []
        self.minval = minval
        self.sequences = sequences

        self._go()

    def __len__(self):
        return len(self.clusters)

    def __iter__(self):
        self.index = 0
        return self

    def next(self):
        if self.index == len(self.clusters):
            raise StopIteration

        self.index += 1

        return self.clusters[self.index - 1]

    def __getitem__(self, index):
        return self.clusters[index]

    def _go(self):
        """Perform tree construction"""
        pass

class UPGMA(Phylogeny):

    """UPGMA tree construction method"""

    def _go(self):

        # Universal set
        Cu = Set()

        # Place each sequence into individual tree node
        for i in range(len(self.sequences)):
            ntree = Tree()
            ntree.setValue((i, 0, self.sequences[i], None))
            Cu.add(ntree)

        n = len(Cu) - 1
        totalNodes = len(Cu)

        for i in range(n):
            min = 10000

            for A in Cu:
                for B in Cu:
                    if A == B:
                        continue

                    Dab = self._distance(A, B)

                    # Choose closest clusters
                    if Dab <= min:
                        min = Dab

                        savex = A.getValue()[0]
                        savey = B.getValue()[0]

                        # Create new root with clusters as children
                        C = Tree()
                        C.setLeft(A)
                        C.setRight(B)

                        A.setParent(C)
                        B.setParent(C)

                        C.setValue((10000 + i, min, None, None))

            totalNodes += 1

            # Remove closest clusters from Cu and add new cluster
            #print "%d,%d = %f" % (savex, savey, min)
            Cu.remove(C.getLeft())
            Cu.remove(C.getRight())
            Cu.add(C)

        self.tree = Cu.pop()

        self._cluster(self.tree)

    def _distance(self, A, B):

        # If both nodes are leaves, return distance
        if A.getIsLeaf() and B.getIsLeaf():
            return self.dmx[A.getValue()[0]][B.getValue()[0]]

        elif A.getIsLeaf():
            d = self._distance(A, B.getLeft()) + self._distance(A, B.getRight())
            return d / 2.0

        else:
            d = self._distance(A.getRight(), B) + self._distance(A.getLeft(), B)
            return d / 2.0

    def _cluster(self, root):

        if root.getIsLeaf():
            return

        if root.getValue()[1] <= self.minval:
            self.clusters.append(root)
            return

        if root.getLeft():
            self._cluster(root.getLeft())

        if root.getRight():
            self._cluster(root.getRight())
