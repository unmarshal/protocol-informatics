
"""
Distance module

Find distance between sequences

Written by Marshall Beddoe <mbeddoe@baselineresearch.net>
Copyright (c) 2004 Baseline Research

Licensed under the LGPL
"""

#
# Note: Gaps are denoted by the integer value 256 as to avoid '_' problems
#

import align, zlib
from Numeric import *

__all__ = [ "Distance", "Entropic", "PairwiseIdentity", "LocalAlignment" ]

class Distance:

    """Implementation of classify base class"""

    def __init__(self, sequences):
        self.sequences = sequences
        self.N = len(sequences)

        # NxN Distance matrix
        self.dmx = zeros((self.N, self.N), Float)

        for i in range(len(sequences)):
            for j in range(len(sequences)):
                self.dmx[i][j] = -1

        self._go()

    def __repr__(self):
        return "%s" % self.dmx

    def __getitem__(self, i):
        return self.dmx[i]

    def __len__(self):
        return len(self.dmx)

    def _go(self):
        """Perform distance calculations"""
        pass

class Entropic(Distance):

    """Distance calculation based off compression ratios"""

    def _go(self):

        # Similarity matrix
        similar = zeros((self.N, self.N), Float)

        for i in range(self.N):
            for j in range(self.N):
                similar[i][j] = -1

        #
        # Do compression ratio calculations
        #
        for i in range(self.N):
            for j in range(self.N):

                if similar[i][j] >= 0:
                    continue

                seq1 = self.sequences[i][1]
                seq2 = self.sequences[j][1]

                # Convert sequences to strings, gaps denoted by '_'
                seq1str = ""
                for x in seq1:
                    if x == 256:
                        seq1str += '_'
                    else:
                        seq1str += chr(x)

                seq2str = ""
                for x in seq2:
                    if x == 256:
                        seq2str += '_'
                    else:
                        seq2str += chr(x)

                comp1 = zlib.compress(seq1str)
                comp2 = zlib.compress(seq2str)

                if len(comp1) > len(comp2):
                    score = len(comp2) * 1.0 / len(comp1) * 1.0
                else:
                    score = len(comp1) * 1.0 / len(comp2) * 1.0

                similar[i][j] = similar[j][i] = score

        #
        # Distance matrix
        #
        for i in range(self.N):
            for j in range(self.N):
                self.dmx[i][j] = similar[i][i] - similar[i][j]


class PairwiseIdentity(Distance):

    """Distance through basic pairwise similarity"""

    def _go(self):

        # Similarity matrix
        similar = zeros((self.N, self.N), Float)

        for i in range(self.N):
            for j in range(self.N):
                similar[i][j] = -1

        #
        # Find pairs
        #
        for i in range(self.N):
            for j in range(self.N):

                if similar[i][j] >= 0:
                    continue

                seq1 = self.sequences[i][1]
                seq2 = self.sequences[j][1]

                minlen = min(len(seq1), len(seq2))

                len1 = len2 = idents = 0

                for x in range(minlen):
                    if seq1[x] != 256:
                        len1 += 1.0

                        if seq1[x] == seq2[x]:
                            idents += 1.0

                    if seq2[x] != 256:
                        len2 += 1.0

                m = max(len1, len2)

                similar[i][j] = idents / m

        #
        # Distance matrix
        #
        for i in range(self.N):
            for j in range(self.N):
                self.dmx[i][j] = similar[i][i] - similar[i][j]

class LocalAlignment(Distance):

    """Distance through local alignment similarity"""

    def __init__(self, sequences, smx=None):
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

        Distance.__init__(self, sequences)

    def _go(self):

        # Similarity matrix
        similar = zeros((self.N, self.N), Float)

        for i in range(self.N):
            for j in range(self.N):
                similar[i][j] = -1

        #
        # Compute similarity matrix of SW scores
        #
        for i in range(self.N):
            for j in range(self.N):

                if similar[i][j] >= 0:
                    continue

                seq1 = self.sequences[i][1]
                seq2 = self.sequences[j][1]

                (nseq1, nseq2, edits1, edits2, score, gaps) = \
                    align.SmithWaterman(seq1, seq2, self.smx, 0, 0)

                similar[i][j] = similar[j][i] = score

        #
        # Compute distance matrix of SW scores
        #
        for i in range(self.N):
            for j in range(self.N):

                if self.dmx[i][j] >= 0:
                    continue

                self.dmx[i][j] = 1 - (similar[i][j] / similar[i][i])
                self.dmx[j][i] = self.dmx[i][j]
