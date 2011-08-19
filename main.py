#!/usr/bin/python -u

#
# Protocol Informatics Prototype
# Written by Marshall Beddoe <mbeddoe@baselineresearch.net>
# Copyright (c) 2004 Baseline Research
#
# Licensed under the LGPL
#

from PI import *
import sys, getopt

def main():

    print "Protocol Informatics Prototype (v0.01 beta)"
    print "Written by Marshall Beddoe <mbeddoe@baselineresearch.net>"
    print "Copyright (c) 2004 Baseline Research\n"

    print "hello world!"

    # Defaults
    format = None
    weight = 1.0
    graph = False

    #
    # Parse command line options and do sanity checking on arguments
    #
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "pagw:")
    except:
        usage()

    for o,a in opts:
        if o in ["-p"]:
            format = "pcap"
        elif o in ["-a"]:
            format = "ascii"
        elif o in ["-w"]:
            weight = float(a)
        elif o in ["-g"]:
            graph = True
        else:
            usage()

    if len(args) == 0:
        usage()

    if weight < 0.0 or weight > 1.0:
        print "FATAL: Weight must be between 0 and 1"
        sys.exit(-1)

    file = sys.argv[len(sys.argv) - 1]

    try:
        file
    except:
        usage()

    #
    # Open file and get sequences
    #
    if format == "pcap":
        try:
            sequences = input.Pcap(file)
        except:
            print "FATAL: Error opening '%s'" % file
            sys.exit(-1)
    elif format == "ascii":
        try:
            sequences = input.ASCII(file)
        except:
            print "FATAL: Error opening '%s'" % file
            sys.exit(-1)
    else:
        print "FATAL: Specify file format"
        sys.exit(-1)

    if len(sequences) == 0:
        print "FATAL: No sequences found in '%s'" % file
        sys.exit(-1)
    else:
        print "Found %d unique sequences in '%s'" % (len(sequences), file)

    #
    # Create distance matrix (LocalAlignment, PairwiseIdentity, Entropic)
    #
    print "Creating distance matrix ..",
    dmx = distance.LocalAlignment(sequences)
    print "complete"

    #
    # Pass distance matrix to phylogenetic creation function
    #
    print "Creating phylogenetic tree ..",
    phylo = phylogeny.UPGMA(sequences, dmx, minval=weight)
    print "complete"

    #
    # Output some pretty graphs of each cluster
    #
    if graph:
        cnum = 1
        for cluster in phylo:
            out = "graph-%d" % cnum
            print "Creating %s .." % out,
            cluster.graph(out)
            print "complete"
            cnum += 1

    print "\nDiscovered %d clusters using a weight of %.02f" % \
        (len(phylo), weight)

    #
    # Perform progressive multiple alignment against clusters
    #
    i = 1
    alist = []
    for cluster in phylo:
        print "Performing multiple alignment on cluster %d .." % i,
        aligned = multialign.NeedlemanWunsch(cluster)
        print "complete"
        alist.append(aligned)
        i += 1
    print ""

    #
    # Display each cluster of aligned sequences
    #
    i = 1
    for seqs in alist:
        print "Output of cluster %d" % i
        output.Ansi(seqs)
        i += 1
        print ""

def usage():
    print "usage: %s [-gpa] [-w <weight>] <sequence file>" % \
        sys.argv[0]
    print "       -g\toutput graphviz of phylogenetic trees"
    print "       -p\tpcap format"
    print "       -a\tascii format"
    print "       -w\tdifference weight for clustering"
    sys.exit(-1)

if __name__ == "__main__":
    main()
