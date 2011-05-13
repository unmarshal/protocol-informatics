
"""
Binary Tree module

Implementation of binary tree class.

Written by Marshall Beddoe <mbeddoe@baselineresearch.net>
Copyright (c) 2004 Baseline Research

Licensed under the LGPL
"""

from pydot import *

class Tree:

    def __init__(self):

        self._parentNode = None
        self._leftNode = None
        self._rightNode = None
        self._value = None
        self._height = 0

    def getIsLeaf(self):
        """Is this tree a leaf node"""

        if self._leftNode == None and self._rightNode == None:
            return True
        else:
            return False

    def getHeight(self):
        """Return height at this root"""
        return self._height

    def getParent(self):
        """Return parent node"""
        return self._parentNode

    def setParent(self, parentNode):
        """Set new parent node"""
        self._parentNode = parentNode

    def getLeft(self):
        """Return left child node"""
        return self._leftNode

    def setLeft(self, leftNode):
        """Set left child node"""
        self._leftNode = leftNode

    def getRight(self):
        """Return right child node"""
        return self._rightNode

    def setRight(self, rightNode):
        """Set right child node"""
        self._rightNode = rightNode

    def getValue(self):
        """Return node value"""
        return self._value

    def setValue(self, value):
        """Set node value"""
        self._value = value

    def graph(self, output, format="raw"):
        """Graph tree from root using graphviz"""

        self.i = 0
        self.graph = Dot(center="TRUE",rankdir="TB")
        #self.graph = Dot(size="3.5,4",page="4.5,6",center="TRUE",rankdir="TB")
        self.subgraph = Subgraph("subG", rank="same")

        self._traverse(self)

        self.graph.add_subgraph(self.subgraph)

        #if format == "raw":
        self.graph.write_raw(output + ".dot", prog="dot")
        self.graph.write_png(output + ".png", prog="dot")
        #elif format == "png":
        #    self.graph.write_png(output + ".png", prog="dot")
        #else:
        #    raise "UnknownFormat"

    def _traverse(self, root):

        if root.getParent():

            if root.getIsLeaf():
                v1 = root.getValue()[2][0]
            else:
                v1 = root.getValue()[0]

            weight = root.getValue()[1]

            v2 = root.getParent().getValue()[0]

            l = "%.02f%%" % (weight * 100.0)

            if v1 >= 10000:
                node1 = Node(v1, shape="plaintext", ratio="auto", label=l)
            else:
                node1 = Node(v1, shape="house", ratio="auto")
                node1.set("style", "filled")
                node1.set("fillcolor", "cyan")

            if v2 >= 10000:
                weight = root.getParent().getValue()[1]
                l = "%.02f%%" % (weight * 100.0)
                node2 = Node(v2, shape="plaintext", ratio="auto", label=l)

            if root.getIsLeaf():
                self.subgraph.add_node(node1)

            self.graph.add_node(node1)
            self.graph.add_node(node2)

            edge = Edge(v2, v1)

            self.graph.add_edge(edge)

        if root.getLeft():
            self._traverse(root.getLeft())

        if root.getRight():
            self._traverse(root.getRight())
