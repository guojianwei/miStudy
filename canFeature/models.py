# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:02:11 2016

@author: G2
"""


class RNA(object):
    def __init__(self, name, seq):
        self._name = name
        self._seq = seq.replace('t','u').replace('T','U')
            
    def __str__(self):
        # Print the name and first 3 and last 3 letters of the sequence
        return "%s, %s...%s" % (self._name, self._seq[0:3], self._seq[-3:])

    def __eq__(self, y):
        return (self._name == y._name and str(self._seq) == str(y._seq))

    def __hash__(self):
        return hash((self._name, str(self._seq)))
