#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

import optimizing

r = 1

print('Gathering deltas')
deltas = optimizing.AllDeltas(r,
                              crosstbl='slicealignq5',
                              intratbl='montagealignq5relhp',
                              edgetbl='montagealignattouchq5relhp')
print('Calculating overall montage positions')
deltas.makemontpos()

print('Creating index')
idx = optimizing.Index(deltas)
print('Creating matrix')
matx = optimizing.Matrix(deltas, idx, 'x', w_cross=100)
maty = optimizing.Matrix(deltas, idx, 'y', w_cross=100)
print('Solving matrix')
soln = optimizing.Solution(deltas, matx, maty)
print('Collecting results')
usol = optimizing.UnifiedSolution(soln)
