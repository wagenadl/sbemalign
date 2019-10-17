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

deltas = optimizing.AllDeltas(r,
                              crosstbl='slicealignq5',
                              intratbl='montagealignq5relhp',
                              edgetbl='montagealignattouchq5relhp')

idx = optimizing.Index(deltas)
mat = optimizing.Matrix(deltas, idx, w_cross=100)
soln = optimizing.Solution(deltas, mat)
usol = optimizing.UnifiedSolution(soln)
