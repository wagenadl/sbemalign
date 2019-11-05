#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import pyqplot as qp
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

db = aligndb.DB()
ri = db.runinfo()

X = Y = 684 # This is correct for our Q5 work
MAXX = X*5  # This is correct for our Q5 work
