#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warp
import pyqplot as qp
import numpy as np

db = aligndb.DB()
ri = db.runinfo()

snrintra = db.vsel('select snrb from montagealignq5relhp')
snredge =  db.vsel('select snrb from montagealignattouchq5relhp where (ix>0 and ix<5) or (iy>0 and iy<5)')
snrcorner =  db.vsel('select snrb from montagealignattouchq5relhp where (ix=0 or ix=5) and (iy=0 or iy=5)')
snrcross =  db.vsel('select snrc from slicealignq5 where (ix1>0 and ix1<5) or (iy1>0 and iy1<5)')

qp.figure('/tmp/snr')
qp.subplot(2,2,1)
yy,xx = np.histogram(snrintra[0], np.arange(100))
xx=(xx[:-1] + xx[1:])/2
qp.brush('777')
qp.bars(xx,yy,1)
qp.xaxis('snr intra', np.arange(0,100,10))
qp.yaxis('counts')
qp.shrink();

qp.subplot(2,2,2)
yy,xx = np.histogram(snredge[0], np.arange(0,50,.25))
xx=(xx[:-1] + xx[1:])/2
qp.brush('777')
qp.bars(xx,yy,.25)
qp.xaxis('snr edge', np.arange(0,50,10))
qp.yaxis('counts')
qp.shrink();


qp.subplot(2,2,3)
yy,xx = np.histogram(snrcross[0], np.arange(100))
xx=(xx[:-1] + xx[1:])/2
qp.brush('777')
qp.bars(xx,yy,1)
qp.xaxis('snr cross', np.arange(0,100,10))
qp.yaxis('counts')
qp.shrink();


qp.subplot(2,2,4)
yy,xx = np.histogram(snrcorner[0], np.arange(0,50,.25))
xx=(xx[:-1] + xx[1:])/2
qp.brush('777')
qp.bars(xx,yy,.25)
qp.xaxis('snr corner', np.arange(0,50,10))
qp.yaxis('counts')
qp.shrink();



