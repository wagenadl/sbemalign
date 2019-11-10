#!/usr/bin/python3

import copro
import config
import rawimage
from pathlib import Path
import numpy as np
import pyqplot as qp

r = 25
m = 0

s = 100
ifn = rawimage.rawtile(r, m, s)
img = rawimage.loadimage(ifn)
img = rawimage.iscale(img, 5)
y100,x100=np.histogram(img, np.arange(0,22e3,1e3))

s = 150
ifn = rawimage.rawtile(r, m, s)
img = rawimage.loadimage(ifn)
img = rawimage.iscale(img, 5)
y150,x150=np.histogram(img, np.arange(0,22e3,1e3))


qp.figure('/tmp/s1')
qp.brush('555')
qp.bars(x100,y100)
qp.brush('r')
qp.bars(x100,y150,500)

qp.xaxis()
qp.shrink()
