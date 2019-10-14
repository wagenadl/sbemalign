#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import swiftir
import numpy as np

import rawimage

from PyQt5.QtWidgets import QApplication, QLabel
from PyQt5.QtGui import QPixmap, QImage                                

r = int(sys.argv[1])
m = int(sys.argv[2])
ix = int(sys.argv[3])
iy = int(sys.argv[4])
tbl = 'montagealignq5a'
if len(sys.argv)>=6:
    tbl = tbl[:-1] + sys.argv[5]

print(f'Working on showing R{r} M{m} ix{ix} iy{iy} from {tbl}')
db = aligndb.DB()
ri = db.runinfo()
S = ri.nslices(r)

app = QApplication([])
lbl = QLabel('hello world')
lbl.show()

for s in range(S):
    img = rawimage.partialq5img(r,m,s,ix,iy)
    Y,X = img.shape
    if tbl.endswith('0'):
        dx = 0
        dy = 0
        snrb = 0
    else:
        dx,dy,snrb = db.sel(f'''select dx+dxb,dy+dyb,snrb from {tbl}
        where r={r} and m={m} and s={s} and ix={ix} and iy={iy}''')[0]
        afm = np.array([[1,0,dx],[0,1,dy]])
        img = swiftir.cv2.warpAffine(img, afm, (Y,X),
                                     flags=swiftir.cv2.INTER_LINEAR,
                                     borderMode=swiftir.cv2.BORDER_CONSTANT,
                                     borderValue=np.mean(img))
    img = img.astype(np.uint8)
    img = QImage(img, X,Y, QImage.Format_Grayscale8)
    lbl.setPixmap(QPixmap(img))
    lbl.resize(X,Y)
    lbl.setWindowTitle('S%i %5.1f %5.1f [%5.1f]' % (s,dx,dy,snrb))
    app.processEvents()
