#!/usr/bin/python3

import copro
import config
import rawimage
from pathlib import Path

def produceq25(r, m, s, ofn=None):
    ifn = rawimage.rawtile(r, m, s)
    img = rawimage.loadimage(ifn)
    img = rawimage.iscale(img, 5)
    img = rawimage.to8bit(img, .1)

    # Create 5x5 subtiles at 1/5th
    Y,X = img.shape
    R = Y//5
    ofntmpl = Path(rawimage.scaledtile(r, m, s, 5))
    for x in range(5):
        for y in range(5):
            imti = img[R*y:R*(y+1), R*x:R*(x+1)]
            ofnti = ofntmpl.with_suffix('.%i%i.tif' % (x, y))
            rawimage.saveimage(imti, ofnti)
    Path(ofntmpl.with_suffix('.complete')).touch()

    img = rawimage.iscale(img, 5)
    if ofn is None:
        ofn = rawimage.scaledtile(r, m, s, 25)
    rawimage.saveimage(img, ofn)

class ProduceQ25(copro.CoPro):
    def filename(self, pars):
        r = pars[0]
        m = pars[1]
        s = pars[2]
        return rawimage.scaledtile(r, m, s, 25)
    def produce(self, pars, ofn):
        r = pars[0]
        m = pars[1]
        s = pars[2]
        produceq25(r, m, s, ofn)

if __name__ == '__main__':
    import sys
    r = int(sys.argv[1])
    m = int(sys.argv[2])
    s = int(sys.argv[3])
    c = ProduceQ25()
    c.ensure((r, m, s))
    
    
