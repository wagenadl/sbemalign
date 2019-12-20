#!/usr/bin/python3

import aligndb
import renderq5utils
import numpy as np
import warp

Q = 5

class Mapper:
    def __init__(self, r, m, s, xywh_model_q5):
        # Following all at scale 1:5
        xl, yt, w, h = xywh_model_q5
        xr = xl + w
        yb = yt + h

        xmdl = np.array([xl, xl, xr, xr])
        ymdl = np.array([yt, yb, yt, yb])
        dx, dy = renderq5utils.interpolatedshifts(r, m, s, xmdl, ymdl)
        x0, y0 = renderq5utils.rigidtileposition(r, m, s)
        xtile = xmdl - x0 - dx
        ytile = ymdl - y0 - dy

        # Following at scale 1:1
        self.xform = warp.getPerspective(xtile*Q, ytile*Q,
                                         xmdl*Q, ymdl*Q)
        self.xtile = xtile * Q # tl, bl, tr, br
        self.ytile = ytile * Q
        self.xl_tile = np.min(self.xtile)
        self.xr_tile = np.max(self.xtile)
        self.yt_tile = np.min(self.ytile)
        self.yb_tile = np.max(self.ytile)

    def __repr__(self):
        return f"[({self.xtile[0]:.0f},{self.ytile[0]:.0f}),({self.xtile[1]:.0f},{self.ytile[1]:.0f}),({self.xtile[2]:.0f},{self.ytile[2]:.0f}),({self.xtile[3]:.0f},{self.ytile[3]:.0f})]"

    def inBBox(self, x, y):
        return x >= self.xl_tile and x <= self.xr_tile \
            and y >= self.yt_tile and y <= self.yb_tile

    def _cross(self, x1, y1, x2, y2):
        return x1*y2 - x2*y1

    def _rightside(self, x1, y1, x2, y2, x, y):
        return self._cross(x-x1, y-y1, x2-x1, y2-y1) <= 0

    def inQuad(self, x, y):
        return self._rightside(self.xtile[0], self.ytile[0],
                               self.xtile[2], self.ytile[2],
                               x, y) \
               and \
               self._rightside(self.xtile[2], self.ytile[2],
                               self.xtile[3], self.ytile[3],
                               x, y) \
               and \
               self._rightside(self.xtile[3], self.ytile[3],
                               self.xtile[1], self.ytile[1],
                               x, y) \
               and \
               self._rightside(self.xtile[1], self.ytile[1],
                               self.xtile[0], self.ytile[0],
                               x, y)
    def map(self, x, y):
        return warp.applyPerspective(self.xform, x, y)

moma = {}
    
def montagemapper(r, m, s):
    '''MONTAGEMAPPER - Returns a list of mappers for all subtiles of montage
    mp = MONTAGEMAPPER(r, m, s) returns a list containing Mappers
    for all the subtiles in the tile (r, m, s).
    The result is cached.'''
    rms = (r,m,s)
    if rms in moma:
        return moma[rms]
    
    mp = []
    xx, yy = renderq5utils.rendergrid(r, m, s)
    IX = len(xx) - 1
    IY = len(yy) - 1
    for ix in range(IX):
        for iy in range(IY):
            xl1 = xx[ix]
            xr1 = xx[ix+1]
            yt1 = yy[iy]
            yb1 = yy[iy+1]
            if xr1 <= xl1 or yb1 <= yt1:
                continue
            xywh = [xl1, yt1, xr1-xl1, yb1-yt1]
            mp.append(Mapper(r, m, s, xywh))

    moma[rms] = mp
    return mp
