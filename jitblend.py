#!/usr/bin/python3

'''This is an attempt to construct imagery for an arbitrary rectangle
in model space based on (scaled) source tiles.'''

import aligndb
import em170428.runinfo
import cv2
import numpy as np

ri = em170428.runinfo.RunInfo()

def involvedtiles(z, xywh, stage):
    '''INVOLVEDTILES - Returns a list of (r,m,s) triplets needed to fill rect.
    lst = INVOLVEDTILES(z, xywh, stage) returns a list of (r,m,s,t) tuples,
    where (r,m,s) specifies the tiles needed to produce imagery at given Z for
    the rectangle specified by XYWH in global model coordinates, at given
    STAGE, and T is the global-to-tile transformation.
    Stages are as for ALIGNDB.GLOBALTOTILE.
    '''
    r,s = ri.zToRunSlice(z)
    x0, y0, w, h = xywh
    M = ri.montageCount(r)
    rms = []
    for m in range(M):
        t = aligndb.globaltotile(r, m, s, stage)
        p_tl = swiftir.applyAffine(t, [x0, y0])
        p_tr = swiftir.applyAffine(t, [x0 + w, y0])
        p_bl = swiftir.applyAffine(t, [x0, y0 + h])
        p_br = swiftir.applyAffine(t, [x0 + w, y0 + h])
        # Does this intersect with [0,0]-[aligndb.X,aligndb.Y] ?
        xmin = np.min((p_tl[0], p_bl[0]))
        xmax = np.max((p_tr[0], p_br[0]))
        ymin = np.min((p_tl[1], p_tr[1]))
        ymax = np.max((p_bl[1], p_br[1]))
        if xmin<aligndb.X and xmax>0:
            if ymin<aligndb.Y and ymax>0:
                mms.append((r,m,s,t))
    return rms

def matmul(lst):
    a = None
    for b in lst:
        if a is None:
            a = b.copy()
        else:
            a = np.dot(a,b,a)
    return a

def noblend(xywh, rms, q, t):
    '''NOBLEND - Extract ROI from transformed image
    img = NOBLEND(xywh, r,m,s,q, t) extracts the rectangle XYWH, specified
    in global model coordinates, from the tile specified by R,M,S,
    given the global-to-tile transform T, at scale Q. Note that whereas
    both the source tile and the final image is produced at scale Q, the
    transform T and rectangle XYWH are given in unscaled coordinates.'''
    t0 = np.eye(3)
    x0,y0,w,h = xywh
    r, m, s = rms
    t0[0,2] = x0
    t0[1,2] = y0
    t1 = np.eye(3)
    t1[0][0] = q
    t1[1][1] = q
    t2 = np.eye(3)
    t2[0][0] = 1./q
    t2[1][1] = 1./q
    img = rawimage.loadimage(rawimage.scaledtile(r,m,s,q))
    img = img.astype(np.float32)
    img[img==0] = 1 # So we can do transparency without going crazy
    t = matmul([t2, t, t0, t1])
    t = t[0:2,0:3]
    img = cv2.warpAffine(img, t, (w//q,h//q),
                         flags=cv2.WARP_INVERSE_MAP + cv2.INTER_LINEAR)
    
def blend(z, q, xywh, stage, mexcl=[]):
    '''BLEND - Return an image based on (scaled) source tiles.
    img = BLEND(z, q, xywh, stage) creates imagery at given Z for the
    global-coordinate rectangle XYWH, at given scale Q and given
    processing STAGE (see ALIGNDB.GLOBALTOTILE), using all the tiles
    that overlap with that rectangle. Pixels with no data are returned
    as 0; pixels with data are guaranteed to not be 0.
    Optional argument MEXCL can be used to specify a list of montages
    to be excluded.'''
    
    tilelist = involvedtiles(z, xywh, stage)
    simg = None
    smsk = None
    n = 0
    for spec in tilelist:
        r,m,s,t = spec
        if m not in mexcl:
            n += 1
    for spec in tilelist:
        r,m,s,t = spec
        if m in mexcl:
            continue
        img = noblend(xywh, (r,m,s), q, t)
        if simg is None:
            simg = img
        else:
            simg += img
        if n>1: # We'll be blending
            msk = (img>=1).astype(img.dtype)
            if smsk is None:
                smsk = msk + 1e-20
            else:
                smsk += msk
    if n>1:
        return simg ./ smsk
    else:
        return simg
