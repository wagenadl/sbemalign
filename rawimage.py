#!/usr/bin/python3

import numpy as np
import cv2
import runinfo
import urllib.request
import re

def loadraw(r, m, s):
    '''LOADRAW - Load raw image
    img = LOADRAW(r, m, s) loads the raw image for given run/montage/slice
    as a 16-bit 17100x17100-pixel image.'''
    root = runinfo.rawroot
    ifn = '%s/Run%i/Montage_%03i/Run%i_OnPoint_%04i.tif'
    ifn = ifn % (root, r, m, r, s)
    img = cv2.imread(ifn, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)
    return img

def loadimage(url):
    '''LOADIMAGE - Download an image from anywhere
    img = LOADIMAGE(url) downloads an image from anywhere.
    Result is grayscale. Bit depth is as from source.
    For convenience, URL may also be a local file name (even without
    file:// prefix).
    For unknown reasons, opencv's imdecode is slower than imread, so
    file:// is much faster than url:// even if the file system is a
    remote sshfs mount.'''
    if url.startswith('file://'):
        url = url[7:]
    if re.compile('^[a-z]+://').match(url):
        resp = urllib.request.Request(url=url)
        with urllib.request.urlopen(resp) as f:
            dat = np.fromfile(f, dtype=np.uint8)
        return cv2.imdecode(dat, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)
    else:
        return cv2.imread(url, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)

def ihist(img):
    '''IHIST - Calculate image histogram
    hst = IHIST(img) where IMG is either an 8-bit or a 16-bit image, calculates
    the image histogram. The result is a vector with 256 or 65536 entries.
    Histogram calculation is surprisingly slow for large images, despite
    the use of opencv.'''
    if img.dtype==np.uint8:
        hst = cv2.calcHist([img], [0], None, [256], [0, 256])
        return hst[:,0]
    elif img.dtype==np.uint16:
        hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
        return hst[:,0]
    else:
        raise ValueError('Only 8 or 16 bit images are supported')

def to8bit(img, stretch=.1):
    '''TO8BIT - Convert 16-bit image to 8 bits
    im8 = TO8BIT(im16) converts a 16-bit image to 8-bit format, clipping
    extreme values.
    Optional argument STRETCH determines what percentage of pixels gets
    clipped. Default is 0.1%.
    Calling TO8BIT on an 8-bit image has no effect; the contrast is not
    stretched in that case.'''
    if img.dtype==np.uint8:
        return img
    elif img.dtype==np.uint16:
        hst = cv2.calcHist([img], [0], None, [65536], [0, 65536])
        hst = np.cumsum(hst)
        hst = hst / hst[-1]
        lowb = np.argmax(hst>=.01*stretch)
        upb = np.argmax(hst>=1-.01*stretch)

        def mklut(lowb, upb):
            lut = np.arange(65536)
            lut.clip(lowb, upb, out=lut)
            lut -= lowb
            np.floor_divide(lut, (upb-lowb+1) / 256, out=lut, casting='unsafe')
            return lut.astype(np.uint8)
        return mklut(lowb, upb)[img]
    else:
        raise ValueError('Only 8 or 16 bit images are supported')

def iscale(img, n):
    '''ISCALE - Downscale an image
    im = ISCALE(img, n) downscales the image IMG by a factor N in both
    x- and y-directions using "area" interpolation.'''
    y,x = img.shape
    return cv2.resize(img, (y//n,x//n), interpolation=cv2.INTER_AREA)



